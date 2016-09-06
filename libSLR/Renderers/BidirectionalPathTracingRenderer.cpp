//
//  BidirectionalPathTracingRenderer.cpp
//
//  Created by 渡部 心 on 2016/01/31.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "BidirectionalPathTracingRenderer.h"

#include "../Core/RenderSettings.h"
#include "../Helper/ThreadPool.h"
#include "../RNGs/XORShiftRNG.h"
#include "../Memory/ArenaAllocator.h"
#include "../Core/ImageSensor.h"
#include "../Core/RandomNumberGenerator.h"
#include "../Core/cameras.h"
#include "../Core/SurfaceObject.h"
#include "../Core/light_path_samplers.h"
#include "../Core/distributions.h"

namespace SLR {
    BidirectionalPathTracingRenderer::BidirectionalPathTracingRenderer(uint32_t spp) : m_samplesPerPixel(spp) {
    }
    
    void BidirectionalPathTracingRenderer::render(const Scene &scene, const RenderSettings &settings) const {
#ifdef DEBUG
        uint32_t numThreads = 1;
#else
        uint32_t numThreads = std::thread::hardware_concurrency();
#endif
        XORShiftRNG topRand(settings.getInt(RenderSettingItem::RNGSeed));
        std::unique_ptr<ArenaAllocator[]> mems = std::unique_ptr<ArenaAllocator[]>(new ArenaAllocator[numThreads]);
        std::unique_ptr<IndependentLightPathSampler[]> samplers = std::unique_ptr<IndependentLightPathSampler[]>(new IndependentLightPathSampler[numThreads]);
        for (int i = 0; i < numThreads; ++i) {
            new (mems.get() + i) ArenaAllocator();
            new (samplers.get() + i) IndependentLightPathSampler(topRand.getUInt());
        }
        std::unique_ptr<IndependentLightPathSampler*[]> samplerRefs = std::unique_ptr<IndependentLightPathSampler*[]>(new IndependentLightPathSampler*[numThreads]);
        for (int i = 0; i < numThreads; ++i)
            samplerRefs[i] = &samplers[i];
        
        const Camera* camera = scene.getCamera();
        ImageSensor* sensor = camera->getSensor();
        
        Job job;
        job.scene = &scene;
        
        job.mems = mems.get();
        job.pathSamplers = samplerRefs.get();
        
        job.camera = camera;
        job.timeStart = settings.getFloat(RenderSettingItem::TimeStart);
        job.timeEnd = settings.getFloat(RenderSettingItem::TimeEnd);
        
        job.sensor = sensor;
        job.imageWidth = settings.getInt(RenderSettingItem::ImageWidth);
        job.imageHeight = settings.getInt(RenderSettingItem::ImageHeight);
        job.numPixelX = sensor->tileWidth();
        job.numPixelY = sensor->tileHeight();
        
        uint32_t exportPass = 1;
        uint32_t imgIdx = 0;
        uint32_t endIdx = 16;
        
        sensor->init(job.imageWidth, job.imageHeight);
        sensor->addSeparatedBuffers(numThreads);
        
        std::chrono::system_clock::time_point start, end;
        start = std::chrono::system_clock::now();
        
        for (int s = 0; s < m_samplesPerPixel; ++s) {
            ThreadPool threadPool(numThreads);
            for (int ty = 0; ty < sensor->numTileY(); ++ty) {
                for (int tx = 0; tx < sensor->numTileX(); ++tx) {
                    job.basePixelX = tx * sensor->tileWidth();
                    job.basePixelY = ty * sensor->tileHeight();
                    threadPool.enqueue(std::bind(&Job::kernel, job, std::placeholders::_1));
                }
            }
            threadPool.wait();
            
            if ((s + 1) == exportPass) {
                char filename[256];
                sprintf(filename, "%03u.bmp", imgIdx);
                end = std::chrono::system_clock::now();
                double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
                sensor->saveImage(filename, settings.getFloat(RenderSettingItem::Brightness) / (s + 1));
                printf("%u samples: %s, %g[s]\n", exportPass, filename, elapsed * 0.001f);
                ++imgIdx;
                if (imgIdx == endIdx)
                    break;
                exportPass += exportPass;
            }
        }
        
        //    sensor.saveImage("output.png", settings.getFloat(RenderSettingItem::SensorResponse) / numSamples);
    }
    
    void BidirectionalPathTracingRenderer::Job::kernel(uint32_t threadID) {
        ArenaAllocator &mem = mems[threadID];
        IndependentLightPathSampler &pathSampler = *pathSamplers[threadID];
        for (int ly = 0; ly < numPixelY; ++ly) {
            for (int lx = 0; lx < numPixelX; ++lx) {
                float time = pathSampler.getTimeSample(timeStart, timeEnd);
                PixelPosition p = pathSampler.getPixelPositionSample(basePixelX + lx, basePixelY + ly);
                
                float selectWLPDF;
                WavelengthSamples wls = WavelengthSamples::createWithEqualOffsets(pathSampler.getWavelengthSample(), pathSampler.getWLSelectionSample(), &selectWLPDF);
                
                // initialize working area for the current pixel.
                curPx = p.x;
                curPy = p.y;
                wlHint = wls.selectedLambda;
                eyeVertices.clear();
                lightVertices.clear();
                
                // light subpath generation
                {
                    // select one light from all the lights in the scene.
                    float lightProb;
                    Light light;
                    scene->selectLight(pathSampler.getLightSelectionSample(), &light, &lightProb);
                    SLRAssert(!std::isnan(lightProb) && !std::isinf(lightProb), "lightProb: unexpected value detected: %f", lightProb);
                    
                    // sample a ray with its radiance (emittance, EDF value) from the selected light.
                    LightPosQuery lightPosQuery(time, wls);
                    LightPosQueryResult lightPosResult;
                    EDFQuery edfQuery;
                    EDFQueryResult edfResult;
                    EDF* edf;
                    SampledSpectrum Le0, Le1;
                    Ray ray = light.sampleRay(lightPosQuery, pathSampler.getLightPosSample(), &lightPosResult, &Le0, &edf,
                                              edfQuery, pathSampler.getEDFSample(), &edfResult, &Le1, mem);
                    
                    // register the first light vertex.
                    float lightAreaPDF = lightProb * lightPosResult.areaPDF;
                    lightVertices.emplace_back(lightPosResult.surfPt, Vector3D::Zero, Normal3D(0, 0, 1), mem.create<EDFProxy>(edf),
                                               Le0 / lightAreaPDF, lightAreaPDF, 1.0f, lightPosResult.posType, WavelengthSamples::Flag(0));
                    
                    // create subsequent light subpath vertices by tracing in the scene.
                    SampledSpectrum alpha = lightVertices.back().alpha * Le1 * (absDot(ray.dir, lightPosResult.surfPt.gNormal) / edfResult.dirPDF);
                    generateSubPath(wls, alpha, ray, edfResult.dirPDF, edfResult.dirType, edfResult.dir_sn.z, true, pathSampler, mem);
                }
                
                // eye subpath generation
                {
                    // sample a ray with its importances (spatial, directional) from the lens and its IDF.
                    LensPosQuery lensQuery(time, wls);
                    LensPosQueryResult lensResult;
                    IDFSample WeSample(p.x / imageWidth, p.y / imageHeight);
                    IDFQueryResult WeResult;
                    IDF* idf;
                    SampledSpectrum We0, We1;
                    Ray ray = camera->sampleRay(lensQuery, pathSampler.getLensPosSample(), &lensResult, &We0, &idf, WeSample, &WeResult, &We1, mem);
                    
                    // register the first eye vertex.
                    eyeVertices.emplace_back(lensResult.surfPt, Vector3D::Zero, Normal3D(0, 0, 1), mem.create<IDFProxy>(idf),
                                             We0 / (lensResult.areaPDF * selectWLPDF), lensResult.areaPDF, 1.0f, lensResult.posType, WavelengthSamples::Flag(0));
                    
                    // create subsequent eye subpath vertices by tracing in the scene.
                    SampledSpectrum alpha = eyeVertices.back().alpha * We1 * (absDot(ray.dir, lensResult.surfPt.gNormal) / WeResult.dirPDF);
                    generateSubPath(wls, alpha, ray, WeResult.dirPDF, WeResult.dirType, WeResult.dirLocal.z, false, pathSampler, mem);
                }
                
                // connection
                for (int t = 1; t <= eyeVertices.size(); ++t) {
                    const BPTVertex &eVtx = eyeVertices[t - 1];
                    for (int s = 1; s <= lightVertices.size(); ++s) {
                        const BPTVertex &lVtx = lightVertices[s - 1];
                        
                        // calculate the remaining factors of the full path
                        // that are not included in the precomputed weights.
                        // ----------------------------------------------------------------
                        float connectDist2;
                        Vector3D connectionVector = lVtx.surfPt.getDirectionFrom(eVtx.surfPt.p, &connectDist2);
                        float cosLightEnd = absDot(connectionVector, lVtx.surfPt.gNormal);
                        float cosEyeEnd = absDot(connectionVector, eVtx.surfPt.gNormal);
                        float G = cosEyeEnd * cosLightEnd / connectDist2;
                        
                        Vector3D lConnectVector = lVtx.surfPt.shadingFrame.toLocal(-connectionVector);
                        DDFQuery queryLightEnd{lVtx.dirIn_sn, lVtx.gNormal_sn, wlHint, true};
                        SampledSpectrum lRevDDF;
                        SampledSpectrum lDDF = lVtx.ddf->evaluate(queryLightEnd, lConnectVector, &lRevDDF);
                        float eExtend2ndDirPDF;
                        float lExtend1stDirPDF = lVtx.ddf->evaluatePDF(queryLightEnd, lConnectVector, &eExtend2ndDirPDF);
                        
                        Vector3D eConnectVector = eVtx.surfPt.shadingFrame.toLocal(connectionVector);
                        DDFQuery queryEyeEnd{eVtx.dirIn_sn, eVtx.gNormal_sn, wlHint, false};
                        SampledSpectrum eRevDDF;
                        SampledSpectrum eDDF = eVtx.ddf->evaluate(queryEyeEnd, eConnectVector, &eRevDDF);
                        float lExtend2ndDirPDF;
                        float eExtend1stDirPDF = eVtx.ddf->evaluatePDF(queryEyeEnd, eConnectVector, &lExtend2ndDirPDF);
                        
                        float wlProb = 1.0f;
                        if ((lVtx.wlFlags | eVtx.wlFlags) & WavelengthSamples::LambdaIsSelected)
                            wlProb = 1.0f / WavelengthSamples::NumComponents;
                        SampledSpectrum connectionTerm = lDDF * (G / wlProb) * eDDF;
                        if (connectionTerm == SampledSpectrum::Zero)
                            continue;
                        
                        if (!scene->testVisibility(eVtx.surfPt, lVtx.surfPt, time))
                            continue;
                        // ----------------------------------------------------------------
                        
                        // calculate the 1st and 2nd subpath extending PDFs and probabilities.
                        // They can't be stored in advance because they depend on the connection.
                        float lExtend1stAreaPDF, lExtend1stRRProb, lExtend2ndAreaPDF, lExtend2ndRRProb;
                        {
                            lExtend1stAreaPDF = lExtend1stDirPDF * cosEyeEnd / connectDist2;
                            lExtend1stRRProb = s > 1 ? std::min((lDDF * cosLightEnd / lExtend1stDirPDF).importance(wlHint), 1.0f) : 1.0f;
                            
                            if (t > 1) {
                                BPTVertex &eVtxNextToEnd = eyeVertices[t - 2];
                                float dist2;
                                Vector3D dir2nd = eVtx.surfPt.getDirectionFrom(eVtxNextToEnd.surfPt.p, &dist2);
                                lExtend2ndAreaPDF = lExtend2ndDirPDF * absDot(eVtxNextToEnd.surfPt.gNormal, dir2nd) / dist2;
                                lExtend2ndRRProb = std::min((eRevDDF * absDot(eVtx.gNormal_sn, eVtx.dirIn_sn) / lExtend2ndDirPDF).importance(wlHint), 1.0f);
                            }
                        }
                        float eExtend1stAreaPDF, eExtend1stRRProb, eExtend2ndAreaPDF, eExtend2ndRRProb;
                        {
                            eExtend1stAreaPDF = eExtend1stDirPDF * cosLightEnd / connectDist2;
                            eExtend1stRRProb = t > 1 ? std::min((eDDF * cosEyeEnd / eExtend1stDirPDF).importance(wlHint), 1.0f) : 1.0f;
                            
                            if (s > 1) {
                                BPTVertex &lVtxNextToEnd = lightVertices[s - 2];
                                float dist2;
                                Vector3D dir2nd = lVtxNextToEnd.surfPt.getDirectionFrom(lVtx.surfPt.p, &dist2);
                                eExtend2ndAreaPDF = eExtend2ndDirPDF * absDot(lVtxNextToEnd.surfPt.gNormal, dir2nd) / dist2;
                                eExtend2ndRRProb = std::min((lRevDDF * absDot(lVtx.gNormal_sn, lVtx.dirIn_sn) / eExtend2ndDirPDF).importance(wlHint), 1.0f);
                            }
                        }
                        
                        // calculate MIS weight and store weighted contribution to a sensor.
                        float MISWeight = calculateMISWeight(lExtend1stAreaPDF, lExtend1stRRProb, lExtend2ndAreaPDF, lExtend2ndRRProb,
                                                             eExtend1stAreaPDF, eExtend1stRRProb, eExtend2ndAreaPDF, eExtend2ndRRProb, s, t);
                        if (std::isinf(MISWeight) || std::isnan(MISWeight))
                            continue;
                        SLRAssert(MISWeight >= 0 && MISWeight <= 1.0f, "invalid MIS weight: %g", MISWeight);
                        SampledSpectrum contribution = MISWeight * lVtx.alpha * connectionTerm * eVtx.alpha;
                        SLRAssert(!contribution.hasNaN() && !contribution.hasInf() && !contribution.hasMinus(),
                                  "Unexpected value detected: %s\n"
                                  "pix: (%f, %f)", contribution.toString().c_str(), px, py);
                        if (t > 1) {
                            sensor->add(p.x, p.y, wls, contribution);
                        }
                        else {
                            const IDF* idf = (const IDF*)eVtx.ddf->getDDF();
                            float hitPx, hitPy;
                            idf->calculatePixel(eConnectVector, &hitPx, &hitPy);
                            sensor->add(threadID, hitPx, hitPy, wls, contribution);
                        }
                    }
                }
                
                mem.reset();
            }
        }
    }
    
    void BidirectionalPathTracingRenderer::Job::generateSubPath(const WavelengthSamples &initWLs, const SampledSpectrum &initAlpha, const SLR::Ray &initRay, float dirPDF, DirectionType sampledType,
                                                                float cosLast, bool adjoint, IndependentLightPathSampler &pathSampler, SLR::ArenaAllocator &mem) {
        std::vector<BPTVertex> &vertices = adjoint ? lightVertices : eyeVertices;
        
        // reject invalid values.
        if (dirPDF == 0.0f)
            return;
        
        WavelengthSamples wls = initWLs;
        Ray ray = initRay;
        SampledSpectrum alpha = initAlpha;
        
        Intersection isect;
        SurfacePoint surfPt;
        float RRProb = 1.0f;
        while (scene->intersect(ray, &isect)) {
            isect.getSurfacePoint(&surfPt);
            float dist2 = squaredDistance(vertices.back().surfPt, surfPt);
            Vector3D dirOut_sn = surfPt.shadingFrame.toLocal(-ray.dir);
            Normal3D gNorm_sn = surfPt.shadingFrame.toLocal(surfPt.gNormal);
            BSDF* bsdf = surfPt.createBSDF(wls, mem);
            
            float areaPDF = dirPDF * absDot(dirOut_sn, gNorm_sn) / dist2;
            vertices.emplace_back(surfPt, dirOut_sn, gNorm_sn, mem.create<BSDFProxy>(bsdf), alpha, areaPDF, RRProb, sampledType, wls.flags);
            
            // implicit path (zero light subpath vertices, s = 0)
            if (!adjoint && surfPt.isEmitting()) {
                EDF* edf = surfPt.createEDF(wls, mem);
                SampledSpectrum Le0 = surfPt.emittance(wls);
                SampledSpectrum Le1 = edf->evaluate(EDFQuery(), dirOut_sn);
                
                float lightProb = scene->evaluateProb(Light(isect.obj));
                float extend1stAreaPDF = lightProb * surfPt.evaluateAreaPDF();
                float extend2ndAreaPDF = edf->evaluatePDF(EDFQuery(), dirOut_sn) * cosLast / dist2;
                
                float MISWeight = calculateMISWeight(extend1stAreaPDF, 1.0f, extend2ndAreaPDF, 1.0f,
                                                     0.0f, 0.0f, 0.0f, 0.0f,
                                                     0, (uint32_t)vertices.size());
                if (!std::isinf(MISWeight) && !std::isnan(MISWeight)) {
                    SampledSpectrum contribution = MISWeight * alpha * Le0 * Le1;
                    SLRAssert(MISWeight >= 0 && MISWeight <= 1.0f, "invalid MIS weight: %g", MISWeight);
                    SLRAssert(!contribution.hasNaN() && !contribution.hasInf() && !contribution.hasMinus(),
                              "Unexpected value detected: %s\n"
                              "pix: (%f, %f)", contribution.toString().c_str(), curPx, curPy);
                    if (wls.flags & WavelengthSamples::LambdaIsSelected)
                        contribution *= WavelengthSamples::NumComponents;
                    sensor->add(curPx, curPy, wls, contribution);
                }
            }
            
            if (surfPt.atInfinity) {
                vertices.pop_back();
                break;
            }
            
            BSDFQuery fsQuery(dirOut_sn, gNorm_sn, wls.selectedLambda, DirectionType::All, adjoint);
            BSDFQueryResult fsResult;
            BSDFReverseInfo revInfo;
            fsResult.reverse = &revInfo;
            SampledSpectrum fs = bsdf->sample(fsQuery, pathSampler.getBSDFSample(), &fsResult);
            if (fs == SampledSpectrum::Zero || fsResult.dirPDF == 0.0f)
                break;
            if (fsResult.dirType.isDispersive())
                wls.flags |= WavelengthSamples::LambdaIsSelected;
            float cosIn = absDot(fsResult.dir_sn, gNorm_sn);
            SampledSpectrum weight = fs * (cosIn / fsResult.dirPDF);
            
            // Russian roulette
            RRProb = std::min(weight.importance(wlHint), 1.0f);
            if (pathSampler.getPathTerminationSample() < RRProb)
                weight /= RRProb;
            else
                break;
            
            alpha *= weight;
            ray = Ray(surfPt.p, surfPt.shadingFrame.fromLocal(fsResult.dir_sn), ray.time, Ray::Epsilon);
            SLRAssert(!weight.hasInf() && !weight.hasNaN(),
                      "weight: unexpected value detected:\nweight: %s\nfs: %s\nlength: %u, cos: %g, dirPDF: %g",
                      weight.toString().c_str(), fs.toString().c_str(), uint32_t(vertices.size()) - 1, cosIn, fsResult.dirPDF);
            
            BPTVertex &vtxNextToLast = vertices[vertices.size() - 2];
            vtxNextToLast.revAreaPDF = revInfo.dirPDF * cosLast / dist2;
            vtxNextToLast.revRRProb = std::min((revInfo.fs * absDot(dirOut_sn, gNorm_sn) / revInfo.dirPDF).importance(wlHint), 1.0f);
            
            cosLast = cosIn;
            dirPDF = fsResult.dirPDF;
            sampledType = fsResult.dirType;
            isect = Intersection();
        }
    }
    
    // calculate power heuristic MIS weight
    float BidirectionalPathTracingRenderer::Job::calculateMISWeight(float lExtend1stAreaPDF, float lExtend1stRRProb, float lExtend2ndAreaPDF, float lExtend2ndRRProb,
                                                                    float eExtend1stAreaPDF, float eExtend1stRRProb, float eExtend2ndAreaPDF, float eExtend2ndRRProb,
                                                                    uint32_t numLVtx, uint32_t numEVtx) const {
        const uint32_t minEyeVertices = 1;
        const uint32_t minLightVertices = 0;
        FloatSum recMISWeight = 1;
        
        // extend/shorten light/eye subpath, not consider implicit light subpath reaching a lens.
        if (numEVtx > minEyeVertices) {
            const BPTVertex &eyeEndVtx = eyeVertices[numEVtx - 1];
            float PDFRatio = lExtend1stAreaPDF * lExtend1stRRProb / (eyeEndVtx.areaPDF * eyeEndVtx.RRProb);
            bool shortenIsDeltaSampled = eyeEndVtx.sampledType.isDelta();
            if (!shortenIsDeltaSampled)
                recMISWeight += PDFRatio * PDFRatio;
            bool prevIsDeltaSampled = shortenIsDeltaSampled;
            if (numEVtx - 1 > minEyeVertices) {
                const BPTVertex &newLightVtx = eyeVertices[numEVtx - 2];
                PDFRatio *= lExtend2ndAreaPDF * lExtend2ndRRProb / (newLightVtx.areaPDF * newLightVtx.RRProb);
                shortenIsDeltaSampled = newLightVtx.sampledType.isDelta();
                if (!shortenIsDeltaSampled && !prevIsDeltaSampled)
                    recMISWeight += PDFRatio * PDFRatio;
                prevIsDeltaSampled = shortenIsDeltaSampled;
                for (int t = numEVtx - 2; t > minEyeVertices; --t) {
                    const BPTVertex &newLightVtx = eyeVertices[t - 1];
                    PDFRatio *= newLightVtx.revAreaPDF * newLightVtx.revRRProb / (newLightVtx.areaPDF * newLightVtx.RRProb);
                    shortenIsDeltaSampled = newLightVtx.sampledType.isDelta();
                    if (!shortenIsDeltaSampled && !prevIsDeltaSampled)
                        recMISWeight += PDFRatio * PDFRatio;
                    prevIsDeltaSampled = shortenIsDeltaSampled;
                }
            }
        }
        
        // extend/shorten eye/light subpath, consider implicit eye subpath reaching a light.
        if (numLVtx > minLightVertices) {
            const BPTVertex &lightEndVtx = lightVertices[numLVtx - 1];
            float PDFRatio = eExtend1stAreaPDF * eExtend1stRRProb / (lightEndVtx.areaPDF * lightEndVtx.RRProb);
            bool shortenIsDeltaSampled = lightEndVtx.sampledType.isDelta();
            if (!shortenIsDeltaSampled)
                recMISWeight += PDFRatio * PDFRatio;
            bool prevIsDeltaSampled = shortenIsDeltaSampled;
            if (numLVtx - 1 > minLightVertices) {
                const BPTVertex &newEyeVtx = lightVertices[numLVtx - 2];
                PDFRatio *= eExtend2ndAreaPDF * eExtend2ndRRProb / (newEyeVtx.areaPDF * newEyeVtx.RRProb);
                shortenIsDeltaSampled = newEyeVtx.sampledType.isDelta();
                if (!shortenIsDeltaSampled && !prevIsDeltaSampled)
                    recMISWeight += PDFRatio * PDFRatio;
                prevIsDeltaSampled = shortenIsDeltaSampled;
                for (int s = numLVtx - 2; s > minLightVertices; --s) {
                    const BPTVertex &newEyeVtx = lightVertices[s - 1];
                    PDFRatio *= newEyeVtx.revAreaPDF * newEyeVtx.revRRProb / (newEyeVtx.areaPDF * newEyeVtx.RRProb);
                    shortenIsDeltaSampled = newEyeVtx.sampledType.isDelta();
                    if (!shortenIsDeltaSampled && !prevIsDeltaSampled)
                        recMISWeight += PDFRatio * PDFRatio;
                    prevIsDeltaSampled = shortenIsDeltaSampled;
                }
            }
        }
        
        return 1.0f / recMISWeight;
    }
}

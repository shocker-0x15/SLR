//
//  BidirectionalPathTracingRenderer.cpp
//  SLR
//
//  Created by 渡部 心 on 2016/01/31.
//  Copyright © 2016年 渡部 心. All rights reserved.
//

#include "BidirectionalPathTracingRenderer.h"

#include "../Core/RenderSettings.h"
#include "../Helper/ThreadPool.h"
#include "../Core/XORShift.h"
#include "../Core/ImageSensor.h"
#include "../Core/RandomNumberGenerator.h"
#include "../Memory/ArenaAllocator.h"
#include "../Core/cameras.h"
#include "../Core/SurfaceObject.h"

namespace SLR {
    BidirectionalPathTracingRenderer::BidirectionalPathTracingRenderer(uint32_t spp) : m_samplesPerPixel(spp) {
    }
    
    void BidirectionalPathTracingRenderer::render(const Scene &scene, const RenderSettings &settings) const {
#ifdef DEBUG
        uint32_t numThreads = 1;
#else
        uint32_t numThreads = std::thread::hardware_concurrency();
#endif
        XORShift topRand(settings.getInt(RenderSettingItem::RNGSeed));
        std::unique_ptr<ArenaAllocator[]> mems = std::unique_ptr<ArenaAllocator[]>(new ArenaAllocator[numThreads]);
        std::unique_ptr<XORShift[]> rngs = std::unique_ptr<XORShift[]>(new XORShift[numThreads]);
        for (int i = 0; i < numThreads; ++i) {
            new (mems.get() + i) ArenaAllocator();
            new (rngs.get() + i) XORShift(topRand.getUInt());
        }
        std::unique_ptr<RandomNumberGenerator*[]> rngRefs = std::unique_ptr<RandomNumberGenerator*[]>(new RandomNumberGenerator*[numThreads]);
        for (int i = 0; i < numThreads; ++i)
            rngRefs[i] = &rngs[i];
        
        const Camera* camera = scene.getCamera();
        ImageSensor* sensor = camera->getSensor();
        
        Job job;
        job.scene = &scene;
        
        job.mems = mems.get();
        job.rngs = rngRefs.get();
        
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
                sensor->saveImage(filename, settings.getFloat(RenderSettingItem::Brightness) / (s + 1));
                printf("%u samples: %s\n", exportPass, filename);
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
        RandomNumberGenerator &rng = *rngs[threadID];
        for (int ly = 0; ly < numPixelY; ++ly) {
            for (int lx = 0; lx < numPixelX; ++lx) {
                float time = timeStart + rng.getFloat0cTo1o() * (timeEnd - timeStart);
                float px = basePixelX + lx + rng.getFloat0cTo1o();
                float py = basePixelY + ly + rng.getFloat0cTo1o();
                
                float selectWLPDF;
                WavelengthSamples wls = WavelengthSamples::createWithEqualOffsets(rng.getFloat0cTo1o(), rng.getFloat0cTo1o(), &selectWLPDF);
                
                // initialize working area for the current pixel.
                curPx = px;
                curPy = py;
                wlHint = wls.selectedLambda;
                eyeVertices.clear();
                lightVertices.clear();
                
                // light subpath generation
                {
                    // select one light from all the lights in the scene.
                    float lightProb;
                    Light light;
                    scene->selectLight(rng.getFloat0cTo1o(), &light, &lightProb);
                    SLRAssert(!std::isnan(lightProb) && !std::isinf(lightProb), "lightProb: unexpected value detected: %f", lightProb);
                    
                    // sample a position (emittance) on the selected light's surface.
                    LightPosQuery lightQuery(time, wls);
                    LightPosSample lightSample(rng.getFloat0cTo1o(), rng.getFloat0cTo1o());
                    LightPosQueryResult lightResult;
                    SampledSpectrum Le0 = light.sample(lightQuery, lightSample, &lightResult);
                    EDF* edf = lightResult.surfPt.createEDF(wls, mem);
                    float lightAreaPDF = lightProb * lightResult.areaPDF;
                    SLRAssert(!std::isnan(lightResult.areaPDF)/* && !std::isinf(lightResult)*/, "areaPDF: unexpected value detected: %f", lightResult.areaPDF);
                    
                    // register the first light vertex.
                    lightVertices.emplace_back(lightResult.surfPt, Vector3D::Zero, Normal3D(0, 0, 1), mem.create<EDFProxy>(edf),
                                               Le0 / lightAreaPDF, lightAreaPDF, 1.0f, WavelengthSamples::Flag(0));
                    
                    // sample a direction from EDF, then create subsequent light subpath vertices by tracing in the scene.
                    EDFQuery edfQuery;
                    EDFSample LeSample(rng.getFloat0cTo1o(), rng.getFloat0cTo1o(), rng.getFloat0cTo1o());
                    EDFQueryResult LeResult;
                    SampledSpectrum Le1 = edf->sample(edfQuery, LeSample, &LeResult);
                    Ray ray(lightResult.surfPt.p + Ray::Epsilon * lightResult.surfPt.gNormal, lightResult.surfPt.shadingFrame.fromLocal(LeResult.dir_sn), time);
                    SampledSpectrum alpha = (Le0 * Le1) * (LeResult.dir_sn.z / (lightProb * lightResult.areaPDF * LeResult.dirPDF));
                    generateSubPath(wls, alpha, ray, LeResult.dirPDF, LeResult.dir_sn.z, true, rng, mem);
                }
                
                // eye subpath generation
                {
                    // sample a position (We0, spatial importance) on the lens surface of the camera.
                    LensPosQuery lensQuery(time, wls);
                    LensPosSample lensSample(rng.getFloat0cTo1o(), rng.getFloat0cTo1o());
                    LensPosQueryResult lensResult;
                    SampledSpectrum We0 = camera->sample(lensQuery, lensSample, &lensResult);
                    IDF* idf = camera->createIDF(lensResult.surfPt, wls, mem);
                    
                    // register the first eye vertex.
                    eyeVertices.emplace_back(lensResult.surfPt, Vector3D::Zero, Normal3D(0, 0, 1), mem.create<IDFProxy>(idf),
                                             We0 / (lensResult.areaPDF * selectWLPDF), lensResult.areaPDF, 1.0f, WavelengthSamples::Flag(0));
                    
                    // sample a direction (directional importance) from IDF, then create subsequent eye subpath vertices by tracing in the scene.
                    IDFSample WeSample(px / imageWidth, py / imageHeight);
                    IDFQueryResult WeResult;
                    SampledSpectrum We1 = idf->sample(WeSample, &WeResult);
                    Ray ray(lensResult.surfPt.p, lensResult.surfPt.shadingFrame.fromLocal(WeResult.dirLocal), time);
                    SampledSpectrum alpha = eyeVertices.back().alpha * We1 * (WeResult.dirLocal.z / WeResult.dirPDF);
                    generateSubPath(wls, alpha, ray, WeResult.dirPDF, WeResult.dirLocal.z, false, rng, mem);
                }
                
                // connection
                for (int t = 1; t <= eyeVertices.size(); ++t) {
                    const BPTVertex &eVtx = eyeVertices[t - 1];
                    for (int s = 1; s <= lightVertices.size(); ++s) {
                        const BPTVertex &lVtx = lightVertices[s - 1];
                        
                        if (!scene->testVisibility(eVtx.surfPt, lVtx.surfPt, time))
                            continue;
                        
                        Vector3D connectionVector = eVtx.surfPt.p - lVtx.surfPt.p;
                        float dist2 = connectionVector.sqLength();
                        connectionVector /= std::sqrt(dist2);
                        
                        Vector3D cVecL = lVtx.surfPt.shadingFrame.toLocal(connectionVector);
                        DDFQuery queryLightEnd{lVtx.dirIn_sn, lVtx.gNormal_sn, wlHint, true};
                        SampledSpectrum ddfL = lVtx.ddf->evaluate(queryLightEnd, cVecL);
                        float cosLightEnd = std::fabs(dot(connectionVector, (Vector3D)lVtx.surfPt.gNormal));

                        Vector3D cVecE = eVtx.surfPt.shadingFrame.toLocal(-connectionVector);
                        DDFQuery queryEyeEnd{eVtx.dirIn_sn, eVtx.gNormal_sn, wlHint, false};
                        SampledSpectrum ddfE = eVtx.ddf->evaluate(queryEyeEnd, cVecE);
                        float cosEyeEnd = std::fabs(dot(connectionVector, (Vector3D)eVtx.surfPt.gNormal));
                        
                        float G = cosEyeEnd * cosLightEnd / dist2;
                        SampledSpectrum connectionTerm = ddfL * G * ddfE;
                        if (connectionTerm == SampledSpectrum::Zero)
                            continue;
                        
                        float lExtendDirPDF = lVtx.ddf->evaluatePDF(queryLightEnd, cVecL);
                        float lExtendAreaPDF = lExtendDirPDF * cosEyeEnd / dist2;
                        float lExtendRRProb = s > 1 ? std::min((ddfL * cosLightEnd / lExtendDirPDF).luminance(), 1.0f) : 1.0f;
                        float eExtendDirPDF = eVtx.ddf->evaluatePDF(queryEyeEnd, cVecE);
                        float eExtendAreaPDF = eExtendDirPDF * cosLightEnd / dist2;
                        float eExtendRRProb = t > 1 ? std::min((ddfE * cosEyeEnd / eExtendDirPDF).luminance(), 1.0f) : 1.0f;
                        
                        float MISWeight = calculateMISWeight(lExtendAreaPDF, lExtendRRProb, eExtendAreaPDF, eExtendRRProb, s, t);
                        SLRAssert(!std::isinf(MISWeight) && !std::isnan(MISWeight), "invalid MIS weight.");
                        SampledSpectrum contribution = MISWeight * lVtx.alpha * connectionTerm * eVtx.alpha;
                        if (t > 1) {
                            sensor->add(px, py, wls, contribution);
                        }
                        else {
                            const IDF* idf = (const IDF*)eVtx.ddf->getDDF();
                            float hitPx, hitPy;
                            idf->calculatePixel(cVecE, &hitPx, &hitPy);
                            sensor->add(threadID, hitPx, hitPy, wls, contribution);
                        }
                    }
                }
                
                mem.reset();
            }
        }
    }
    
    void BidirectionalPathTracingRenderer::Job::generateSubPath(const WavelengthSamples &initWLs, const SampledSpectrum &initAlpha, const SLR::Ray &initRay, float dirPDF, float cosLast,
                                                                bool adjoint, RandomNumberGenerator &rng, SLR::ArenaAllocator &mem) {
        std::vector<BPTVertex> &vertices = adjoint ? lightVertices : eyeVertices;
        
        WavelengthSamples wls = initWLs;
        Ray ray = initRay;
        SampledSpectrum alpha = initAlpha;
        
        Intersection isect;
        SurfacePoint surfPt;
        float RRProb = 1.0f;
        while (scene->intersect(ray, &isect)) {
            isect.getSurfacePoint(&surfPt);
            Vector3D dirOut_sn = surfPt.shadingFrame.toLocal(-ray.dir);
            Normal3D gNorm_sn = surfPt.shadingFrame.toLocal(surfPt.gNormal);
            BSDF* bsdf = surfPt.createBSDF(wls, mem);
            
            float areaPDF = dirPDF * std::fabs(dirOut_sn.z) / (isect.dist * isect.dist);
            vertices.emplace_back(surfPt, dirOut_sn, gNorm_sn, mem.create<BSDFProxy>(bsdf), alpha, areaPDF, RRProb, wls.flags);
            
            // implicit path (zero light subpath vertices, s = 0)
            if (!adjoint && surfPt.isEmitting()) {
                EDF* edf = surfPt.createEDF(wls, mem);
                SampledSpectrum Le = surfPt.emittance(wls) * edf->evaluate(EDFQuery(), dirOut_sn);
                
                float lightProb = scene->evaluateProb(Light(isect.obj));
                float lightPDF = lightProb * surfPt.evaluateAreaPDF();
                float lExtendAreaPDF = lightProb * lightPDF;
                float MISWeight = calculateMISWeight(lExtendAreaPDF, 1.0f, 0.0f, 0.0f, 0, (uint32_t)vertices.size());
                SLRAssert(!std::isinf(MISWeight) && !std::isnan(MISWeight), "invalid MIS weight.");
                SampledSpectrum contribution = MISWeight * alpha * Le;
                sensor->add(curPx, curPy, wls, contribution);
            }
            
            if (surfPt.atInfinity)
                break;
            
            BSDFQuery fsQuery(dirOut_sn, gNorm_sn, wls.selectedLambda, DirectionType::All, adjoint);
            BSDFSample fsSample(rng.getFloat0cTo1o(), rng.getFloat0cTo1o(), rng.getFloat0cTo1o());
            BSDFQueryResult fsResult;
            BSDFReverseInfo revInfo;
            fsResult.reverse = &revInfo;
            SampledSpectrum fs = bsdf->sample(fsQuery, fsSample, &fsResult);
            if (fs == SampledSpectrum::Zero || fsResult.dirPDF == 0.0f)
                break;
            if (fsResult.dirType.isDispersive())
                wls.flags |= WavelengthSamples::LambdaIsSelected;
            dirPDF = fsResult.dirPDF;
            cosLast = std::fabs(dot(fsResult.dir_sn, (Vector3D)gNorm_sn));
            SampledSpectrum weight = fs * (cosLast / dirPDF);
            alpha *= weight;
            SLRAssert(!alpha.hasInf() && !alpha.hasNaN(),
                      "alpha: unexpected value detected:\nalpha: %s\nfs: %s\nlength: %u, cos: %g, dirPDF: %g",
                      alpha.toString().c_str(), fs.toString().c_str(), uint32_t(vertices.size()) - 1, std::fabs(fsResult.dir_sn.z), fsResult.dirPDF);
            
            Vector3D dirIn = surfPt.shadingFrame.fromLocal(fsResult.dir_sn);
            ray = Ray(surfPt.p + Ray::Epsilon * dirIn, dirIn, ray.time);
            
            // Russian roulette
            RRProb = std::min(weight.luminance(), 1.0f);
            if (rng.getFloat0cTo1o() < RRProb)
                alpha /= RRProb;
            else
                break;
            
            BPTVertex &vtxNextToLast = vertices[vertices.size() - 2];
            vtxNextToLast.revAreaPDF = revInfo.dirPDF * cosLast / (isect.dist * isect.dist);
            vtxNextToLast.revRRProb = std::min((revInfo.fs * std::fabs(dot(dirOut_sn, (Vector3D)gNorm_sn)) / revInfo.dirPDF).luminance(), 1.0f);
            SLRAssert(!std::isnan(vtxNextToLast.revAreaPDF) && !std::isinf(vtxNextToLast.revAreaPDF), "invalid reverse area PDF.");
            SLRAssert(!std::isnan(vtxNextToLast.revRRProb) && !std::isinf(vtxNextToLast.revRRProb), "invalid reverse RR probability.");
        }
    }
    
    // calculate power heuristic MIS weight
    float BidirectionalPathTracingRenderer::Job::calculateMISWeight(float lExtendAreaPDF, float lExtendRRProb, float eExtendAreaPDF, float eExtendRRProb,
                                                                    uint32_t numLVtx, uint32_t numEVtx) const {
        FloatSum recMISWeight = 1;
        float PDFRatio;
        
        // extend/shorten light/eye subpath, not consider implicit light subpath reaching a lens.
        if (numEVtx > 1) {
            const BPTVertex &eyeEndVtx = eyeVertices[numEVtx - 1];
            PDFRatio = lExtendAreaPDF * lExtendRRProb / (eyeEndVtx.areaPDF * eyeEndVtx.RRProb);
            recMISWeight += PDFRatio * PDFRatio;
            for (int s = numEVtx - 1; s >= 2; --s) {
                const BPTVertex &newLightVtx = eyeVertices[s - 1];
                PDFRatio *= newLightVtx.revAreaPDF * newLightVtx.revRRProb / (newLightVtx.areaPDF * newLightVtx.RRProb);
                recMISWeight += PDFRatio * PDFRatio;
            }
        }
        
        // extend/shorten eye/light subpath, consider implicit eye subpath reaching a light.
        if (numLVtx > 0) {
            const BPTVertex &lightEndVtx = lightVertices[numLVtx - 1];
            PDFRatio = eExtendAreaPDF * eExtendRRProb / (lightEndVtx.areaPDF * lightEndVtx.RRProb);
            recMISWeight += PDFRatio * PDFRatio;
            for (int t = numLVtx - 1; t >= 1; --t) {
                const BPTVertex &newEyeVtx = lightVertices[t - 1];
                PDFRatio *= newEyeVtx.revAreaPDF * newEyeVtx.revRRProb / (newEyeVtx.areaPDF * newEyeVtx.RRProb);
                recMISWeight += PDFRatio * PDFRatio;
            }
        }
        
        return 1.0f / recMISWeight;
    }
}

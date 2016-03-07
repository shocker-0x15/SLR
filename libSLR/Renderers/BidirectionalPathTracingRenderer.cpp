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
                WavelengthSamples wls = WavelengthSamples::sampleUniform(rng.getFloat0cTo1o(), &selectWLPDF);
                
                // initialize working area for the current pixel.
                heroIndex = wls.heroIdx;
                wlPDFs = selectWLPDF / WavelengthSamples::NumComponents;
                curPx = px;
                curPy = py;
                eyeVertices.clear();
                lightVertices.clear();
                
                // light subpath generation
                {
                    // select one light from all the lights in the scene.
                    float lightProb;
                    Light light;
                    scene->selectLight(rng.getFloat0cTo1o(), &light, &lightProb);
                    SLRAssert(!std::isnan(lightProb) && !std::isinf(lightProb), "lightProb: unexpected value detected: %f", lightProb);
                    
                    // sample a ray with its radiance (emittance, EDF value) from the selected light.
                    LightPosQuery lightPosQuery(time, wls);
                    LightPosSample lightPosSample(rng.getFloat0cTo1o(), rng.getFloat0cTo1o());
                    LightPosQueryResult lightPosResult;
                    EDFQuery edfQuery(wls.heroIdx);
                    EDFSample edfSample(rng.getFloat0cTo1o(), rng.getFloat0cTo1o(), rng.getFloat0cTo1o());
                    EDFQueryResult edfResult;
                    EDF* edf;
                    SampledSpectrum Le0, Le1;
                    Ray ray = light.sampleRay(lightPosQuery, lightPosSample, &lightPosResult, &Le0, &edf, edfQuery, edfSample, &edfResult, &Le1, mem);
                    
                    // register the first light vertex.
                    SampledSpectrum lightAreaPDF = lightProb * lightPosResult.areaPDF;
                    SampledSpectrum wlLightPDFs = lightAreaPDF;
                    lightVertices.emplace_back(lightPosResult.surfPt, Vector3D::Zero, Normal3D(0, 0, 1), mem.create<EDFProxy>(edf),
                                               wlLightPDFs, Le0 / lightAreaPDF[wls.heroIdx], lightAreaPDF, SampledSpectrum::One, lightPosResult.posType);
                    
                    // create subsequent light subpath vertices by tracing in the scene.
                    SampledSpectrum alpha = lightVertices.back().alpha * Le1 * (absDot(ray.dir, lightPosResult.surfPt.gNormal) / edfResult.dirPDF[wls.heroIdx]);
                    generateSubPath(wls, wlLightPDFs, alpha, ray, edfResult.dirPDF, edfResult.dirType, edfResult.dir_sn.z, true, rng, mem);
                }
                
                // eye subpath generation
                {
                    // sample a ray with its importances (spatial, directional) from the lens and its IDF.
                    LensPosQuery lensQuery(time, wls);
                    LensPosSample lensSample(rng.getFloat0cTo1o(), rng.getFloat0cTo1o());
                    LensPosQueryResult lensResult;
                    IDFSample WeSample(px / imageWidth, py / imageHeight);
                    IDFQueryResult WeResult;
                    IDF* idf;
                    SampledSpectrum We0, We1;
                    Ray ray = camera->sampleRay(lensQuery, lensSample, &lensResult, &We0, &idf, WeSample, &WeResult, &We1, mem);
                    
                    // register the first eye vertex.
                    SampledSpectrum wlEyePDFs = lensResult.areaPDF;
                    eyeVertices.emplace_back(lensResult.surfPt, Vector3D::Zero, Normal3D(0, 0, 1), mem.create<IDFProxy>(idf),
                                             wlEyePDFs, We0 / lensResult.areaPDF[wls.heroIdx], lensResult.areaPDF, SampledSpectrum::One, lensResult.posType);
                    
                    // create subsequent eye subpath vertices by tracing in the scene.
                    SampledSpectrum alpha = eyeVertices.back().alpha * We1 * (absDot(ray.dir, lensResult.surfPt.gNormal) / WeResult.dirPDF[wls.heroIdx]);
                    generateSubPath(wls, wlEyePDFs, alpha, ray, WeResult.dirPDF, WeResult.dirType, WeResult.dirLocal.z, false, rng, mem);
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
                        DDFQuery queryLightEnd{lVtx.dirIn_sn, lVtx.gNormal_sn, wls.heroIdx, true};
                        SampledSpectrum lRevDDF, eExtend2ndDirPDF;
                        SampledSpectrum lDDF = lVtx.ddf->evaluate(queryLightEnd, lConnectVector, &lRevDDF);
                        SampledSpectrum lExtend1stDirPDF = lVtx.ddf->evaluatePDF(queryLightEnd, lConnectVector, &eExtend2ndDirPDF);
                        
                        Vector3D eConnectVector = eVtx.surfPt.shadingFrame.toLocal(connectionVector);
                        DDFQuery queryEyeEnd{eVtx.dirIn_sn, eVtx.gNormal_sn, wls.heroIdx, false};
                        SampledSpectrum eRevDDF, lExtend2ndDirPDF;
                        SampledSpectrum eDDF = eVtx.ddf->evaluate(queryEyeEnd, eConnectVector, &eRevDDF);
                        SampledSpectrum eExtend1stDirPDF = eVtx.ddf->evaluatePDF(queryEyeEnd, eConnectVector, &lExtend2ndDirPDF);
                        
                        SampledSpectrum connectionTerm = lDDF * G * eDDF;
                        if (connectionTerm == SampledSpectrum::Zero)
                            continue;
                        
                        if (!scene->testVisibility(eVtx.surfPt, lVtx.surfPt, time))
                            continue;
                        // ----------------------------------------------------------------
                        
                        // calculate the 1st and 2nd subpath extending PDFs and probabilities.
                        // They can't be stored in advance because they depend on the connection.
                        SampledSpectrum lExtend1stAreaPDF, lExtend1stRRProb, lExtend2ndAreaPDF, lExtend2ndRRProb;
                        {
                            lExtend1stAreaPDF = lExtend1stDirPDF * cosEyeEnd / connectDist2;
                            lExtend1stRRProb = s > 1 ? min(lDDF * cosLightEnd / lExtend1stDirPDF, SampledSpectrum::One) : SampledSpectrum::One;
                            
                            if (t > 1) {
                                BPTVertex &eVtxNextToEnd = eyeVertices[t - 2];
                                float dist2;
                                Vector3D dir2nd = eVtx.surfPt.getDirectionFrom(eVtxNextToEnd.surfPt.p, &dist2);
                                lExtend2ndAreaPDF = lExtend2ndDirPDF * absDot(eVtxNextToEnd.surfPt.gNormal, dir2nd) / dist2;
                                lExtend2ndRRProb = min(eRevDDF * absDot(eVtx.gNormal_sn, eVtx.dirIn_sn) / lExtend2ndDirPDF, SampledSpectrum::One);
                            }
                        }
                        SampledSpectrum eExtend1stAreaPDF, eExtend1stRRProb, eExtend2ndAreaPDF, eExtend2ndRRProb;
                        {
                            eExtend1stAreaPDF = eExtend1stDirPDF * cosLightEnd / connectDist2;
                            eExtend1stRRProb = t > 1 ? min(eDDF * cosEyeEnd / eExtend1stDirPDF, SampledSpectrum::One) : SampledSpectrum::One;
                            
                            if (s > 1) {
                                BPTVertex &lVtxNextToEnd = lightVertices[s - 2];
                                float dist2;
                                Vector3D dir2nd = lVtxNextToEnd.surfPt.getDirectionFrom(lVtx.surfPt.p, &dist2);
                                eExtend2ndAreaPDF = eExtend2ndDirPDF * absDot(lVtxNextToEnd.surfPt.gNormal, dir2nd) / dist2;
                                eExtend2ndRRProb = min(lRevDDF * absDot(lVtx.gNormal_sn, lVtx.dirIn_sn) / eExtend2ndDirPDF, SampledSpectrum::One);
                            }
                        }
                        
                        // calculate MIS weight and store weighted contribution to a sensor.
                        float MISWeight = calculateMISWeight(lExtend1stAreaPDF, lExtend1stRRProb, lExtend2ndAreaPDF, lExtend2ndRRProb,
                                                             eExtend1stAreaPDF, eExtend1stRRProb, eExtend2ndAreaPDF, eExtend2ndRRProb, s, t);
                        if (std::isinf(MISWeight) || std::isnan(MISWeight))
                            continue;
                        SLRAssert(MISWeight >= 0 && MISWeight <= 1.0f, "invalid MIS weight: %g", MISWeight);
                        SampledSpectrum contribution = MISWeight * lVtx.alpha * connectionTerm * eVtx.alpha / wlPDFs[heroIndex];
                        SLRAssert(contribution.hasNaN() == false && contribution.hasInf() == false && contribution.hasMinus() == false,
                                  "Unexpected value detected: %s\n"
                                  "pix: (%f, %f)", contribution.toString().c_str(), px, py);
                        if (t > 1) {
                            sensor->add(px, py, wls, contribution);
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
    
    void BidirectionalPathTracingRenderer::Job::generateSubPath(const WavelengthSamples &initWLs, const SampledSpectrum &initWLPDFs,
                                                                const SampledSpectrum &initAlpha, const SLR::Ray &initRay, const SampledSpectrum &initDirPDF, DirectionType sampledType,
                                                                float cosLast, bool adjoint, RandomNumberGenerator &rng, SLR::ArenaAllocator &mem) {
        std::vector<BPTVertex> &vertices = adjoint ? lightVertices : eyeVertices;
        
        // reject invalid values.
        if (initDirPDF[initWLs.heroIdx] == 0.0f)
            return;
        
        WavelengthSamples wls = initWLs;
        SampledSpectrum subpathWLPDFs = initWLPDFs;
        SampledSpectrum alpha = initAlpha;
        Ray ray = initRay;
        SampledSpectrum dirPDF = initDirPDF;
        
        Intersection isect;
        SurfacePoint surfPt;
        SampledSpectrum RRProb = 1.0f;
        while (scene->intersect(ray, &isect)) {
            isect.getSurfacePoint(&surfPt);
            float dist2 = squaredDistance(vertices.back().surfPt, surfPt);
            Vector3D dirOut_sn = surfPt.shadingFrame.toLocal(-ray.dir);
            Normal3D gNorm_sn = surfPt.shadingFrame.toLocal(surfPt.gNormal);
            BSDF* bsdf = surfPt.createBSDF(wls, mem);
            
            SampledSpectrum areaPDF = dirPDF * absDot(dirOut_sn, gNorm_sn) / dist2;
            subpathWLPDFs *= areaPDF * RRProb;
            vertices.emplace_back(surfPt, dirOut_sn, gNorm_sn, mem.create<BSDFProxy>(bsdf), subpathWLPDFs, alpha, areaPDF, RRProb, sampledType);
            
            // implicit path (zero light subpath vertices, s = 0)
            if (!adjoint && surfPt.isEmitting()) {
                EDF* edf = surfPt.createEDF(wls, mem);
                SampledSpectrum Le0 = surfPt.emittance(wls);
                SampledSpectrum Le1 = edf->evaluate(EDFQuery(wls.heroIdx), dirOut_sn);
                
                float lightProb = scene->evaluateProb(Light(isect.obj));
                SampledSpectrum extend1stAreaPDF = lightProb * surfPt.evaluateAreaPDF();
                SampledSpectrum extend2ndAreaPDF = edf->evaluatePDF(EDFQuery(wls.heroIdx), dirOut_sn) * cosLast / dist2;
                
                float MISWeight = calculateMISWeight(extend1stAreaPDF, SampledSpectrum::One, extend2ndAreaPDF, SampledSpectrum::One,
                                                     SampledSpectrum::Zero, SampledSpectrum::Zero, SampledSpectrum::Zero, SampledSpectrum::Zero,
                                                     0, (uint32_t)vertices.size());
                if (!std::isinf(MISWeight) && !std::isnan(MISWeight)) {
                    SLRAssert(MISWeight >= 0 && MISWeight <= 1.0f, "invalid MIS weight: %g", MISWeight);
                    SampledSpectrum contribution = MISWeight * alpha * Le0 * Le1 / wlPDFs[heroIndex];
                    SLRAssert(contribution.hasNaN() == false && contribution.hasInf() == false && contribution.hasMinus() == false,
                              "Unexpected value detected: %s\n"
                              "pix: (%f, %f)", contribution.toString().c_str(), curPx, curPy);
                    sensor->add(curPx, curPy, wls, contribution);
                }
            }
            
            if (surfPt.atInfinity) {
                vertices.pop_back();
                break;
            }
            
            BSDFQuery fsQuery(wls.heroIdx, dirOut_sn, gNorm_sn, DirectionType::All, adjoint);
            BSDFSample fsSample(rng.getFloat0cTo1o(), rng.getFloat0cTo1o(), rng.getFloat0cTo1o());
            BSDFQueryResult fsResult;
            BSDFReverseInfo revInfo;
            fsResult.reverse = &revInfo;
            SampledSpectrum fs = bsdf->sample(fsQuery, fsSample, &fsResult);
            if (fs == SampledSpectrum::Zero || fsResult.dirPDF[wls.heroIdx] == 0.0f)
                break;
            float cosIn = absDot(fsResult.dir_sn, gNorm_sn);
            
            // Russian roulette
            RRProb = min(positiveMask(fs / fsResult.dirPDF * cosIn, fs), SampledSpectrum::One);
            if (rng.getFloat0cTo1o() >= RRProb[wls.heroIdx])
                break;
            
            SampledSpectrum weight = fs * (cosIn / fsResult.dirPDF[wls.heroIdx]);
            alpha *= weight / RRProb[wls.heroIdx];
            ray = Ray(surfPt.p, surfPt.shadingFrame.fromLocal(fsResult.dir_sn), ray.time, Ray::Epsilon);
            SLRAssert(!weight.hasInf() && !weight.hasNaN(),
                      "weight: unexpected value detected:\nweight: %s\nfs: %s\nlength: %u, cos: %g, dirPDF: %s",
                      weight.toString().c_str(), fs.toString().c_str(), uint32_t(vertices.size()) - 1, cosIn, fsResult.dirPDF.toString().c_str());
            
            BPTVertex &vtxNextToLast = vertices[vertices.size() - 2];
            vtxNextToLast.revAreaPDF = revInfo.dirPDF * cosLast / dist2;
            vtxNextToLast.revRRProb = min(positiveMask(revInfo.fs * absDot(dirOut_sn, gNorm_sn) / revInfo.dirPDF, revInfo.fs), SampledSpectrum::One);
            
            cosLast = cosIn;
            dirPDF = fsResult.dirPDF;
            sampledType = fsResult.dirType;
            isect = Intersection();
        }
    }
    
    // calculate power heuristic MIS weight
    float BidirectionalPathTracingRenderer::Job::calculateMISWeight(const SampledSpectrum &lExtend1stAreaPDF, const SampledSpectrum &lExtend1stRRProb,
                                                                    const SampledSpectrum &lExtend2ndAreaPDF, const SampledSpectrum &lExtend2ndRRProb,
                                                                    const SampledSpectrum &eExtend1stAreaPDF, const SampledSpectrum &eExtend1stRRProb,
                                                                    const SampledSpectrum &eExtend2ndAreaPDF, const SampledSpectrum &eExtend2ndRRProb,
                                                                    uint32_t numLVtx, uint32_t numEVtx) const {
        const uint32_t minEyeVertices = 1;
        const uint32_t minLightVertices = 0;
        FloatSum denomMISWeight = 0;
        
        const auto accumulateWLVariants = [](FloatSum &sum, const SampledSpectrum &PDFs) {
            for (int i = 0; i < WavelengthSamples::NumComponents; ++i)
                sum += PDFs[i] * PDFs[i];
        };
        
        const SampledSpectrum &eyeEndWLPDFs = numEVtx > 1 ? eyeVertices[numEVtx - 1].wlPDFs : SampledSpectrum::One;
        const SampledSpectrum &lightEndWLPDFs = numLVtx > 1 ? lightVertices[numLVtx - 1].wlPDFs : SampledSpectrum::One;
        const SampledSpectrum weightingPathPDFs = wlPDFs * lightEndWLPDFs * eyeEndWLPDFs;
        float numMISWeight = weightingPathPDFs[heroIndex] * weightingPathPDFs[heroIndex];
        accumulateWLVariants(denomMISWeight, weightingPathPDFs);
        
        // extend/shorten light/eye subpath, not consider implicit light subpath reaching a lens.
        if (numEVtx > minEyeVertices) {
            const BPTVertex &eyeEndVtx = eyeVertices[numEVtx - 1];
            SampledSpectrum extendedPathPDF = positiveMask(weightingPathPDFs * lExtend1stAreaPDF * lExtend1stRRProb / (eyeEndVtx.areaPDF * eyeEndVtx.RRProb), weightingPathPDFs);
            bool shortenIsDeltaSampled = eyeEndVtx.sampledType.isDelta();
            if (!shortenIsDeltaSampled)
                accumulateWLVariants(denomMISWeight, extendedPathPDF);
            bool prevIsDeltaSampled = shortenIsDeltaSampled;
            
            if (numEVtx - 1 > minEyeVertices) {
                const BPTVertex &newLightVtx = eyeVertices[numEVtx - 2];
                extendedPathPDF *= positiveMask(lExtend2ndAreaPDF * lExtend2ndRRProb / (newLightVtx.areaPDF * newLightVtx.RRProb), weightingPathPDFs);
                shortenIsDeltaSampled = newLightVtx.sampledType.isDelta();
                if (!shortenIsDeltaSampled && !prevIsDeltaSampled)
                    accumulateWLVariants(denomMISWeight, extendedPathPDF);
                prevIsDeltaSampled = shortenIsDeltaSampled;
                
                for (int t = numEVtx - 2; t > minEyeVertices; --t) {
                    const BPTVertex &newLightVtx = eyeVertices[t - 1];
                    extendedPathPDF *= positiveMask(newLightVtx.revAreaPDF * newLightVtx.revRRProb / (newLightVtx.areaPDF * newLightVtx.RRProb), weightingPathPDFs);
                    shortenIsDeltaSampled = newLightVtx.sampledType.isDelta();
                    if (!shortenIsDeltaSampled && !prevIsDeltaSampled)
                        accumulateWLVariants(denomMISWeight, extendedPathPDF);
                    prevIsDeltaSampled = shortenIsDeltaSampled;
                }
            }
        }
        
        // extend/shorten eye/light subpath, consider implicit eye subpath reaching a light.
        if (numLVtx > minLightVertices) {
            const BPTVertex &lightEndVtx = lightVertices[numLVtx - 1];
            SampledSpectrum extendedPathPDF = positiveMask(weightingPathPDFs * eExtend1stAreaPDF * eExtend1stRRProb / (lightEndVtx.areaPDF * lightEndVtx.RRProb), weightingPathPDFs);
            bool shortenIsDeltaSampled = lightEndVtx.sampledType.isDelta();
            if (!shortenIsDeltaSampled)
                accumulateWLVariants(denomMISWeight, extendedPathPDF);
            bool prevIsDeltaSampled = shortenIsDeltaSampled;
            
            if (numLVtx - 1 > minLightVertices) {
                const BPTVertex &newEyeVtx = lightVertices[numLVtx - 2];
                extendedPathPDF *= positiveMask(eExtend2ndAreaPDF * eExtend2ndRRProb / (newEyeVtx.areaPDF * newEyeVtx.RRProb), weightingPathPDFs);
                shortenIsDeltaSampled = newEyeVtx.sampledType.isDelta();
                if (!shortenIsDeltaSampled && !prevIsDeltaSampled)
                    accumulateWLVariants(denomMISWeight, extendedPathPDF);
                prevIsDeltaSampled = shortenIsDeltaSampled;
                
                for (int s = numLVtx - 2; s > minLightVertices; --s) {
                    const BPTVertex &newEyeVtx = lightVertices[s - 1];
                    extendedPathPDF *= positiveMask(newEyeVtx.revAreaPDF * newEyeVtx.revRRProb / (newEyeVtx.areaPDF * newEyeVtx.RRProb), weightingPathPDFs);
                    shortenIsDeltaSampled = newEyeVtx.sampledType.isDelta();
                    if (!shortenIsDeltaSampled && !prevIsDeltaSampled)
                        accumulateWLVariants(denomMISWeight, extendedPathPDF);
                    prevIsDeltaSampled = shortenIsDeltaSampled;
                }
            }
        }
        
        return numMISWeight / denomMISWeight;
    }
}

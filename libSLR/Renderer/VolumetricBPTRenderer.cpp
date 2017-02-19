//
//  VolumetricBPTRenderer.cpp
//
//  Created by 渡部 心 on 2017/02/16.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "VolumetricBPTRenderer.h"

#include "../Core/RenderSettings.h"
#include "../Core/ImageSensor.h"
#include "../Core/random_number_generator.h"
#include "../Core/camera.h"
#include "../Core/light_path_sampler.h"
#include "../Core/ProgressReporter.h"
#include "../Scene/Scene.h"
#include "../MemoryAllocators/ArenaAllocator.h"
#include "../Helper/ThreadPool.h"
#include "../RNG/XORShiftRNG.h"

namespace SLR {
    VolumetricBPTRenderer::VolumetricBPTRenderer(uint32_t spp) : m_samplesPerPixel(spp) {
    }
    
    void VolumetricBPTRenderer::render(const Scene &scene, const RenderSettings &settings) const {
        uint32_t numThreads = settings.getInt(RenderSettingItem::NumThreads);
        XORShiftRNG topRand(settings.getInt(RenderSettingItem::RNGSeed));
        ArenaAllocator* mems = new ArenaAllocator[numThreads];
        IndependentLightPathSampler* samplers = new IndependentLightPathSampler[numThreads];
        for (int i = 0; i < numThreads; ++i) {
            new (mems + i) ArenaAllocator();
            new (samplers + i) IndependentLightPathSampler(topRand.getUInt());
        }
        
        const Camera* camera = scene.getCamera();
        ImageSensor* sensor = camera->getSensor();
        
        Job job;
        job.scene = &scene;
        
        job.mems = mems;
        job.pathSamplers = samplers;
        
        job.camera = camera;
        job.sensor = sensor;
        job.timeStart = settings.getFloat(RenderSettingItem::TimeStart);
        job.timeEnd = settings.getFloat(RenderSettingItem::TimeEnd);
        job.imageWidth = settings.getInt(RenderSettingItem::ImageWidth);
        job.imageHeight = settings.getInt(RenderSettingItem::ImageHeight);
        job.numPixelX = sensor->tileWidth();
        job.numPixelY = sensor->tileHeight();
        
        sensor->init(job.imageWidth, job.imageHeight);
        sensor->addSeparatedBuffers(numThreads);
        
        printf("Volumetric Bidirectional Path Tracing: %u[spp]\n", m_samplesPerPixel);
        ProgressReporter reporter;
        job.reporter = &reporter;
        
        reporter.pushJob("Rendering", m_samplesPerPixel * sensor->numTileX() * sensor->numTileY());
        char nextTitle[32];
        snprintf(nextTitle, sizeof(nextTitle), "To %5uspp", 1);
        reporter.pushJob(nextTitle, 1 * sensor->numTileX() * sensor->numTileY());
        uint32_t imgIdx = 0;
        uint32_t exportPass = 1;
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
                reporter.popJob();
                
                reporter.beginOtherThreadPrint();
                char filename[256];
                sprintf(filename, "%03u.bmp", imgIdx);
                double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(reporter.elapsed()).count();
                sensor->saveImage(filename, settings.getFloat(RenderSettingItem::Brightness) / (s + 1));
                printf("%u samples: %s, %g[s]\n", exportPass, filename, elapsed * 0.001f);
                reporter.endOtherThreadPrint();
                
                ++imgIdx;
                if ((s + 1) == m_samplesPerPixel)
                    break;
                exportPass += exportPass;
                snprintf(nextTitle, sizeof(nextTitle), "To %5uspp", exportPass);
                reporter.pushJob(nextTitle, (exportPass >> 1) * sensor->numTileX() * sensor->numTileY());
            }
        }
        reporter.popJob();
        reporter.finish();
        
        delete[] samplers;
        delete[] mems;
    }
    
    void VolumetricBPTRenderer::Job::kernel(uint32_t threadID) {
        ArenaAllocator &mem = mems[threadID];
        IndependentLightPathSampler &pathSampler = pathSamplers[threadID];
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
                    Light* light;
                    scene->selectLight(pathSampler.getLightSelectionSample(), time, mem, &light, &lightProb);
                    SLRAssert(std::isfinite(lightProb), "lightProb: unexpected value detected: %f", lightProb);
                    
                    // sample a ray with its radiance (emittance, EDF value) from the selected light.
                    LightPosQuery lightPosQuery(time, wls);
                    LightPosQueryResult* lightPosResult;
                    EDFQuery edfQuery;
                    EDFQueryResult edfResult;
                    EDF* edf;
                    SampledSpectrum Le0, Le1;
                    Ray ray = light->sampleRay(lightPosQuery, pathSampler, edfQuery, mem,
                                               &lightPosResult, &Le0, &edf, &edfResult, &Le1);
                    InteractionPoint* interPt = lightPosResult->getInteractionPoint();
                    
                    // register the first light vertex.
                    float lightSpatialPDF = lightProb * lightPosResult->spatialPDF();
                    lightVertices.emplace_back(interPt, mem.create<EDFProxy>(edf, edfQuery),
                                               Le0 / lightSpatialPDF, 0.0f, lightSpatialPDF, 1.0f, lightPosResult->posType, false);
                    
                    // create subsequent light subpath vertices by tracing in the scene.
                    SampledSpectrum alpha = lightVertices.back().alpha * Le1 * (interPt->calcCosTerm(ray.dir) / edfResult.dirPDF);
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
                    InteractionPoint* interPt = mem.create<SurfacePoint>(lensResult.surfPt);
                    
                    // register the first eye vertex.
                    eyeVertices.emplace_back(interPt, mem.create<IDFProxy>(idf),
                                             We0 / (lensResult.areaPDF * selectWLPDF), 0.0f, lensResult.areaPDF, 1.0f, lensResult.posType, false);
                    
                    // create subsequent eye subpath vertices by tracing in the scene.
                    SampledSpectrum alpha = eyeVertices.back().alpha * We1 * (interPt->calcCosTerm(ray.dir) / WeResult.dirPDF);
                    generateSubPath(wls, alpha, ray, WeResult.dirPDF, WeResult.dirType, WeResult.dirLocal.z, false, pathSampler, mem);
                }
                
                // connection
                for (int t = 1; t <= eyeVertices.size(); ++t) {
                    const VBPTVertex &eVtx = eyeVertices[t - 1];
                    for (int s = 1; s <= lightVertices.size(); ++s) {
                        const VBPTVertex &lVtx = lightVertices[s - 1];
                        
                        // calculate the remaining factors of the full path
                        // that are not included in the precomputed weights.
                        // ----------------------------------------------------------------
                        float connectDist2;
                        Vector3D connectionVector = lVtx.interPt->getDirectionFrom(eVtx.interPt->getPosition(), &connectDist2);
                        float cosLightEnd = lVtx.interPt->calcCosTerm(connectionVector);
                        float cosEyeEnd = eVtx.interPt->calcCosTerm(connectionVector);
                        float G = cosEyeEnd * cosLightEnd / connectDist2;
                        
                        Vector3D lConnectVector = lVtx.interPt->toLocal(-connectionVector);
                        SampledSpectrum lRevDDF;
                        SampledSpectrum lDDF = lVtx.ddf->evaluate(lConnectVector, &lRevDDF);
                        float eExtend2ndDirPDF;
                        float lExtend1stDirPDF = lVtx.ddf->evaluatePDF(lConnectVector, &eExtend2ndDirPDF);
                        
                        Vector3D eConnectVector = eVtx.interPt->toLocal(connectionVector);
                        SampledSpectrum eRevDDF;
                        SampledSpectrum eDDF = eVtx.ddf->evaluate(eConnectVector, &eRevDDF);
                        float lExtend2ndDirPDF;
                        float eExtend1stDirPDF = eVtx.ddf->evaluatePDF(eConnectVector, &lExtend2ndDirPDF);
                        
                        float wlProb = 1.0f;
                        if (lVtx.lambdaSelected || eVtx.lambdaSelected)
                            wlProb = 1.0f / WavelengthSamples::NumComponents;
                        SampledSpectrum connectionTerm = lDDF * (G / wlProb) * eDDF;
                        if (connectionTerm == SampledSpectrum::Zero)
                            continue;
                        
                        SampledSpectrum visibility;
                        bool singleWavelength;
                        if (!scene->testVisibility(eVtx.interPt, lVtx.interPt, time, wls, pathSampler, &visibility, &singleWavelength))
                            continue;
                        if (singleWavelength && !wls.lambdaSelected())
                            visibility[wls.selectedLambda] *= WavelengthSamples::NumComponents;
                        connectionTerm *= visibility;
                        // ----------------------------------------------------------------
                        
                        // calculate the 1st and 2nd subpath extending PDFs and probabilities.
                        // They can't be stored in advance because they depend on the connection.
                        float lExtend1stSpatialPDF, lExtend1stRRProb, lExtend2ndSpatialPDF, lExtend2ndRRProb;
                        {
                            lExtend1stSpatialPDF = lExtend1stDirPDF * cosEyeEnd / connectDist2;
                            lExtend1stRRProb = s > 1 ? std::min((lDDF * cosLightEnd / lExtend1stDirPDF).importance(wlHint), 1.0f) : 1.0f;
                            
                            if (t > 1) {
                                VBPTVertex &eVtxNextToEnd = eyeVertices[t - 2];
                                float dist2;
                                Vector3D dir2nd = eVtx.interPt->getDirectionFrom(eVtxNextToEnd.interPt->getPosition(), &dist2);
                                lExtend2ndSpatialPDF = lExtend2ndDirPDF * eVtxNextToEnd.interPt->calcCosTerm(dir2nd) / dist2;
                                lExtend2ndRRProb = std::min((eRevDDF * eVtx.cosIn / lExtend2ndDirPDF).importance(wlHint), 1.0f);
                            }
                        }
                        float eExtend1stSpatialPDF, eExtend1stRRProb, eExtend2ndSpatialPDF, eExtend2ndRRProb;
                        {
                            eExtend1stSpatialPDF = eExtend1stDirPDF * cosLightEnd / connectDist2;
                            eExtend1stRRProb = t > 1 ? std::min((eDDF * cosEyeEnd / eExtend1stDirPDF).importance(wlHint), 1.0f) : 1.0f;
                            
                            if (s > 1) {
                                VBPTVertex &lVtxNextToEnd = lightVertices[s - 2];
                                float dist2;
                                Vector3D dir2nd = lVtxNextToEnd.interPt->getDirectionFrom(lVtx.interPt->getPosition(), &dist2);
                                eExtend2ndSpatialPDF = eExtend2ndDirPDF * lVtxNextToEnd.interPt->calcCosTerm(dir2nd) / dist2;
                                eExtend2ndRRProb = std::min((lRevDDF * lVtx.cosIn / eExtend2ndDirPDF).importance(wlHint), 1.0f);
                            }
                        }
                        
                        // calculate MIS weight and store weighted contribution to a sensor.
                        float MISWeight = calculateMISWeight(lExtend1stSpatialPDF, lExtend1stRRProb, lExtend2ndSpatialPDF, lExtend2ndRRProb,
                                                             eExtend1stSpatialPDF, eExtend1stRRProb, eExtend2ndSpatialPDF, eExtend2ndRRProb, s, t);
                        if (std::isinf(MISWeight) || std::isnan(MISWeight))
                            continue;
                        SLRAssert(MISWeight >= 0 && MISWeight <= 1.0f, "invalid MIS weight: %g", MISWeight);
                        SampledSpectrum contribution = MISWeight * lVtx.alpha * connectionTerm * eVtx.alpha;
                        SLRAssert(contribution.allFinite() && !contribution.hasMinus(),
                                  "Unexpected value detected: %s\n"
                                  "pix: (%f, %f)", contribution.toString().c_str(), p.x, p.y);
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
        reporter->update();
    }
    
    void VolumetricBPTRenderer::Job::generateSubPath(const WavelengthSamples &initWLs, const SampledSpectrum &initAlpha, const Ray &initRay, float dirPDF, DirectionType sampledType,
                                                     float cosLast, bool adjoint, IndependentLightPathSampler &pathSampler, ArenaAllocator &mem) {
        std::vector<VBPTVertex> &vertices = adjoint ? lightVertices : eyeVertices;
        
        // reject invalid values.
        if (dirPDF == 0.0f)
            return;
        
        WavelengthSamples wls = initWLs;
        Ray ray = initRay;
        SampledSpectrum alpha = initAlpha;
        
        Interaction* interact;
        SampledSpectrum medThroughput;
        bool singleWavelength;
        InteractionPoint* interPt;
        
        float RRProb = 1.0f;
        while (scene->interact(ray, wls, pathSampler, mem, &interact, &medThroughput, &singleWavelength)) {
            if (singleWavelength && !wls.lambdaSelected()) {
                medThroughput[wls.selectedLambda] *= WavelengthSamples::NumComponents;
                wls.flags |= WavelengthSamples::LambdaIsSelected;
            }
            alpha *= medThroughput;
            
            interPt = interact->createInteractionPoint(mem);
            alpha *= interPt->evaluateInteractance();
            
            float dist2 = squaredDistance(vertices.back().interPt, interPt);
            Vector3D dirOut_sn = interPt->toLocal(-ray.dir);
            float cosOut = interPt->calcCosTerm(-ray.dir);
            AbstractBDF* abdf = interPt->createAbstractBDF(wls, mem);
            ABDFQuery* abdfQuery = interPt->createABDFQuery(dirOut_sn, wls.selectedLambda, DirectionType::All, true, adjoint, mem);
            
            float spatialPDF = dirPDF * cosOut / dist2;
            vertices.emplace_back(interPt, mem.create<ABDFProxy>(abdf, abdfQuery), alpha, cosOut, spatialPDF, RRProb, sampledType, wls.lambdaSelected());
            
            // implicit path (zero light subpath vertices, s = 0)
            if (!adjoint && interPt->isEmitting()) {
                EDF* edf = interPt->createEDF(wls, mem);
                SampledSpectrum Le0 = interPt->emittance(wls);
                SampledSpectrum Le1 = edf->evaluate(EDFQuery(), dirOut_sn);
                
                float extend1stSpatialPDF = interact->getLightProb() * interPt->evaluateSpatialPDF();
                float extend2ndSpatialPDF = edf->evaluatePDF(EDFQuery(), dirOut_sn) * cosLast / dist2;
                
                float MISWeight = calculateMISWeight(extend1stSpatialPDF, 1.0f, extend2ndSpatialPDF, 1.0f,
                                                     0.0f, 0.0f, 0.0f, 0.0f,
                                                     0, (uint32_t)vertices.size());
                if (!std::isinf(MISWeight) && !std::isnan(MISWeight)) {
                    SampledSpectrum contribution = MISWeight * alpha * Le0 * Le1;
                    SLRAssert(MISWeight >= 0 && MISWeight <= 1.0f, "invalid MIS weight: %g", MISWeight);
                    SLRAssert(contribution.allFinite() && !contribution.hasMinus(),
                              "Unexpected value detected: %s\n"
                              "pix: (%f, %f)", contribution.toString().c_str(), curPx, curPy);
                    if (wls.lambdaSelected())
                        contribution *= WavelengthSamples::NumComponents;
                    sensor->add(curPx, curPy, wls, contribution);
                }
            }
            
            if (interPt->atInfinity()) {
                vertices.pop_back();
                break;
            }
            
            ABDFQueryResult* abdfResult;
            SampledSpectrum abdfValue = abdf->sample(abdfQuery, pathSampler, mem, &abdfResult);
            if (abdfValue == SampledSpectrum::Zero || abdfResult->dirPDF == 0.0f)
                break;
            if (abdfResult->sampledType.isDispersive())
                wls.flags |= WavelengthSamples::LambdaIsSelected;
            Vector3D vecIn = interPt->fromLocal(abdfResult->dirLocal);
            float cosIn = interPt->calcCosTerm(vecIn);
            SampledSpectrum weight = abdfValue * (cosIn / abdfResult->dirPDF);
            
            // Russian roulette
            RRProb = std::min(weight.importance(wlHint), 1.0f);
            if (pathSampler.getPathTerminationSample() < RRProb)
                weight /= RRProb;
            else
                break;
            
            alpha *= weight;
            ray = Ray(interPt->getPosition(), vecIn, ray.time, Ray::Epsilon);
            SLRAssert(weight.allFinite(),
                      "weight: unexpected value detected:\nweight: %s\nfs: %s\nlength: %u, cos: %g, dirPDF: %g",
                      weight.toString().c_str(), abdfValue.toString().c_str(), uint32_t(vertices.size()) - 1, cosIn, abdfResult->dirPDF);
            
            VBPTVertex &vtxNextToLast = vertices[vertices.size() - 2];
            vtxNextToLast.revSpatialPDF = abdfResult->reverseDirPDF() * cosLast / dist2;
            vtxNextToLast.revRRProb = std::min((abdfResult->reverseValue() * cosOut / abdfResult->reverseDirPDF()).importance(wlHint), 1.0f);
            
            cosLast = cosIn;
            dirPDF = abdfResult->dirPDF;
            sampledType = abdfResult->sampledType;
        }
    }
    
    // calculate power heuristic MIS weight
    float VolumetricBPTRenderer::Job::calculateMISWeight(float lExtend1stSpatialPDF, float lExtend1stRRProb, float lExtend2ndSpatialPDF, float lExtend2ndRRProb,
                                                         float eExtend1stSpatialPDF, float eExtend1stRRProb, float eExtend2ndSpatialPDF, float eExtend2ndRRProb,
                                                         uint32_t numLVtx, uint32_t numEVtx) const {
        const auto extendAndShorten = [](float extend1stSpatialPDF, float extend1stRRProb, float extend2ndSpatialPDF, float extend2ndRRProb,
                                         const std::vector<VBPTVertex> &subPathToShorten, uint32_t numVertices, uint32_t minNumVertices,
                                         FloatSum* recMISWeight) {
            if (numVertices > minNumVertices) {
                const VBPTVertex &endVtx = subPathToShorten[numVertices - 1];
                float PDFRatio = extend1stSpatialPDF * extend1stRRProb / (endVtx.spatialPDF * endVtx.RRProb);
                bool shortenIsDeltaSampled = endVtx.sampledType.isDelta();
                if (!shortenIsDeltaSampled)
                    *recMISWeight += PDFRatio * PDFRatio;
                bool prevIsDeltaSampled = shortenIsDeltaSampled;
                if (numVertices - 1 > minNumVertices) {
                    const VBPTVertex &newVtx = subPathToShorten[numVertices - 2];
                    PDFRatio *= extend2ndSpatialPDF * extend2ndRRProb / (newVtx.spatialPDF * newVtx.RRProb);
                    shortenIsDeltaSampled = newVtx.sampledType.isDelta();
                    if (!shortenIsDeltaSampled && !prevIsDeltaSampled)
                        *recMISWeight += PDFRatio * PDFRatio;
                    prevIsDeltaSampled = shortenIsDeltaSampled;
                    for (int i = numVertices - 2; i > minNumVertices; --i) {
                        const VBPTVertex &newVtx = subPathToShorten[i - 1];
                        PDFRatio *= newVtx.revSpatialPDF * newVtx.revRRProb / (newVtx.spatialPDF * newVtx.RRProb);
                        shortenIsDeltaSampled = newVtx.sampledType.isDelta();
                        if (!shortenIsDeltaSampled && !prevIsDeltaSampled)
                            *recMISWeight += PDFRatio * PDFRatio;
                        prevIsDeltaSampled = shortenIsDeltaSampled;
                    }
                }
            }
        };
        
        // initialize the reciprocal of MISWeight by 1. This corresponds to the current strategy (numLVtx, numEVtx).
        FloatSum recMISWeight = 1;
        
        // extend/shorten light/eye subpath, not consider implicit light subpath reaching a lens.
        const uint32_t minEyeVertices = 1;
        extendAndShorten(lExtend1stSpatialPDF, lExtend1stRRProb, lExtend2ndSpatialPDF, lExtend2ndRRProb,
                         eyeVertices, numEVtx, minEyeVertices, &recMISWeight);
        
        // extend/shorten eye/light subpath, consider implicit eye subpath reaching a light.
        const uint32_t minLightVertices = 0;
        extendAndShorten(eExtend1stSpatialPDF, eExtend1stRRProb, eExtend2ndSpatialPDF, eExtend2ndRRProb,
                         lightVertices, numLVtx, minLightVertices, &recMISWeight);
        
        return 1.0f / recMISWeight;
    }
}

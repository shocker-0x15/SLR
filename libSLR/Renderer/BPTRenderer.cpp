//
//  BPTRenderer.cpp
//
//  Created by 渡部 心 on 2016/01/31.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "BPTRenderer.h"

#include "../MemoryAllocators/ArenaAllocator.h"
#include "../Core/random_number_generator.h"
#include "../Core/camera.h"
#include "../Core/light_path_sampler.h"
#include "../Core/ImageSensor.h"
#include "../Core/RenderSettings.h"
#include "../Core/ProgressReporter.h"
#include "../RNG/XORShiftRNG.h"
#include "../Scene/Scene.h"
#include "../Helper/ThreadPool.h"

namespace SLR {
    BPTRenderer::BPTRenderer(uint32_t spp) : m_samplesPerPixel(spp) {
    }
    
    void BPTRenderer::render(const Scene &scene, const RenderSettings &settings) const {
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
        
        printf("Bidirectional Path Tracing: %u[spp]\n", m_samplesPerPixel);
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
    
    void BPTRenderer::Job::kernel(uint32_t threadID) {
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
                wlHint = wls.selectedLambdaIndex;
                eyeVertices.clear();
                lightVertices.clear();
                
                // light subpath generation
                {
                    // select one light from all the lights in the scene.
                    float lightProb;
                    SurfaceLight light;
                    scene->selectSurfaceLight(pathSampler.getLightSelectionSample(), time, &light, &lightProb);
                    SLRAssert(std::isfinite(lightProb), "lightProb: unexpected value detected: %f", lightProb);
                    
                    // sample a ray with its radiance (emittance, EDF value) from the selected light.
                    LightPosQuery lightPosQuery(time, wls);
                    SurfaceLightPosQueryResult lightPosResult;
                    EDFQuery edfQuery;
                    EDFQueryResult edfResult;
                    EDF* edf;
                    SampledSpectrum Le0, Le1;
                    Ray ray;
                    float epsilon;
                    light.sampleRay(lightPosQuery, pathSampler.getSurfaceLightPosSample(), edfQuery, pathSampler.getEDFSample(), mem,
                                    &lightPosResult, &Le0, &edf, &edfResult, &Le1, &ray, &epsilon);
                    
                    // register the first light vertex.
                    float lightAreaPDF = lightProb * lightPosResult.areaPDF;
                    lightVertices.emplace_back(lightPosResult.surfPt, mem.create<EDFProxy>(edf, edfQuery),
                                               Le0 / lightAreaPDF, 0.0f, lightAreaPDF, 1.0f, lightPosResult.posType, false);
                    
                    // create subsequent light subpath vertices by tracing in the scene.
                    SampledSpectrum alpha = lightVertices.back().alpha * Le1 * (lightPosResult.surfPt.calcCosTerm(ray.dir) / edfResult.dirPDF);
                    generateSubPath(wls, alpha, ray, epsilon, edfResult.dirPDF, edfResult.dirType, edfResult.dir_sn.z, true, pathSampler, mem);
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
                    Ray ray;
                    float epsilon;
                    camera->sampleRay(lensQuery, pathSampler.getLensPosSample(), WeSample, mem, 
                                      &lensResult, &We0, &idf, &WeResult, &We1, &ray, &epsilon);
                    
                    // register the first eye vertex.
                    eyeVertices.emplace_back(lensResult.surfPt, mem.create<IDFProxy>(idf),
                                             We0 / (lensResult.areaPDF * selectWLPDF), 0.0f, lensResult.areaPDF, 1.0f, lensResult.posType, false);
                    
                    // create subsequent eye subpath vertices by tracing in the scene.
                    SampledSpectrum alpha = eyeVertices.back().alpha * We1 * (lensResult.surfPt.calcCosTerm(ray.dir) / WeResult.dirPDF);
                    generateSubPath(wls, alpha, ray, epsilon, WeResult.dirPDF, WeResult.dirType, WeResult.dirLocal.z, false, pathSampler, mem);
                }
                
                // connection
                for (int t = 1; t <= eyeVertices.size(); ++t) {
                    const BPTVertex &eVtx = eyeVertices[t - 1];
                    for (int s = 1; s <= lightVertices.size(); ++s) {
                        const BPTVertex &lVtx = lightVertices[s - 1];
                        
                        // ----------------------------------------------------------------
                        // calculate the remaining factors of the full path
                        // that are not included in the precomputed weights.
                        
                        float connectDist2;
                        Vector3D connectionVector = lVtx.surfPt.getDirectionFrom(eVtx.surfPt.getPosition(), &connectDist2);
                        float cosLightEnd = lVtx.surfPt.calcCosTerm(connectionVector);
                        float cosEyeEnd = eVtx.surfPt.calcCosTerm(connectionVector);
                        float G = cosEyeEnd * cosLightEnd / connectDist2;
                        
                        Vector3D lConnectVector = lVtx.surfPt.toLocal(-connectionVector);
                        SampledSpectrum lRevDDF;
                        SampledSpectrum lDDF = lVtx.ddf->evaluate(lConnectVector, &lRevDDF);
                        float eExtend2ndDirPDF;
                        float lExtend1stDirPDF = lVtx.ddf->evaluatePDF(lConnectVector, &eExtend2ndDirPDF);
                        
                        Vector3D eConnectVector = eVtx.surfPt.toLocal(connectionVector);
                        SampledSpectrum eRevDDF;
                        SampledSpectrum eDDF = eVtx.ddf->evaluate(eConnectVector, &eRevDDF);
                        float lExtend2ndDirPDF;
                        float eExtend1stDirPDF = eVtx.ddf->evaluatePDF(eConnectVector, &lExtend2ndDirPDF);
                        
                        SampledSpectrum connectionTerm = lDDF * G * eDDF;
                        if (connectionTerm == SampledSpectrum::Zero)
                            continue;
                        
                        if (!scene->testVisibility(eVtx.surfPt, lVtx.surfPt, time))
                            continue;
                        
                        if (lVtx.lambdaSelected || eVtx.lambdaSelected)
                            connectionTerm[wls.selectedLambdaIndex] *= WavelengthSamples::NumComponents;
                        
                        // ----------------------------------------------------------------
                        
                        
                        
                        // ----------------------------------------------------------------
                        // calculate the 1st and 2nd subpath extending PDFs and probabilities.
                        // They can't be stored in advance because they depend on the connection.
                        
                        float lExtend1stAreaPDF, lExtend1stRRProb, lExtend2ndAreaPDF = 0.0f, lExtend2ndRRProb = 0.0f;
                        {
                            lExtend1stAreaPDF = lExtend1stDirPDF * cosEyeEnd / connectDist2;
                            lExtend1stRRProb = s > 1 ? std::min((lDDF * cosLightEnd / lExtend1stDirPDF).importance(wlHint), 1.0f) : 1.0f;
                            
                            if (t > 1) {
                                BPTVertex &eVtxNextToEnd = eyeVertices[t - 2];
                                float dist2;
                                Vector3D dir2nd = eVtx.surfPt.getDirectionFrom(eVtxNextToEnd.surfPt.getPosition(), &dist2);
                                lExtend2ndAreaPDF = lExtend2ndDirPDF * eVtxNextToEnd.surfPt.calcCosTerm(dir2nd) / dist2;
                                lExtend2ndRRProb = std::min((eRevDDF * eVtx.cosIn / lExtend2ndDirPDF).importance(wlHint), 1.0f);
                            }
                        }
                        float eExtend1stAreaPDF, eExtend1stRRProb, eExtend2ndAreaPDF = 0.0f, eExtend2ndRRProb = 0.0f;
                        {
                            eExtend1stAreaPDF = eExtend1stDirPDF * cosLightEnd / connectDist2;
                            eExtend1stRRProb = t > 1 ? std::min((eDDF * cosEyeEnd / eExtend1stDirPDF).importance(wlHint), 1.0f) : 1.0f;
                            
                            if (s > 1) {
                                BPTVertex &lVtxNextToEnd = lightVertices[s - 2];
                                float dist2;
                                Vector3D dir2nd = lVtxNextToEnd.surfPt.getDirectionFrom(lVtx.surfPt.getPosition(), &dist2);
                                eExtend2ndAreaPDF = eExtend2ndDirPDF * lVtxNextToEnd.surfPt.calcCosTerm(dir2nd) / dist2;
                                eExtend2ndRRProb = std::min((lRevDDF * lVtx.cosIn / eExtend2ndDirPDF).importance(wlHint), 1.0f);
                            }
                        }
                        
                        // ----------------------------------------------------------------
                        
                        
                        
                        // ----------------------------------------------------------------
                        // calculate MIS weight and store weighted contribution to a sensor.
                        
                        float MISWeight = calculateMISWeight(lExtend1stAreaPDF, lExtend1stRRProb, lExtend2ndAreaPDF, lExtend2ndRRProb,
                                                             eExtend1stAreaPDF, eExtend1stRRProb, eExtend2ndAreaPDF, eExtend2ndRRProb, s, t);
                        if (std::isinf(MISWeight) || std::isnan(MISWeight))
                            continue;
                        SLRAssert(MISWeight >= 0 && MISWeight <= 1.0f, "invalid MIS weight: %g", MISWeight);
                        SampledSpectrum contribution = MISWeight * lVtx.alpha * connectionTerm * eVtx.alpha;
                        SLRAssert(contribution.allFinite() && !contribution.hasNegative(),
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
                        
                        // ----------------------------------------------------------------
                    }
                }
                
                mem.reset();
            }
        }
        reporter->update();
    }
    
    void BPTRenderer::Job::generateSubPath(const WavelengthSamples &initWLs, const SampledSpectrum &initAlpha, const Ray &initRay, float initEpsilon, float dirPDF, DirectionType sampledType,
                                           float cosLast, bool adjoint, IndependentLightPathSampler &pathSampler, ArenaAllocator &mem) {
        std::vector<BPTVertex> &vertices = adjoint ? lightVertices : eyeVertices;
        
        // reject invalid values.
        if (dirPDF == 0.0f)
            return;
        
        WavelengthSamples wls = initWLs;
        Ray ray = initRay;
        RaySegment segment(initEpsilon);
        SampledSpectrum alpha = initAlpha;
        
        SurfaceInteraction si;
        SurfacePoint surfPt;
        float RRProb = 1.0f;
        while (scene->intersect(ray, segment, &si)) {
            si.calculateSurfacePoint(&surfPt);
            
            float dist2 = squaredDistance(vertices.back().surfPt, surfPt);
            Vector3D dirOut_sn = surfPt.toLocal(-ray.dir);
            Normal3D gNorm_sn = surfPt.getLocalGeometricNormal();
            float cosOut = surfPt.calcCosTerm(-ray.dir);
            BSDF* bsdf = surfPt.createBSDF(wls, mem);
            BSDFQuery fsQuery(dirOut_sn, gNorm_sn, wls.selectedLambdaIndex, DirectionType::All, true, adjoint);
            
            float areaPDF = dirPDF * cosOut / dist2;
            vertices.emplace_back(surfPt, mem.create<BSDFProxy>(bsdf, fsQuery), alpha, cosOut, 
                                  areaPDF, RRProb, sampledType, wls.wavelengthSelected());
            
            // implicit path (zero light subpath vertices, s = 0)
            if (!adjoint && surfPt.isEmitting()) {
                EDF* edf = surfPt.createEDF(wls, mem);
                SampledSpectrum Le0 = surfPt.emittance(wls);
                SampledSpectrum Le1 = edf->evaluate(EDFQuery(), dirOut_sn);
                
                float extend1stAreaPDF = si.getLightProb() * surfPt.evaluateAreaPDF();
                float extend2ndAreaPDF = edf->evaluatePDF(EDFQuery(), dirOut_sn) * cosLast / dist2;
                
                float MISWeight = calculateMISWeight(extend1stAreaPDF, 1.0f, extend2ndAreaPDF, 1.0f,
                                                     0.0f, 0.0f, 0.0f, 0.0f,
                                                     0, (uint32_t)vertices.size());
                if (!std::isinf(MISWeight) && !std::isnan(MISWeight)) {
                    SampledSpectrum contribution = MISWeight * alpha * Le0 * Le1;
                    SLRAssert(MISWeight >= 0 && MISWeight <= 1.0f, "invalid MIS weight: %g", MISWeight);
                    SLRAssert(contribution.allFinite() && !contribution.hasNegative(),
                              "Unexpected value detected: %s\n"
                              "pix: (%f, %f)", contribution.toString().c_str(), curPx, curPy);
                    if (wls.wavelengthSelected())
                        contribution[wls.selectedLambdaIndex] *= WavelengthSamples::NumComponents;
                    sensor->add(curPx, curPy, wls, contribution);
                }
            }
            
            if (surfPt.atInfinity()) {
                vertices.pop_back();
                break;
            }

            BSDFQueryResult fsResult;
            SampledSpectrum fs = bsdf->sample(fsQuery, pathSampler.getBSDFSample(), &fsResult);
            if (fs == SampledSpectrum::Zero || fsResult.dirPDF == 0.0f)
                break;
            if (fsResult.sampledType.isDispersive() && !wls.wavelengthSelected())
                wls.flags |= WavelengthSamples::WavelengthIsSelected;
            Vector3D vecIn = surfPt.fromLocal(fsResult.dirLocal);
            float cosIn = surfPt.calcCosTerm(vecIn);
            SampledSpectrum weight = fs * (cosIn / fsResult.dirPDF);
            
            // Russian roulette
            RRProb = std::min(weight.importance(wlHint), 1.0f);
            if (pathSampler.getPathTerminationSample() < RRProb)
                weight /= RRProb;
            else
                break;
            
            alpha *= weight;
            ray = Ray(surfPt.getPosition(), vecIn, ray.time);
            segment = RaySegment(Ray::Epsilon);
            SLRAssert(weight.allFinite(),
                      "weight: unexpected value detected:\nweight: %s\nfs: %s\nlength: %u, cos: %g, dirPDF: %g",
                      weight.toString().c_str(), fs.toString().c_str(), uint32_t(vertices.size()) - 1, cosIn, fsResult.dirPDF);
            
            BPTVertex &vtxNextToLast = vertices[vertices.size() - 2];
            vtxNextToLast.revAreaPDF = fsResult.reverse.dirPDF * cosLast / dist2;
            vtxNextToLast.revRRProb = std::min((fsResult.reverse.value * cosOut / fsResult.reverse.dirPDF).importance(wlHint), 1.0f);
            
            cosLast = cosIn;
            dirPDF = fsResult.dirPDF;
            sampledType = fsResult.sampledType;
            si = SurfaceInteraction();
        }
    }
    
    // calculate power heuristic MIS weight
    float BPTRenderer::Job::calculateMISWeight(float lExtend1stAreaPDF, float lExtend1stRRProb, float lExtend2ndAreaPDF, float lExtend2ndRRProb,
                                               float eExtend1stAreaPDF, float eExtend1stRRProb, float eExtend2ndAreaPDF, float eExtend2ndRRProb,
                                               uint32_t numLVtx, uint32_t numEVtx) const {
        const auto extendAndShorten = [](float extend1stAreaPDF, float extend1stRRProb, float extend2ndAreaPDF, float extend2ndRRProb,
                                         const std::vector<BPTVertex> &subPathToShorten, uint32_t numVertices, uint32_t minNumVertices,
                                         FloatSum* recMISWeight) {
            if (numVertices > minNumVertices) {
                const BPTVertex &endVtx = subPathToShorten[numVertices - 1];
                float PDFRatio = extend1stAreaPDF * extend1stRRProb / (endVtx.areaPDF * endVtx.RRProb);
                bool shortenIsDeltaSampled = endVtx.sampledType.isDelta();
                if (!shortenIsDeltaSampled)
                    *recMISWeight += PDFRatio * PDFRatio;
                bool prevIsDeltaSampled = shortenIsDeltaSampled;
                
                if (numVertices - 1 > minNumVertices) {
                    const BPTVertex &newVtx = subPathToShorten[numVertices - 2];
                    PDFRatio *= extend2ndAreaPDF * extend2ndRRProb / (newVtx.areaPDF * newVtx.RRProb);
                    shortenIsDeltaSampled = newVtx.sampledType.isDelta();
                    if (!shortenIsDeltaSampled && !prevIsDeltaSampled)
                        *recMISWeight += PDFRatio * PDFRatio;
                    prevIsDeltaSampled = shortenIsDeltaSampled;
                    
                    for (int i = numVertices - 2; i > minNumVertices; --i) {
                        const BPTVertex &newVtx = subPathToShorten[i - 1];
                        PDFRatio *= newVtx.revAreaPDF * newVtx.revRRProb / (newVtx.areaPDF * newVtx.RRProb);
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
        extendAndShorten(lExtend1stAreaPDF, lExtend1stRRProb, lExtend2ndAreaPDF, lExtend2ndRRProb,
                         eyeVertices, numEVtx, minEyeVertices, &recMISWeight);
        
        // extend/shorten eye/light subpath, consider implicit eye subpath reaching a light.
        const uint32_t minLightVertices = 0;
        extendAndShorten(eExtend1stAreaPDF, eExtend1stRRProb, eExtend2ndAreaPDF, eExtend2ndRRProb,
                         lightVertices, numLVtx, minLightVertices, &recMISWeight);
        
        return 1.0f / recMISWeight;
    }
}

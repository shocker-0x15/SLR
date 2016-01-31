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
        uint32_t numThreads = std::thread::hardware_concurrency();
#else
        uint32_t numThreads = 1;
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
        sensor->addSeparatedBuffers(numThreads);
        
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
        
        for (int s = 0; s < m_samplesPerPixel; ++s) {
#ifdef DEBUG
            ThreadPool threadPool(1);
#else
            ThreadPool threadPool;
#endif
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
        std::vector<BPTVertex> eyeVertices;
        std::vector<BPTVertex> lightVertices;
        for (int ly = 0; ly < numPixelY; ++ly) {
            for (int lx = 0; lx < numPixelX; ++lx) {
                float time = timeStart + rng.getFloat0cTo1o() * (timeEnd - timeStart);
                float px = basePixelX + lx + rng.getFloat0cTo1o();
                float py = basePixelY + ly + rng.getFloat0cTo1o();
                
                float selectWLPDF;
                WavelengthSamples wls = WavelengthSamples::createWithEqualOffsets(rng.getFloat0cTo1o(), rng.getFloat0cTo1o(), &selectWLPDF);
                
                // initialize working area for the current pixel.
                wlHint = wls.selectedLambda;
                eyeVertices.clear();
                lightVertices.clear();
                
                // light subpath generation
                {
                    float lightProb;
                    Light light;
                    scene->selectLight(rng.getFloat0cTo1o(), &light, &lightProb);
                    SLRAssert(!std::isnan(lightProb) && !std::isinf(lightProb), "lightProb: unexpected value detected: %f", lightProb);
                    
                    LightPosQuery lightQuery(time, wls);
                    LightPosSample lightSample(rng.getFloat0cTo1o(), rng.getFloat0cTo1o());
                    LightPosQueryResult lightResult;
                    SampledSpectrum Le0 = light.sample(lightQuery, lightSample, &lightResult);
                    EDF* edf = lightResult.surfPt.createEDF(wls, mem);
                    float lightAreaPDF = lightProb * lightResult.areaPDF;
                    SLRAssert(!std::isnan(lightResult.areaPDF)/* && !std::isinf(lightResult)*/, "areaPDF: unexpected value detected: %f", lightResult.areaPDF);
                    
                    eyeVertices.emplace_back(lightResult.surfPt, mem.create<EDFProxy>(edf), Le0 / lightAreaPDF, lightAreaPDF, 1.0f, WavelengthSamples::Flag(0));
                    
                    EDFQuery edfQuery;
                    EDFSample LeSample(rng.getFloat0cTo1o(), rng.getFloat0cTo1o(), rng.getFloat0cTo1o());
                    EDFQueryResult LeResult;
                    SampledSpectrum Le1 = edf->sample(edfQuery, LeSample, &LeResult);
                    Ray ray(lightResult.surfPt.p, lightResult.surfPt.shadingFrame.fromLocal(LeResult.dir_sn), time);
                    SampledSpectrum alpha = (Le0 * Le1) * (LeResult.dir_sn.z / (lightProb * lightResult.areaPDF * LeResult.dirPDF));
                    generateSubPath(wls, alpha, ray, LeResult.dirPDF, LeResult.dir_sn.z, true, rng, mem);
                }
                
                // eye subpath generation
                {
                    LensPosQuery lensQuery(time, wls);
                    LensPosSample lensSample(rng.getFloat0cTo1o(), rng.getFloat0cTo1o());
                    LensPosQueryResult lensResult;
                    SampledSpectrum We0 = camera->sample(lensQuery, lensSample, &lensResult);
                    IDF* idf = camera->createIDF(lensResult.surfPt, wls, mem);
                    
                    eyeVertices.emplace_back(lensResult.surfPt, mem.create<IDFProxy>(idf), We0 / lensResult.areaPDF, lensResult.areaPDF, 1.0f, WavelengthSamples::Flag(0));
                    
                    IDFSample WeSample(px / imageWidth, py / imageHeight);
                    IDFQueryResult WeResult;
                    SampledSpectrum We1 = idf->sample(WeSample, &WeResult);
                    Ray ray(lensResult.surfPt.p, lensResult.surfPt.shadingFrame.fromLocal(WeResult.dirLocal), time);
                    SampledSpectrum alpha = (We0 * We1) * (WeResult.dirLocal.z / (lensResult.areaPDF * WeResult.dirPDF * selectWLPDF));
                    generateSubPath(wls, alpha, ray, WeResult.dirPDF, WeResult.dirLocal.z, false, rng, mem);
                }
                
                // connection
                for (int t = 1; t <= eyeVertices.size(); ++t) {
                    const BPTVertex &eVtx = eyeVertices[t - 1];
                    
                    // implicit path (zero light subpath vertices, s = 0)
                    if (t > 1 && eVtx.surfPt.isEmitting()) {
                        
                    }
                    
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
                        float cosLightEnd = std::fabs(cVecL.z);

                        Vector3D cVecE = eVtx.surfPt.shadingFrame.toLocal(-connectionVector);
                        DDFQuery queryEyeEnd{eVtx.dirIn_sn, eVtx.gNormal_sn, wlHint, false};
                        SampledSpectrum ddfE = eVtx.ddf->evaluate(queryEyeEnd, cVecE);
                        float cosEyeEnd = std::fabs(cVecE.z);
                        
                        float G = cosEyeEnd * cosLightEnd / dist2;
                        SampledSpectrum connectionTerm = ddfL * G * ddfE;
                        
                        float lExtendDirPDF = lVtx.ddf->evaluatePDF(queryLightEnd, cVecL);
                        float lExtendAreaPDF = lExtendDirPDF * cosEyeEnd / dist2;
                        float lExtendRRProb = std::min((ddfL * cosLightEnd / lExtendDirPDF).luminance(), 1.0f);
                        float eExtendDirPDF = eVtx.ddf->evaluatePDF(queryEyeEnd, cVecE);
                        float eExtendAreaPDF = eExtendDirPDF * cosLightEnd / dist2;
                        float eExtendRRProb = std::min((ddfE * cosEyeEnd / eExtendDirPDF).luminance(), 1.0f);
                        
                        float MISWeight = calculateMISWeight(lExtendAreaPDF, lExtendRRProb, eExtendAreaPDF, eExtendRRProb, s, t);
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
            BSDF* bsdf = surfPt.createBSDF(wls, mem);
            
            float areaPDF = dirPDF * std::fabs(dirOut_sn.z) / (isect.dist * isect.dist);
            vertices.emplace_back(surfPt, mem.create<BSDFProxy>(bsdf), alpha, areaPDF, RRProb, wls.flags);
            
            if (surfPt.atInfinity)
                break;
            
            Normal3D gNorm_sn = surfPt.shadingFrame.toLocal(surfPt.gNormal);
            
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
            cosLast = std::fabs(fsResult.dir_sn.z);
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
            vtxNextToLast.revRRProb = std::min((revInfo.fs * std::fabs(dirOut_sn.z) / revInfo.dirPDF).luminance(), 1.0f);
        }
    }
    
    float BidirectionalPathTracingRenderer::Job::calculateMISWeight(float lExtendAreaPDF, float lExtendRRProb, float eExtendAreaPDF, float eExtendRRProb,
                                                                    uint32_t numLVtx, uint32_t numEVtx) const {
        FloatSum recMISWeight = 1;
        float PDFRatio;
        
        // extend light subpath
        const BPTVertex &eyeEndVtx = eyeVertices[numEVtx - 1];
        PDFRatio = lExtendAreaPDF * lExtendRRProb / (eyeEndVtx.areaPDF * eyeEndVtx.RRProb);
        recMISWeight += PDFRatio * PDFRatio;
        for (int s = numEVtx - 1; s >= 1; --s) {
            const BPTVertex &eyeVtx = eyeVertices[s - 1];
            PDFRatio *= eyeVertices[s].revAreaPDF * eyeVertices[s].revRRProb / (eyeVtx.areaPDF * eyeVtx.RRProb);
            recMISWeight += PDFRatio * PDFRatio;
        }
        
        // extend eye subpath
        const BPTVertex &lightEndVtx = lightVertices[numLVtx - 1];
        PDFRatio = eExtendAreaPDF * eExtendRRProb / (lightEndVtx.areaPDF * lightEndVtx.RRProb);
        recMISWeight += PDFRatio * PDFRatio;
        for (int t = numLVtx; t >= 0; --t) {
            const BPTVertex &lightVtx = eyeVertices[t - 1];
            PDFRatio *= lightVertices[t].revAreaPDF * lightVertices[t].revRRProb / (lightVtx.areaPDF * lightVtx.RRProb);
            recMISWeight += PDFRatio * PDFRatio;
        }
        
        return recMISWeight;
    }
}

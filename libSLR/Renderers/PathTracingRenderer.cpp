//
//  PathTracingRenderer.cpp
//
//  Created by 渡部 心 on 2015/07/11.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "PathTracingRenderer.h"

#include "../Core/RenderSettings.h"
#include "../Helper/ThreadPool.h"
#include "../Core/XORShiftRNG.h"
#include "../Core/ImageSensor.h"
#include "../Core/RandomNumberGenerator.h"
#include "../Memory/ArenaAllocator.h"
#include "../Core/cameras.h"
#include "../Core/geometry.h"
#include "../Core/SurfaceObject.h"
#include "../Core/directional_distribution_functions.h"

namespace SLR {
    PathTracingRenderer::PathTracingRenderer(uint32_t spp) : m_samplesPerPixel(spp) {
        
    }
    
    void PathTracingRenderer::render(const Scene &scene, const RenderSettings &settings) const {
#ifdef DEBUG
        uint32_t numThreads = 1;
#else
        uint32_t numThreads = std::thread::hardware_concurrency();
#endif
        XORShiftRNG topRand(settings.getInt(RenderSettingItem::RNGSeed));
        std::unique_ptr<ArenaAllocator[]> mems = std::unique_ptr<ArenaAllocator[]>(new ArenaAllocator[numThreads]);
        std::unique_ptr<XORShiftRNG[]> rngs = std::unique_ptr<XORShiftRNG[]>(new XORShiftRNG[numThreads]);
        for (int i = 0; i < numThreads; ++i) {
            new (mems.get() + i) ArenaAllocator();
            new (rngs.get() + i) XORShiftRNG(topRand.getUInt());
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
    
    void PathTracingRenderer::Job::kernel(uint32_t threadID) {
        ArenaAllocator &mem = mems[threadID];
        RandomNumberGenerator &rng = *rngs[threadID];
        for (int ly = 0; ly < numPixelY; ++ly) {
            for (int lx = 0; lx < numPixelX; ++lx) {
                float time = timeStart + rng.getFloat0cTo1o() * (timeEnd - timeStart);
                float px = basePixelX + lx + rng.getFloat0cTo1o();
                float py = basePixelY + ly + rng.getFloat0cTo1o();
                
                float selectWLPDF;
                WavelengthSamples wls = WavelengthSamples::createWithEqualOffsets(rng.getFloat0cTo1o(), rng.getFloat0cTo1o(), &selectWLPDF);
                
                LensPosQuery lensQuery(time, wls);
                LensPosSample lensSample(rng.getFloat0cTo1o(), rng.getFloat0cTo1o());
                LensPosQueryResult lensResult;
                SampledSpectrum We0 = camera->sample(lensQuery, lensSample, &lensResult);
                
                IDFSample WeSample(px / imageWidth, py / imageHeight);
                IDFQueryResult WeResult;
                IDF* idf = camera->createIDF(lensResult.surfPt, wls, mem);
                SampledSpectrum We1 = idf->sample(WeSample, &WeResult);
                
                Ray ray(lensResult.surfPt.p, lensResult.surfPt.shadingFrame.fromLocal(WeResult.dirLocal), time);
                SampledSpectrum C = contribution(*scene, wls, ray, rng, mem);
                SLRAssert(C.hasNaN() == false && C.hasInf() == false && C.hasMinus() == false,
                          "Unexpected value detected: %s\n"
                          "pix: (%f, %f)", C.toString().c_str(), px, py);
                
                SampledSpectrum weight = (We0 * We1) * (absDot(ray.dir, lensResult.surfPt.gNormal) / (lensResult.areaPDF * WeResult.dirPDF * selectWLPDF));
                SLRAssert(weight.hasNaN() == false && weight.hasInf() == false && weight.hasMinus() == false,
                          "Unexpected value detected: %s\n"
                          "pix: (%f, %f)", weight.toString().c_str(), px, py);
                sensor->add(px, py, wls, weight * C);
                
                mem.reset();
            }
        }
    }
    
    SampledSpectrum PathTracingRenderer::Job::contribution(const Scene &scene, const WavelengthSamples &initWLs, const Ray &initRay, RandomNumberGenerator &rng, ArenaAllocator &mem) const {
        WavelengthSamples wls = initWLs;
        Ray ray = initRay;
        SurfacePoint surfPt;
        SampledSpectrum alpha = SampledSpectrum::One;
        float initY = alpha.importance(wls.selectedLambda);
        SampledSpectrumSum sp(SampledSpectrum::Zero);
        uint32_t pathLength = 0;
        
        Intersection isect;
        if (!scene.intersect(ray, &isect))
            return SampledSpectrum::Zero;
        isect.getSurfacePoint(&surfPt);
        
        Vector3D dirOut_sn = surfPt.shadingFrame.toLocal(-ray.dir);
        if (surfPt.isEmitting()) {
            EDF* edf = surfPt.createEDF(wls, mem);
            SampledSpectrum Le = surfPt.emittance(wls) * edf->evaluate(EDFQuery(), dirOut_sn);
            sp += alpha * Le;
        }
        if (surfPt.atInfinity)
            return sp;

        while (true) {
            ++pathLength;
            if (pathLength >= 100)
                break;
            Normal3D gNorm_sn = surfPt.shadingFrame.toLocal(surfPt.gNormal);
            BSDF* bsdf = surfPt.createBSDF(wls, mem);
            BSDFQuery fsQuery(dirOut_sn, gNorm_sn, wls.selectedLambda);
            
            // Next Event Estimation (explicit light sampling)
            if (bsdf->hasNonDelta()) {
                float lightProb;
                Light light;
                scene.selectLight(rng.getFloat0cTo1o(), &light, &lightProb);
                SLRAssert(!std::isnan(lightProb) && !std::isinf(lightProb), "lightProb: unexpected value detected: %f", lightProb);
                
                LightPosQuery lpQuery(ray.time, wls);
                LightPosSample lpSample(rng.getFloat0cTo1o(), rng.getFloat0cTo1o());
                LightPosQueryResult lpResult;
                SampledSpectrum M = light.sample(lpQuery, lpSample, &lpResult);
                SLRAssert(!std::isnan(lpResult.areaPDF)/* && !std::isinf(xpResult.areaPDF)*/, "areaPDF: unexpected value detected: %f", lpResult.areaPDF);
                
                if (scene.testVisibility(surfPt, lpResult.surfPt, ray.time)) {
                    float dist2;
                    Vector3D shadowDir = lpResult.surfPt.getDirectionFrom(surfPt.p, &dist2);
                    Vector3D shadowDir_l = lpResult.surfPt.shadingFrame.toLocal(-shadowDir);
                    Vector3D shadowDir_sn = surfPt.shadingFrame.toLocal(shadowDir);
                    
                    EDF* edf = lpResult.surfPt.createEDF(wls, mem);
                    SampledSpectrum Le = M * edf->evaluate(EDFQuery(), shadowDir_l);
                    float lightPDF = lightProb * lpResult.areaPDF;
                    SLRAssert(!Le.hasNaN() && !Le.hasInf(), "Le: unexpected value detected: %s", Le.toString().c_str());
                    
                    SampledSpectrum fs = bsdf->evaluate(fsQuery, shadowDir_sn);
                    float cosLight = absDot(-shadowDir, lpResult.surfPt.gNormal);
                    float bsdfPDF = bsdf->evaluatePDF(fsQuery, shadowDir_sn) * cosLight / dist2;
                    
                    float MISWeight = 1.0f;
                    if (!lpResult.posType.isDelta() && !std::isinf(lpResult.areaPDF))
                        MISWeight = (lightPDF * lightPDF) / (lightPDF * lightPDF + bsdfPDF * bsdfPDF);
                    SLRAssert(MISWeight <= 1.0f, "Invalid MIS weight: %g", MISWeight);
                    
                    float G = absDot(shadowDir_sn, gNorm_sn) * cosLight / dist2;
                    sp += alpha * Le * fs * (G * MISWeight / lightPDF);
                    SLRAssert(!std::isnan(G) && !std::isinf(G), "G: unexpected value detected: %f", G);
                }
            }
            
            // get a next direction by sampling BSDF.
            BSDFSample fsSample(rng.getFloat0cTo1o(), rng.getFloat0cTo1o(), rng.getFloat0cTo1o());
            BSDFQueryResult fsResult;
            SampledSpectrum fs = bsdf->sample(fsQuery, fsSample, &fsResult);
            if (fs == SampledSpectrum::Zero || fsResult.dirPDF == 0.0f)
                break;
            if (fsResult.dirType.isDispersive()) {
                fsResult.dirPDF /= WavelengthSamples::NumComponents;
                wls.flags |= WavelengthSamples::LambdaIsSelected;
            }
            alpha *= fs * absDot(fsResult.dir_sn, gNorm_sn) / fsResult.dirPDF;
            SLRAssert(!alpha.hasInf() && !alpha.hasNaN(),
                      "alpha: %s\nlength: %u, cos: %g, dirPDF: %g",
                      alpha.toString().c_str(), pathLength, absDot(fsResult.dir_sn, gNorm_sn), fsResult.dirPDF);
            
            Vector3D dirIn = surfPt.shadingFrame.fromLocal(fsResult.dir_sn);
            ray = Ray(surfPt.p, dirIn, ray.time, Ray::Epsilon);
            
            // find a next intersection point.
            isect = Intersection();
            if (!scene.intersect(ray, &isect))
                break;
            isect.getSurfacePoint(&surfPt);
            
            dirOut_sn = surfPt.shadingFrame.toLocal(-ray.dir);
            
            // implicit light sampling
            if (surfPt.isEmitting()) {
                float bsdfPDF = fsResult.dirPDF;
                
                EDF* edf = surfPt.createEDF(wls, mem);
                SampledSpectrum Le = surfPt.emittance(wls) * edf->evaluate(EDFQuery(), dirOut_sn);
                float lightProb = scene.evaluateProb(Light(isect.obj));
                float dist2 = surfPt.getSquaredDistance(ray.org);
                float lightPDF = lightProb * surfPt.evaluateAreaPDF() * dist2 / absDot(ray.dir, surfPt.gNormal);
                SLRAssert(!Le.hasNaN() && !Le.hasInf(), "Le: unexpected value detected: %s", Le.toString().c_str());
                SLRAssert(!std::isnan(lightPDF)/* && !std::isinf(lightPDF)*/, "lightPDF: unexpected value detected: %f", lightPDF);
                
                float MISWeight = 1.0f;
                if (!fsResult.dirType.isDelta())
                    MISWeight = (bsdfPDF * bsdfPDF) / (lightPDF * lightPDF + bsdfPDF * bsdfPDF);
                SLRAssert(MISWeight <= 1.0f, "Invalid MIS weight: %g", MISWeight);
                
                sp += alpha * Le * MISWeight;
            }
            if (surfPt.atInfinity)
                break;
            
            // Russian roulette
            float continueProb = std::min(alpha.importance(wls.selectedLambda) / initY, 1.0f);
            if (rng.getFloat0cTo1o() < continueProb)
                alpha /= continueProb;
            else
                break;
        }
        
        return sp;
    }
}

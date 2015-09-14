//
//  PathTracingRenderer.cpp
//
//  Created by 渡部 心 on 2015/07/11.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "PathTracingRenderer.h"

#include "../Core/RenderSettings.h"
#include "../Helper/ThreadPool.h"
#include "../Core/XORShift.h"
#include "../Core/ImageSensor.h"
#include "../Core/RandomNumberGenerator.h"
#include "../Memory/ArenaAllocator.h"
#include "../SceneGraph/Scene.h"
#include "../Core/cameras.h"
#include "../Core/geometry.h"
#include "../Core/SurfaceObject.h"
#include "../Core/directional_distribution_functions.h"
#include "../Helper/StopWatch.h"

PathTracingRenderer::PathTracingRenderer() {
    
}

void PathTracingRenderer::render(const Scene &scene, const RenderSettings &settings) const {
    StopWatch sw;
    
    sw.start();
    
    uint32_t numSamples = settings.getInt(RenderSettingItem::NumSamples);
    uint32_t imageWidth = settings.getInt(RenderSettingItem::ImageWidth);
    uint32_t imageHeight = settings.getInt(RenderSettingItem::ImageHeight);
    ImageSensor sensor(imageWidth, imageHeight);
    
    uint32_t numThreads = std::thread::hardware_concurrency();
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
    
    Job job;
    job.sensor = &sensor;
    job.camera = scene.getMainCamera();
    job.scene = &scene;
    job.imageWidth = sensor.width();
    job.imageHeight = sensor.height();
    job.timeStart = settings.getFloat(RenderSettingItem::TimeStart);
    job.timeEnd = settings.getFloat(RenderSettingItem::TimeEnd);
    job.numPixelX = sensor.tileWidth();
    job.numPixelY = sensor.tileHeight();
    job.mems = mems.get();
    job.rngs = rngRefs.get();
    
    uint32_t imgIdx = 0;
    uint32_t endIdx = 30;
    uint32_t exportTime = 30000;
    
    for (int s = 0; s < numSamples; ++s) {
#ifdef DEBUG
        ThreadPool threadPool(1);
#else
        ThreadPool threadPool;
#endif
        for (int ty = 0; ty < sensor.numTileY(); ++ty) {
            for (int tx = 0; tx < sensor.numTileX(); ++tx) {
                job.basePixelX = tx * sensor.tileWidth();
                job.basePixelY = ty * sensor.tileHeight();
                threadPool.enqueue(std::bind(&Job::kernel, job, std::placeholders::_1));
            }
        }
        threadPool.wait();
        
        if (s == 0) {
            sensor.saveImage("initial.bmp", settings.getFloat(RenderSettingItem::SensorResponse) / (s + 1));
            printf("initial image\n");
        }
        
        uint64_t elapsed;
        if ((elapsed = sw.elapsed()) > exportTime) {
            exportTime += 30000;
            if (exportTime == endIdx * 30000)
                exportTime -= 3000;
            char filename[256];
            sprintf(filename, "%03u.bmp", imgIdx);
            sensor.saveImage(filename, settings.getFloat(RenderSettingItem::SensorResponse) / (s + 1));
            printf("%g sec, %u samples: %s\n", elapsed * 1e-3f, s, filename);
            ++imgIdx;
            if (imgIdx == endIdx)
                break;
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
            
            WavelengthSamples wls(rng.getFloat0cTo1o());
            
            LensPosQuery lensQuery(time, wls);
            LensPosSample lensSample(rng.getFloat0cTo1o(), rng.getFloat0cTo1o());
            LensPosQueryResult lensResult;
            Spectrum We0 = camera->sample(lensQuery, lensSample, &lensResult) * 29.375f;
            
            IDFSample WeSample(px / imageWidth, py / imageHeight);
            IDFQueryResult WeResult;
            IDF* idf = camera->createIDF(lensResult.surfPt, wls, mem);
            Spectrum We1 = idf->sample(WeSample, &WeResult);
            
            Ray ray(lensResult.surfPt.p, lensResult.surfPt.shadingFrame.fromLocal(WeResult.dirLocal), time);
            Spectrum C = contribution(*scene, wls, ray, rng, mem);
            SLRAssert(C.hasNaN() == false && C.hasInf() == false,
                     "Unexpected value detected: %s\n"
                     "pix: (%f, %f)", C.toString().c_str(), px, py);
            
            Spectrum weight = (We0 * We1) / (lensResult.areaPDF * WeResult.dirPDF);
            SLRAssert(weight.hasNaN() == false && weight.hasInf() == false,
                     "Unexpected value detected: %s\n"
                     "pix: (%f, %f)", weight.toString().c_str(), px, py);
            sensor->add(px, py, wls, weight * C);
            
            mem.reset();
        }
    }
}

Spectrum PathTracingRenderer::Job::contribution(const Scene &scene, const WavelengthSamples &wls, const Ray &initRay, RandomNumberGenerator &rng, ArenaAllocator &mem) const {
    Ray ray = initRay;
    SurfacePoint surfPt;
    Spectrum alpha = Spectrum::One;
    float initY = alpha.luminance();
    SpectrumSum sp(Spectrum::Zero);
    uint32_t pathLength = 0;
    
    Intersection isect;
    if (!scene.intersect(ray, &isect))
        return Spectrum::Zero;
    isect.getSurfacePoint(&surfPt);
    
    Vector3D dirOut_sn = surfPt.shadingFrame.toLocal(-ray.dir);
    if (surfPt.isEmitting()) {
        EDF* edf = surfPt.createEDF(wls, mem);
        Spectrum Le = surfPt.emittance(wls) * edf->evaluate(EDFQuery(), dirOut_sn);
        sp += alpha * Le;
    }
    if (surfPt.atInfinity)
        return sp;
    
    while (true) {
        ++pathLength;
        Normal3D gNorm_sn = surfPt.shadingFrame.toLocal(surfPt.gNormal);
        BSDF* bsdf = surfPt.createBSDF(wls, mem);
        
        // Next Event Estimation (explicit light sampling)
        if (bsdf->hasNonDelta()) {
            float lightProb;
            Light light;
            scene.selectLight(rng.getFloat0cTo1o(), &light, &lightProb);
            SLRAssert(!std::isnan(lightProb) && !std::isinf(lightProb), "lightProb: unexpected value detected: %f", lightProb);
            
            LightPosQuery xpQuery(ray.time, wls);
            LightPosSample xpSample(rng.getFloat0cTo1o(), rng.getFloat0cTo1o());
            LightPosQueryResult xpResult;
            Spectrum M = light.sample(xpQuery, xpSample, &xpResult);
            SLRAssert(!std::isnan(xpResult.areaPDF)/* && !std::isinf(xpResult.areaPDF)*/, "areaPDF: unexpected value detected: %f", xpResult.areaPDF);
            
            float dist2;
            Vector3D shadowDir = xpResult.surfPt.getShadowDirection(surfPt, &dist2);
            
            if (scene.testVisiblility(surfPt, xpResult.surfPt, ray.time)) {
                EDF* edf = xpResult.surfPt.createEDF(wls, mem);
                Vector3D shadowDir_l = xpResult.surfPt.shadingFrame.toLocal(-shadowDir);
                Spectrum Le = M * edf->evaluate(EDFQuery(), shadowDir_l);
                float lightPDF = lightProb * xpResult.areaPDF;
                SLRAssert(!Le.hasNaN() && !Le.hasInf(), "Le: unexpected value detected: %s", Le.toString().c_str());
                
                Vector3D shadowDir_sn = surfPt.shadingFrame.toLocal(shadowDir);
                BSDFQuery queryBSDF(dirOut_sn, gNorm_sn);
                Spectrum fs = bsdf->evaluate(queryBSDF, shadowDir_sn);
                float bsdfPDF = bsdf->evaluatePDF(queryBSDF, shadowDir_sn) * std::fabs(shadowDir_l.z) / dist2;
                SLRAssert(!std::isnan(bsdfPDF) && !std::isinf(bsdfPDF), "bsdfPDF: unexpected value detected: %f", bsdfPDF);
                SLRAssert(!fs.hasNaN() && !fs.hasInf(), "fs: unexpected value detected: %s", fs.toString().c_str());
                
                float MISWeight = 1.0f;
                if (!xpResult.isDeltaPos && !std::isinf(xpResult.areaPDF))
                    MISWeight = (lightPDF * lightPDF) / (lightPDF * lightPDF + bsdfPDF * bsdfPDF);
                
                float G = std::fabs(shadowDir_sn.z) * std::fabs(shadowDir_l.z) / dist2;
                sp += alpha * Le * fs * (G * MISWeight / lightPDF);
                SLRAssert(!std::isnan(G) && !std::isinf(G), "G: unexpected value detected: %f", G);
            }
        }
        
        BSDFQuery fsQuery(dirOut_sn, gNorm_sn);
        BSDFSample fsSample(rng.getFloat0cTo1o(), rng.getFloat0cTo1o(), rng.getFloat0cTo1o());
        BSDFQueryResult fsResult;
        Spectrum fs = bsdf->sample(fsQuery, fsSample, &fsResult);
        if (fs == Spectrum::Zero || fsResult.dirPDF == 0.0f)
            break;
        alpha *= fs * (std::fabs(fsResult.dir_sn.z) / fsResult.dirPDF);
        SLRAssert(!alpha.hasInf() && !alpha.hasNaN(), "alpha: unexpected value detected: %s", alpha.toString().c_str());
        
        Vector3D dirIn = surfPt.shadingFrame.fromLocal(fsResult.dir_sn);
        ray = Ray(surfPt.p + Ray::Epsilon * dirIn, dirIn, ray.time);
        
        isect = Intersection();
        if (!scene.intersect(ray, &isect))
            break;
        isect.getSurfacePoint(&surfPt);
        
        dirOut_sn = surfPt.shadingFrame.toLocal(-ray.dir);
        
        // implicit light sampling
        if (surfPt.isEmitting()) {
            float bsdfPDF = fsResult.dirPDF;
            SLRAssert(!std::isnan(bsdfPDF) && !std::isinf(bsdfPDF), "bsdfPDF: unexpected value detected: %f", bsdfPDF);
            
            EDF* edf = surfPt.createEDF(wls, mem);
            Spectrum Le = surfPt.emittance(wls) * edf->evaluate(EDFQuery(), dirOut_sn);
            float lightProb = scene.evaluateProb(Light(isect.obj));
            float dist2 = surfPt.atInfinity ? 1.0f : sqDistance(ray.org, surfPt.p);
            float lightPDF = lightProb * surfPt.evaluateAreaPDF() * dist2 / std::fabs(dirOut_sn.z);
            SLRAssert(!Le.hasNaN() && !Le.hasInf(), "Le: unexpected value detected: %s", Le.toString().c_str());
            SLRAssert(!std::isnan(lightProb) && !std::isinf(lightProb), "lightProb: unexpected value detected: %f", lightProb);
            SLRAssert(!std::isnan(lightPDF)/* && !std::isinf(lightPDF)*/, "lightPDF: unexpected value detected: %f", lightPDF);
            
            float MISWeight = 1.0f;
            if (!fsResult.dirType.isDelta())
                MISWeight = (bsdfPDF * bsdfPDF) / (lightPDF * lightPDF + bsdfPDF * bsdfPDF);
            
            sp += alpha * Le * MISWeight;
        }
        if (surfPt.atInfinity)
            break;
        
        // Russian roulette
        float continueProb = std::min(alpha.luminance() / initY, 1.0f);
        if (rng.getFloat0cTo1o() < continueProb)
            alpha /= continueProb;
        else
            break;
    }
    
    return sp;
}
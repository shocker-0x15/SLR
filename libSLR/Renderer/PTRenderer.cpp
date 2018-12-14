//
//  PTRenderer.cpp
//
//  Created by 渡部 心 on 2015/07/11.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "PTRenderer.h"

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
    PTRenderer::PTRenderer(uint32_t spp) : m_samplesPerPixel(spp) {
        
    }
    
    void PTRenderer::render(const Scene &scene, const RenderSettings &settings) const {
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
        
        slrprintf("Path Tracing: %u[spp]\n", m_samplesPerPixel);
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
                slrprintf("%u samples: %s, %g[s]\n", exportPass, filename, elapsed * 0.001f);
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
    
    void PTRenderer::Job::kernel(uint32_t threadID) {
        ArenaAllocator &mem = mems[threadID];
        IndependentLightPathSampler &pathSampler = pathSamplers[threadID];
        for (int ly = 0; ly < numPixelY; ++ly) {
            for (int lx = 0; lx < numPixelX; ++lx) {
                float time = pathSampler.getTimeSample(timeStart, timeEnd);
                PixelPosition p = pathSampler.getPixelPositionSample(basePixelX + lx, basePixelY + ly);
                
                float selectWLPDF;
                WavelengthSamples wls = WavelengthSamples::createWithEqualOffsets(pathSampler.getWavelengthSample(), pathSampler.getWLSelectionSample(), &selectWLPDF);
                
                LensPosQuery lensQuery(time, wls);
                LensPosQueryResult lensResult;
                SampledSpectrum We0 = camera->sample(lensQuery, pathSampler.getLensPosSample(), &lensResult);
                
                IDFSample WeSample(p.x / imageWidth, p.y / imageHeight);
                IDFQueryResult WeResult;
                IDF* idf = camera->createIDF(lensResult.surfPt, wls, mem);
                SampledSpectrum We1 = idf->sample(WeSample, &WeResult);
                
                Ray ray(lensResult.surfPt.getPosition(), lensResult.surfPt.fromLocal(WeResult.dirLocal), time);
                SampledSpectrum C = contribution(*scene, wls, ray, pathSampler, mem);
                SLRAssert(C.hasNaN() == false && C.hasInf() == false && C.hasNegative() == false,
                          "Unexpected value detected: %s\n"
                          "pix: (%f, %f)", C.toString().c_str(), p.x, p.y);
                
                SampledSpectrum weight = (We0 * We1) * (lensResult.surfPt.calcCosTerm(ray.dir) / (lensResult.areaPDF * WeResult.dirPDF * selectWLPDF));
                SLRAssert(weight.hasNaN() == false && weight.hasInf() == false && weight.hasNegative() == false,
                          "Unexpected value detected: %s\n"
                          "pix: (%f, %f)", weight.toString().c_str(), p.x, p.y);
                
                C *= weight;
                if (C.allFinite())
                    sensor->add(p.x, p.y, wls, C);
                else
                    slrDevPrintf("Unexpected value detected: %s at (%f, %f)\n", C.toString().c_str(), p.x, p.y);
                
                mem.reset();
            }
        }
        reporter->update();
    }
    
    SampledSpectrum PTRenderer::Job::contribution(const Scene &scene, const WavelengthSamples &initWLs, const Ray &initRay, IndependentLightPathSampler &pathSampler, ArenaAllocator &mem) const {
        WavelengthSamples wls = initWLs;
        Ray ray = initRay;
        RaySegment segment;
        SurfacePoint surfPt;
        SampledSpectrum alpha = SampledSpectrum::One;
        float initY = alpha.importance(wls.selectedLambdaIndex);
        SampledSpectrumSum sp(SampledSpectrum::Zero);
        uint32_t pathLength = 0;
        
        SurfaceInteraction si;
        if (!scene.intersect(ray, segment, pathSampler, &si))
            return SampledSpectrum::Zero;
        si.calculateSurfacePoint(&surfPt);
        
        Vector3D dirOut_sn = surfPt.toLocal(-ray.dir);
        if (surfPt.isEmitting()) {
            EDF* edf = surfPt.createEDF(wls, mem);
            SampledSpectrum Le = surfPt.emittance(wls) * edf->evaluate(EDFQuery(), dirOut_sn);
            sp += alpha * Le;
        }
        if (surfPt.atInfinity())
            return sp;

        while (true) {
            ++pathLength;
            if (pathLength >= 100)
                break;
            Normal3D gNorm_sn = surfPt.getLocalGeometricNormal();
            BSDF* bsdf = surfPt.createBSDF(wls, mem);
            BSDFQuery fsQuery(dirOut_sn, gNorm_sn, wls.selectedLambdaIndex);
            
            // Next Event Estimation (explicit light sampling)
            if (bsdf->hasNonDelta()) {
                SurfaceLight light;
                float lightProb;
                scene.selectSurfaceLight(pathSampler.getLightSelectionSample(), ray.time, &light, &lightProb);
                SLRAssert(std::isfinite(lightProb), "lightProb: unexpected value detected: %f", lightProb);
                
                LightPosQuery lpQuery(ray.time, wls);
                SurfaceLightPosQueryResult lpResult;
                SampledSpectrum M = light.sample(lpQuery, pathSampler.getSurfaceLightPosSample(), &lpResult);
                SLRAssert(!std::isnan(lpResult.areaPDF)/* && !std::isinf(xpResult.areaPDF)*/, "areaPDF: unexpected value detected: %f", lpResult.areaPDF);
                
                float fractionalVisibility;
                if (scene.testVisibility(surfPt, lpResult.surfPt, ray.time, &fractionalVisibility)) {
                    float dist2;
                    Vector3D shadowDir = lpResult.surfPt.getDirectionFrom(surfPt.getPosition(), &dist2);
                    Vector3D shadowDir_l = lpResult.surfPt.toLocal(-shadowDir);
                    Vector3D shadowDir_sn = surfPt.toLocal(shadowDir);
                    float cosShading = surfPt.calcCosTerm(shadowDir);
                    float cosLight = lpResult.surfPt.calcCosTerm(shadowDir);
                    
                    if (cosShading * cosLight > 0) {
                        EDF* edf = lpResult.surfPt.createEDF(wls, mem);
                        SampledSpectrum Le = M * edf->evaluate(EDFQuery(), shadowDir_l);
                        float lightPDF = lightProb * lpResult.areaPDF;
                        SLRAssert(Le.allFinite(), "Le: unexpected value detected: %s", Le.toString().c_str());
                        
                        SampledSpectrum fs = bsdf->evaluate(fsQuery, shadowDir_sn);
                        float bsdfPDF = bsdf->evaluatePDF(fsQuery, shadowDir_sn) * cosLight / dist2;
                        
                        float MISWeight = 1.0f;
                        if (!lpResult.posType.isDelta() && !std::isinf(lpResult.areaPDF))
                            MISWeight = (lightPDF * lightPDF) / (lightPDF * lightPDF + bsdfPDF * bsdfPDF);
                        SLRAssert(MISWeight <= 1.0f, "Invalid MIS weight: %g", MISWeight);
                        
                        float G = fractionalVisibility * cosShading * cosLight / dist2;
                        sp += alpha * Le * fs * (G * MISWeight / lightPDF);
                        SLRAssert(std::isfinite(G), "G: unexpected value detected: %f", G);   
                    }
                }
            }
            
            // get a next direction by sampling BSDF.
            BSDFQueryResult fsResult;
            SampledSpectrum fs = bsdf->sample(fsQuery, pathSampler.getBSDFSample(), &fsResult);
            if (fs == SampledSpectrum::Zero || fsResult.dirPDF == 0.0f)
                break;
            if (fsResult.sampledType.isDispersive() && !wls.wavelengthSelected()) {
                fsResult.dirPDF /= WavelengthSamples::NumComponents;
                wls.flags |= WavelengthSamples::WavelengthIsSelected;
            }
            alpha *= fs * absDot(fsResult.dirLocal, gNorm_sn) / fsResult.dirPDF;
            SLRAssert(alpha.allFinite(),
                      "alpha: %s\nlength: %u, cos: %g, dirPDF: %g",
                      alpha.toString().c_str(), pathLength, absDot(fsResult.dirLocal, gNorm_sn), fsResult.dirPDF);
            
            Vector3D dirIn = surfPt.fromLocal(fsResult.dirLocal);
            ray = Ray(surfPt.getPosition(), dirIn, ray.time);
            segment = RaySegment(Ray::Epsilon);
            
            // find a next intersection point.
            si = SurfaceInteraction();
            if (!scene.intersect(ray, segment, pathSampler, &si))
                break;
            si.calculateSurfacePoint(&surfPt);
            
            dirOut_sn = surfPt.toLocal(-ray.dir);
            
            // implicit light sampling
            if (surfPt.isEmitting()) {
                float bsdfPDF = fsResult.dirPDF;
                
                EDF* edf = surfPt.createEDF(wls, mem);
                SampledSpectrum Le = surfPt.emittance(wls) * edf->evaluate(EDFQuery(), dirOut_sn);
                float dist2 = surfPt.getSquaredDistance(ray.org);
                float lightPDF = si.getLightProb() * surfPt.evaluateAreaPDF() * dist2 / surfPt.calcCosTerm(ray.dir);
                SLRAssert(Le.allFinite(), "Le: unexpected value detected: %s", Le.toString().c_str());
                SLRAssert(!std::isnan(lightPDF)/* && !std::isinf(lightPDF)*/, "lightPDF: unexpected value detected: %f", lightPDF);
                
                float MISWeight = 1.0f;
                if (!fsResult.sampledType.isDelta())
                    MISWeight = (bsdfPDF * bsdfPDF) / (lightPDF * lightPDF + bsdfPDF * bsdfPDF);
                SLRAssert(MISWeight <= 1.0f, "Invalid MIS weight: %g", MISWeight);
                
                sp += alpha * Le * MISWeight;
            }
            if (surfPt.atInfinity())
                break;
            
            // Russian roulette
            float continueProb = std::min(alpha.importance(wls.selectedLambdaIndex) / initY, 1.0f);
            if (pathSampler.getPathTerminationSample() < continueProb)
                alpha /= continueProb;
            else
                break;
        }
        
        return sp;
    }
}

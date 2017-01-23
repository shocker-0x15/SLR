//
//  VolumetricPathTracingRenderer.cpp
//
//  Created by 渡部 心 on 2016/09/19.
//  Copyright © 2016年 渡部 心. All rights reserved.
//

#include "VolumetricPathTracingRenderer.h"

#include "../Core/RenderSettings.h"
#include "../Core/ImageSensor.h"
#include "../Core/RandomNumberGenerator.h"
#include "../Core/cameras.h"
#include "../Core/light_path_samplers.h"
#include "../Core/ProgressReporter.h"
#include "../Scene/Scene.h"
#include "../Memory/ArenaAllocator.h"
#include "../Helper/ThreadPool.h"
#include "../RNGs/XORShiftRNG.h"

namespace SLR {
    VolumetricPathTracingRenderer::VolumetricPathTracingRenderer(uint32_t spp) : m_samplesPerPixel(spp) {
        
    }
    
    void VolumetricPathTracingRenderer::render(const Scene &scene, const RenderSettings &settings) const {
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
        
        printf("Volumetric Path Tracing: %u[spp]\n", m_samplesPerPixel);
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
    
    void VolumetricPathTracingRenderer::Job::kernel(uint32_t threadID) {
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
                SLRAssert(C.hasNaN() == false && C.hasInf() == false && C.hasMinus() == false,
                          "Unexpected value detected: %s\n"
                          "pix: (%f, %f)", C.toString().c_str(), p.x, p.y);
                
                SampledSpectrum weight = (We0 * We1) * (lensResult.surfPt.calcCosTerm(ray.dir) / (lensResult.areaPDF * WeResult.dirPDF * selectWLPDF));
                SLRAssert(weight.hasNaN() == false && weight.hasInf() == false && weight.hasMinus() == false,
                          "Unexpected value detected: %s\n"
                          "pix: (%f, %f)", weight.toString().c_str(), p.x, p.y);
                sensor->add(p.x, p.y, wls, weight * C);
                
                mem.reset();
            }
        }
        reporter->update();
    }
    
    SampledSpectrum VolumetricPathTracingRenderer::Job::contribution(const Scene &scene, const WavelengthSamples &initWLs, const Ray &initRay,
                                                                     IndependentLightPathSampler &pathSampler, ArenaAllocator &mem) const {
        WavelengthSamples wls = initWLs;
        Ray ray = initRay;
        SurfacePoint surfPt;
        MediumPoint medPt;
        SampledSpectrum alpha = SampledSpectrum::One;
        float initY = alpha.importance(wls.selectedLambda);
        SampledSpectrumSum sp(SampledSpectrum::Zero);
        uint32_t pathLength = 0;
        
        Interaction* interact;
        SampledSpectrum medThroughput;
        bool singleWavelength;
        InteractionPoint* interPt;
        
        if (!scene.interact(ray, wls, pathSampler, mem, &interact, &medThroughput, &singleWavelength))
            return SampledSpectrum::Zero;
        
        if (singleWavelength && !wls.lambdaSelected()) {
            medThroughput[wls.selectedLambda] *= WavelengthSamples::NumComponents;
            wls.flags |= WavelengthSamples::LambdaIsSelected;
        }
        alpha *= medThroughput;
        
        interPt = interact->createInteractionPoint(mem);
        
        Vector3D dirOut_local = interPt->toLocal(-ray.dir);
        if (interPt->isEmitting()) {
            EDF* edf = interPt->createEDF(wls, mem);
            SampledSpectrum Le = interPt->fluxDensity(wls) * edf->evaluate(EDFQuery(), dirOut_local);
            sp += alpha * Le;
        }
        if (interPt->atInfinity())
            return sp;
        
        while (true) {
            ++pathLength;
            if (pathLength >= 100)
                break;
            
            AbstractBDF* abdf = interPt->createAbstractBDF(wls, mem);
            ABDFQuery* abdfQuery = interPt->createABDFQuery(dirOut_local, wls.selectedLambda, DirectionType::All, false, mem);
            
            // on Surface: probability: 1.0 or probability density delta function
            // in Medium: probability density: extinction coefficient
            alpha *= interPt->evaluateInteractance(wls);
            
            // Next Event Estimation (explicit light sampling)
            if (abdf->hasNonDelta()) {
                Light* light;
                float lightProb;
                scene.selectLight(pathSampler.getLightSelectionSample(), ray.time, mem, &light, &lightProb);
                SLRAssert(std::isfinite(lightProb), "lightProb: unexpected value detected: %f", lightProb);
                
                LightPosQuery lpQuery(ray.time, wls);
                LightPosQueryResult* lpResult;
                SampledSpectrum emittance = light->sample(lpQuery, pathSampler, mem, &lpResult);
                SLRAssert(!std::isnan(lpResult->spatialPDF())/* && !std::isinf(xpResult.spatialPDF)*/,
                          "spatialPDF: unexpected value detected: %f", lpResult->spatialPDF());
                
                InteractionPoint* lightPt = lpResult->getInteractionPoint();
                SampledSpectrum visibility;
                if (scene.testVisibility(interPt, lightPt, ray.time, wls, pathSampler, &visibility, &singleWavelength)) {
                    if (singleWavelength && !wls.lambdaSelected())
                        visibility[wls.selectedLambda] *= WavelengthSamples::NumComponents;
                    
                    float dist2;
                    Vector3D shadowDir = lightPt->getDirectionFrom(interPt->getPosition(), &dist2);
                    Vector3D shadowDir_l = lightPt->toLocal(-shadowDir);
                    Vector3D shadowDir_sn = interPt->toLocal(shadowDir);
                    
                    EDF* edf = lightPt->createEDF(wls, mem);
                    SampledSpectrum Le = emittance * edf->evaluate(EDFQuery(), shadowDir_l);
                    float lightPDF = lightProb * lpResult->spatialPDF();
                    SLRAssert(Le.allFinite(), "Le: unexpected value detected: %s", Le.toString().c_str());
                    
                    SampledSpectrum abdfValue = abdf->evaluate(abdfQuery, shadowDir_sn);
                    float cosLight = lightPt->calcCosTerm(-shadowDir);
                    float abdfPDF = abdf->evaluatePDF(abdfQuery, shadowDir_sn) * cosLight / dist2;
                    
                    float MISWeight = 1.0f;
                    if (!lpResult->sampledPositionType().isDelta() && !std::isinf(lpResult->spatialPDF()))
                        MISWeight = (lightPDF * lightPDF) / (lightPDF * lightPDF + abdfPDF * abdfPDF);
                    SLRAssert(MISWeight <= 1.0f, "Invalid MIS weight: %g", MISWeight);
                    
                    float G = interPt->calcCosTerm(shadowDir) * cosLight / dist2;
                    sp += alpha * visibility * Le * abdfValue * (G * MISWeight / lightPDF);
                    SLRAssert(std::isfinite(G), "G: unexpected value detected: %f", G);
                }
            }
            
            // get a next direction by sampling AbstractBDF.
            ABDFQueryResult* abdfResult;
            SampledSpectrum abdfValue = abdf->sample(abdfQuery, pathSampler, mem, &abdfResult);
            if (abdfValue == SampledSpectrum::Zero || abdfResult->dirPDF == 0.0f)
                break;
            if (abdfResult->sampledType.isDispersive()) {
                abdfResult->dirPDF /= WavelengthSamples::NumComponents;
                wls.flags |= WavelengthSamples::LambdaIsSelected;
            }
            Vector3D dirIn = interPt->fromLocal(abdfResult->dirLocal);
            alpha *= abdfValue * interPt->calcCosTerm(dirIn) / abdfResult->dirPDF;
            SLRAssert(alpha.allFinite(),
                      "alpha: %s\nlength: %u, cos: %g, dirPDF: %g",
                      alpha.toString().c_str(), pathLength, interPt->calcCosTerm(dirIn), abdfResult->dirPDF);
            
            ray = Ray(interPt->getPosition(), dirIn, ray.time, Ray::Epsilon); // No need to adding the epsilon for medium.
            
            // find a next intersection point.
            if (!scene.interact(ray, wls, pathSampler, mem, &interact, &medThroughput, &singleWavelength))
                break;
            
            if (singleWavelength && !wls.lambdaSelected()) {
                medThroughput[wls.selectedLambda] *= WavelengthSamples::NumComponents;
                wls.flags |= WavelengthSamples::LambdaIsSelected;
            }
            alpha *= medThroughput;
            
            interPt = interact->createInteractionPoint(mem);
            dirOut_local = interPt->toLocal(-ray.dir);
            
            // implicit light sampling
            if (interPt->isEmitting()) {
                float abdfPDF = abdfResult->dirPDF;
                
                EDF* edf = interPt->createEDF(wls, mem);
                SampledSpectrum Le = interPt->fluxDensity(wls) * edf->evaluate(EDFQuery(), dirOut_local);
                float dist2 = interPt->getSquaredDistance(ray.org);
                float lightPDF = interact->getLightProb() * interPt->evaluateSpatialPDF() * dist2 / interPt->calcCosTerm(ray.dir);
                SLRAssert(Le.allFinite(), "Le: unexpected value detected: %s", Le.toString().c_str());
                SLRAssert(!std::isnan(lightPDF)/* && !std::isinf(lightPDF)*/, "lightPDF: unexpected value detected: %f", lightPDF);
                
                float MISWeight = 1.0f;
                if (!abdfResult->sampledType.isDelta())
                    MISWeight = (abdfPDF * abdfPDF) / (lightPDF * lightPDF + abdfPDF * abdfPDF);
                SLRAssert(MISWeight <= 1.0f, "Invalid MIS weight: %g", MISWeight);
                
                sp += alpha * Le * MISWeight;
            }
            if (interPt->atInfinity())
                break;
            
            // Russian roulette
            float continueProb = std::min(alpha.importance(wls.selectedLambda) / initY, 1.0f);
            if (pathSampler.getPathTerminationSample() < continueProb)
                alpha /= continueProb;
            else
                break;
        }
        
        return sp;
    }
}

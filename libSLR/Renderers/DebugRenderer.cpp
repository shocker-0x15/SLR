//
//  DebugRenderer.cpp
//
//  Created by 渡部 心 on 2016/02/06.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "DebugRenderer.h"

#include "../Core/RenderSettings.h"
#include "../Helper/ThreadPool.h"
#include "../Memory/ArenaAllocator.h"
#include "../RNGs/XORShiftRNG.h"
#include "../Core/Image.h"
#include "../Core/ImageSensor.h"
#include "../Core/RandomNumberGenerator.h"
#include "../Core/cameras.h"
#include "../Core/geometry.h"
#include "../Core/SurfaceObject.h"
#include "../Core/directional_distribution_functions.h"
#include "../Core/light_path_samplers.h"

namespace SLR {
    DebugRenderer::DebugRenderer(bool channelFlags[(int)ExtraChannel::NumChannels]) {
        for (int i = 0; i < m_channels.size(); ++i)
            m_channels[i] = channelFlags[i];
    }
    
    void DebugRenderer::render(const Scene &scene, const RenderSettings &settings) const {
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
        job.renderer = this;
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
        
        sensor->init(job.imageWidth, job.imageHeight);
        
        DefaultAllocator &defMem = DefaultAllocator::instance();
        std::array<std::unique_ptr<Image2D, Allocator::DeleterType>, (int)ExtraChannel::NumChannels> chImages;
        for (int i = 0; i < m_channels.size(); ++i) {
            if (!m_channels[i]) {
                chImages[i] = nullptr;
                continue;
            }
            switch ((ExtraChannel)i) {
                case ExtraChannel::GeometricNormal:
                    chImages[i] = defMem.createUnique<TiledImage2D>(job.imageWidth, job.imageHeight, ColorFormat::RGB8x3, &defMem);
                    break;
                case ExtraChannel::ShadingNormal:
                    chImages[i] = defMem.createUnique<TiledImage2D>(job.imageWidth, job.imageHeight, ColorFormat::RGB8x3, &defMem);
                    break;
                case ExtraChannel::ShadingTangent:
                    chImages[i] = defMem.createUnique<TiledImage2D>(job.imageWidth, job.imageHeight, ColorFormat::RGB8x3, &defMem);
                    break;
                case ExtraChannel::Distance:
                    chImages[i] = defMem.createUnique<TiledImage2D>(job.imageWidth, job.imageHeight, ColorFormat::Gray8, &defMem);
                    break;
                default:
                    break;
            }
        }
        job.chImages = &chImages;
        
        ThreadPool threadPool(numThreads);
        for (int ty = 0; ty < sensor->numTileY(); ++ty) {
            for (int tx = 0; tx < sensor->numTileX(); ++tx) {
                job.basePixelX = tx * sensor->tileWidth();
                job.basePixelY = ty * sensor->tileHeight();
                threadPool.enqueue(std::bind(&Job::kernel, job, std::placeholders::_1));
            }
        }
        threadPool.wait();
        
//        char filename[256];
//        sprintf(filename, "output.bmp");
//        sensor->saveImage(filename, 1);
        
        for (int i = 0; i < chImages.size(); ++i) {
            if (!chImages[i])
                continue;
            std::string filename;
            switch ((ExtraChannel)i) {
                case ExtraChannel::GeometricNormal:
                    filename = "geometric_normal.bmp";
                    break;
                case ExtraChannel::ShadingNormal:
                    filename = "shading_normal.bmp";
                    break;
                case ExtraChannel::ShadingTangent:
                    filename = "shading_tangent.bmp";
                    break;
                case ExtraChannel::Distance:
                    filename = "distance.bmp";
                    break;
                default:
                    break;
            }
            chImages[i]->saveImage(filename, false);
        }
    }
    
    void DebugRenderer::Job::kernel(uint32_t threadID) {
        ArenaAllocator &mem = mems[threadID];
        IndependentLightPathSampler &pathSampler = *pathSamplers[threadID];
        for (int ly = 0; ly < numPixelY; ++ly) {
            for (int lx = 0; lx < numPixelX; ++lx) {
                float time = pathSampler.getTimeSample(timeStart, timeEnd);
                PixelPosition p = pathSampler.getPixelPositionSample(basePixelX + lx, basePixelY + ly);
                
                float selectWLPDF;
                WavelengthSamples wls = WavelengthSamples::createWithEqualOffsets(pathSampler.getWavelengthSample(), pathSampler.getWLSelectionSample(), &selectWLPDF);
                
                LensPosQuery lensQuery(time, wls);
                LensPosQueryResult lensResult;
                SampledSpectrum We0 = camera->sample(lensQuery, pathSampler.getLensPosSample(), &lensResult);
                (void)We0;
                
                IDFSample WeSample(p.x / imageWidth, p.y / imageHeight);
                IDFQueryResult WeResult;
                IDF* idf = camera->createIDF(lensResult.surfPt, wls, mem);
                SampledSpectrum We1 = idf->sample(WeSample, &WeResult);
                (void)We1;
                
                Ray ray(lensResult.surfPt.p, lensResult.surfPt.shadingFrame.fromLocal(WeResult.dirLocal), time);
                DebugInfo info = contribution(*scene, wls, ray, pathSampler, mem);
                
                for (int i = 0; i < renderer->m_channels.size(); ++i) {
                    if (!renderer->m_channels[i])
                        continue;
                    auto &chImg = (*chImages)[i];
                    switch ((ExtraChannel)i) {
                        case ExtraChannel::GeometricNormal: {
                            RGB8x3 val;
                            val.r = (uint8_t)std::clamp((0.5f * info.surfPt.gNormal.x + 0.5f) * 255, 0.0f, 255.0f);
                            val.g = (uint8_t)std::clamp((0.5f * info.surfPt.gNormal.y + 0.5f) * 255, 0.0f, 255.0f);
                            val.b = (uint8_t)std::clamp((0.5f * info.surfPt.gNormal.z + 0.5f) * 255, 0.0f, 255.0f);
                            chImg->set((int)p.x, (int)p.y, val);
                            break;
                        }
                        case ExtraChannel::ShadingNormal: {
                            RGB8x3 val;
                            val.r = (uint8_t)std::clamp((0.5f * info.surfPt.shadingFrame.z.x + 0.5f) * 255, 0.0f, 255.0f);
                            val.g = (uint8_t)std::clamp((0.5f * info.surfPt.shadingFrame.z.y + 0.5f) * 255, 0.0f, 255.0f);
                            val.b = (uint8_t)std::clamp((0.5f * info.surfPt.shadingFrame.z.z + 0.5f) * 255, 0.0f, 255.0f);
                            chImg->set((int)p.x, (int)p.y, val);
                            break;
                        }
                        case ExtraChannel::ShadingTangent: {
                            RGB8x3 val;
                            val.r = (uint8_t)std::clamp((0.5f * info.surfPt.shadingFrame.x.x + 0.5f) * 255, 0.0f, 255.0f);
                            val.g = (uint8_t)std::clamp((0.5f * info.surfPt.shadingFrame.x.y + 0.5f) * 255, 0.0f, 255.0f);
                            val.b = (uint8_t)std::clamp((0.5f * info.surfPt.shadingFrame.x.z + 0.5f) * 255, 0.0f, 255.0f);
                            chImg->set((int)p.x, (int)p.y, val);
                            break;
                        }
                        case ExtraChannel::Distance: {
                            Gray8 val;
                            SLRAssert_NotImplemented();
                            chImg->set((int)p.x, (int)p.y, val);
                            break;
                        }
                        default:
                            break;
                    }
                }
                
                mem.reset();
            }
        }
    }
    
    DebugRenderer::Job::DebugInfo DebugRenderer::Job::contribution(const Scene &scene, const WavelengthSamples &initWLs, const Ray &initRay,
                                                                   IndependentLightPathSampler &pathSampler, ArenaAllocator &mem) const {
        Ray ray = initRay;
        SurfacePoint surfPt;
        
        DebugInfo ret;
        
        Intersection isect;
        if (!scene.intersect(ray, &isect))
            return DebugInfo();
        isect.getSurfacePoint(&surfPt);
        ret.surfPt = surfPt;
        
        return ret;
    }
}

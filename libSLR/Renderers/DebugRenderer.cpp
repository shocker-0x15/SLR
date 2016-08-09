//
//  DebugRenderer.cpp
//
//  Created by 渡部 心 on 2016/02/06.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "DebugRenderer.h"

#include "../Core/RenderSettings.h"
#include "../Helper/ThreadPool.h"
#include "../Core/XORShiftRNG.h"
#include "../Core/Image.h"
#include "../Core/ImageSensor.h"
#include "../Core/RandomNumberGenerator.h"
#include "../Memory/ArenaAllocator.h"
#include "../Core/cameras.h"
#include "../Core/geometry.h"
#include "../Core/SurfaceObject.h"
#include "../Core/directional_distribution_functions.h"

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
        job.renderer = this;
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
        RandomNumberGenerator &rng = *rngs[threadID];
        for (int ly = 0; ly < numPixelY; ++ly) {
            for (int lx = 0; lx < numPixelX; ++lx) {
                float time = timeStart + rng.getFloat0cTo1o() * (timeEnd - timeStart);
                float px = basePixelX + lx + rng.getFloat0cTo1o();
                float py = basePixelY + ly + rng.getFloat0cTo1o();
                
                float selectWLPDF;
                WavelengthSamples wls = WavelengthSamples::sampleUniform(rng.getFloat0cTo1o(), &selectWLPDF);
                
                LensPosQuery lensQuery(time, wls);
                LensPosSample lensSample(rng.getFloat0cTo1o(), rng.getFloat0cTo1o());
                LensPosQueryResult lensResult;
                SampledSpectrum We0 = camera->sample(lensQuery, lensSample, &lensResult);
                (void)We0;
                
                IDFSample WeSample(px / imageWidth, py / imageHeight);
                IDFQueryResult WeResult;
                IDF* idf = camera->createIDF(lensResult.surfPt, wls, mem);
                SampledSpectrum We1 = idf->sample(WeSample, &WeResult);
                (void)We1;
                
                Ray ray(lensResult.surfPt.p, lensResult.surfPt.shadingFrame.fromLocal(WeResult.dirLocal), time);
                DebugInfo info = contribution(*scene, wls, ray, rng, mem);
                
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
                            chImg->set((int)px, (int)py, val);
                            break;
                        }
                        case ExtraChannel::ShadingNormal: {
                            RGB8x3 val;
                            val.r = (uint8_t)std::clamp((0.5f * info.surfPt.shadingFrame.z.x + 0.5f) * 255, 0.0f, 255.0f);
                            val.g = (uint8_t)std::clamp((0.5f * info.surfPt.shadingFrame.z.y + 0.5f) * 255, 0.0f, 255.0f);
                            val.b = (uint8_t)std::clamp((0.5f * info.surfPt.shadingFrame.z.z + 0.5f) * 255, 0.0f, 255.0f);
                            chImg->set((int)px, (int)py, val);
                            break;
                        }
                        case ExtraChannel::ShadingTangent: {
                            RGB8x3 val;
                            val.r = (uint8_t)std::clamp((0.5f * info.surfPt.shadingFrame.x.x + 0.5f) * 255, 0.0f, 255.0f);
                            val.g = (uint8_t)std::clamp((0.5f * info.surfPt.shadingFrame.x.y + 0.5f) * 255, 0.0f, 255.0f);
                            val.b = (uint8_t)std::clamp((0.5f * info.surfPt.shadingFrame.x.z + 0.5f) * 255, 0.0f, 255.0f);
                            chImg->set((int)px, (int)py, val);
                            break;
                        }
                        case ExtraChannel::Distance: {
                            Gray8 val;
                            SLRAssert_NotImplemented();
                            chImg->set((int)px, (int)py, val);
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
    
    DebugRenderer::Job::DebugInfo DebugRenderer::Job::contribution(const Scene &scene, const WavelengthSamples &initWLs, const Ray &initRay, RandomNumberGenerator &rng, ArenaAllocator &mem) const {
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

//
//  PathTracingRenderer.h
//
//  Created by 渡部 心 on 2015/07/11.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__PathTracingRenderer__
#define __SLR__PathTracingRenderer__

#include "../defines.h"
#include "../references.h"
#include "../Core/Renderer.h"

namespace SLR {
    class SLR_API PathTracingRenderer : public Renderer {
        struct Job {
            const Scene* scene;
            
            ArenaAllocator* mems;
            IndependentLightPathSampler* pathSamplers;
            
            const Camera* camera;
            ImageSensor* sensor;
            float timeStart;
            float timeEnd;
            uint32_t imageWidth;
            uint32_t imageHeight;
            uint32_t numPixelX;
            uint32_t numPixelY;
            uint32_t basePixelX;
            uint32_t basePixelY;
            
            ProgressReporter* reporter;
            
            void kernel(uint32_t threadID);
            SampledSpectrum contribution(const Scene &scene, const WavelengthSamples &initWLs, const Ray &initRay, IndependentLightPathSampler &pathSampler, ArenaAllocator &mem) const;
        };
        
        uint32_t m_samplesPerPixel;
    public:
        PathTracingRenderer(uint32_t spp);
        void render(const Scene &scene, const RenderSettings &settings) const override;
    };    
}

#endif

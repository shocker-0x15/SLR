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

class PathTracingRenderer : public Renderer {
    struct Job {
        ImageSensor* sensor;
        Camera* camera;
        const Scene* scene;
        uint32_t imageWidth;
        uint32_t imageHeight;
        float timeStart;
        float timeEnd;
        ArenaAllocator* mems;
        RandomNumberGenerator** rngs;
        
        uint32_t numPixelX;
        uint32_t numPixelY;
        uint32_t basePixelX;
        uint32_t basePixelY;
        
        void kernel(uint32_t threadID);
        Spectrum contribution(const Scene &scene, const WavelengthSamples &wls, const Ray &initRay, RandomNumberGenerator &rng, ArenaAllocator &mem) const;
    };
public:
    PathTracingRenderer();
    void render(const Scene &scene, const RenderSettings &settings) const override;
};

#endif

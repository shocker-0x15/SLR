//
//
//  Created by 渡部 心 on 2015/05/30.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__cameras__
#define __SLR__cameras__

#include "../defines.h"
#include "../references.h"
#include "geometry.h"

namespace SLR {
    struct SLR_API LensPosQuery {
        float time;
        WavelengthSamples wls;
        LensPosQuery(float t, const WavelengthSamples &lambdas) : time(t), wls(lambdas) { };
    };
    
    struct SLR_API LensPosSample {
        float uPos[2];
        LensPosSample(float up0, float up1) : uPos{up0, up1} { };
    };
    
    struct SLR_API LensPosQueryResult {
        SurfacePoint surfPt;
        float areaPDF;
        bool isDeltaPos;
    };
    
    class SLR_API Camera {
    protected:
        const Transform* m_transform;
    public:
        Camera(const Transform* l2w) : m_transform(l2w) { };
        virtual ~Camera() { };
        
        void setTransform(const Transform* t);
        
        virtual ImageSensor* getSensor() const = 0;
        
        virtual SampledSpectrum sample(const LensPosQuery &query, const LensPosSample &smp, LensPosQueryResult* result) const = 0;
        virtual IDF* createIDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const = 0;
    };    
}

#endif

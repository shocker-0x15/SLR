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
#include "directional_distribution_functions.h"

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
        SampledSpectrum areaPDF;
        DirectionType posType;
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
        virtual Ray sampleRay(const LensPosQuery &lensQuery, const LensPosSample &lensSample, LensPosQueryResult* lensResult, SampledSpectrum* We0, IDF** idf,
                              const IDFSample &WeSample, IDFQueryResult* WeResult, SampledSpectrum* We1,
                              ArenaAllocator &mem) const {
            // sample a position (We0, spatial importance) on the lens surface of the camera.
            *We0 = sample(lensQuery, lensSample, lensResult);
            *idf = createIDF(lensResult->surfPt, lensQuery.wls, mem);

            // sample a direction (directional importance) from IDF, then create subsequent eye subpath vertices by tracing in the scene.
            *We1 = (*idf)->sample(WeSample, WeResult);
            return Ray(lensResult->surfPt.p, lensResult->surfPt.shadingFrame.fromLocal(WeResult->dirLocal), lensQuery.time);
        }
    };    
}

#endif

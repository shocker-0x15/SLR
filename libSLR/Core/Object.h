//
//  Object.h
//
//  Created by 渡部 心 on 2016/09/19.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef __SLR_Object_h__
#define __SLR_Object_h__

#include "../defines.h"
#include "../references.h"
#include "geometry.h"
#include "directional_distribution_functions.h"

namespace SLR {
    class SLR_API Object {
    public:
        virtual bool isEmitting() const = 0;
        virtual float importance() const = 0;
        virtual void selectLight(float u, Light* light, float* prob) const = 0;
        virtual float evaluateProb(const Light &light) const = 0;
        
        virtual SampledSpectrum sample(const Light &light, const LightPosQuery &query, const LightPosSample &smp, LightPosQueryResult* result) const = 0;
        virtual Ray sampleRay(const Light &light,
                              const LightPosQuery &lightPosQuery, const LightPosSample &lightPosSample, LightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                              const EDFQuery &edfQuery, const EDFSample &edfSample, EDFQueryResult* edfResult, SampledSpectrum* Le1,
                              ArenaAllocator &mem) const = 0;
    };
    
    
    struct SLR_API LightPosQuery {
        float time;
        WavelengthSamples wls;
        LightPosQuery(float t, const WavelengthSamples &lambdas) : time(t), wls(lambdas) { }
    };
    
    struct SLR_API LightPosSample {
        float uPos[2];
        LightPosSample(float up0, float up1) : uPos{up0, up1} { }
    };
    
    struct SLR_API LightPosQueryResult {
        SurfacePoint surfPt;
        float areaPDF;
        DirectionType posType;
    };
    
    class SLR_API Light {
        mutable std::vector<const Object*> m_hierarchy;
    public:
        Light() { }
        Light(const std::vector<const Object*> &hierarchy) : m_hierarchy(hierarchy) { }
        Light(const std::vector<const SurfaceObject*> &hierarchy);
        
        void push(const Object* obj) const { m_hierarchy.push_back(obj); }
        void pop() const { m_hierarchy.pop_back(); }
        const Object* top() const { return m_hierarchy.back(); }
        
        SampledSpectrum sample(const LightPosQuery &query, const LightPosSample &smp, LightPosQueryResult* result) const;
        Ray sampleRay(const LightPosQuery &lightPosQuery, const LightPosSample &lightPosSample, LightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                      const EDFQuery &edfQuery, const EDFSample &edfSample, EDFQueryResult* edfResult, SampledSpectrum* Le1,
                      ArenaAllocator &mem) const;
    };
}

#endif /* Object_h */

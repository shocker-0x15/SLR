//
//  Object.cpp
//
//  Created by 渡部 心 on 2016/09/19.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "Object.h"
#include "SurfaceObject.h"

namespace SLR {
    Light::Light(const std::vector<const SurfaceObject*> &hierarchy) {
        m_hierarchy.resize(hierarchy.size());
        std::copy(hierarchy.begin(), hierarchy.end(), m_hierarchy.begin());
    }
    
    SampledSpectrum Light::sample(const LightPosQuery &query, const LightPosSample &smp, LightPosQueryResult* result) const {
        return m_hierarchy.back()->sample(*this, query, smp, result);
    }
    
    Ray Light::sampleRay(const LightPosQuery &lightPosQuery, const LightPosSample &lightPosSample, LightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                         const EDFQuery &edfQuery, const EDFSample &edfSample, EDFQueryResult* edfResult, SampledSpectrum* Le1,
                         ArenaAllocator &mem) const {
        return m_hierarchy.back()->sampleRay(*this, lightPosQuery, lightPosSample, lightPosResult, Le0, edf, edfQuery, edfSample, edfResult, Le1, mem);
    }
}

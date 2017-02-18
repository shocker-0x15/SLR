//
//  IBLEDF.h
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_IBLEDF__
#define __SLR_IBLEDF__

#include "../defines.h"
#include "../references.h"
#include "../Core/directional_distribution_functions.h"

namespace SLR {
    class SLR_API IBLEDF : public EDF {
        float m_worldDiscArea;
    public:
        IBLEDF(float worldDiscArea) : EDF(DirectionType::Emission | DirectionType::LowFreq), m_worldDiscArea(worldDiscArea) { };
        
        SampledSpectrum sample(const EDFQuery &query, const EDFSample &smp, EDFQueryResult* result) const override;
        SampledSpectrum evaluate(const EDFQuery &query, const Vector3D &dir) const override;
        float evaluatePDF(const EDFQuery &query, const Vector3D &dir) const override;
        float weight(const EDFQuery &query) const override { return 1.0f; };
    };
}

#endif /* __SLR_IBLEDF__ */

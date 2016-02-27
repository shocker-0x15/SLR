//
//  IBLEDF.cpp
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "IBLEDF.h"

namespace SLR {
    SampledSpectrum IBLEDF::sample(const EDFQuery &query, const EDFSample &smp, EDFQueryResult* result) const {
        result->dir_sn = Vector3D(0, 0, 1);
        result->dirPDF = 1.0f / m_worldDiscArea;
        result->dirType = m_type;
        return 1.0f / M_PI;
    }
    
    SampledSpectrum IBLEDF::evaluate(const EDFQuery &query, const Vector3D &dir) const {
        if (!query.flags.matches(m_type))
            return SampledSpectrum::Zero;
        return 1.0f / M_PI;
    }
    
    float IBLEDF::evaluatePDF(const EDFQuery &query, const Vector3D &dir) const {
        if (!query.flags.matches(m_type))
            return 0;
        // The true value is: lim_{l to inf} l^2 / (pi * r_w^2) => inf
        return 1.0f / m_worldDiscArea;
    }    
}

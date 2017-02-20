//
//  basic_edfs.h
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_basic_edfs__
#define __SLR_basic_edfs__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/directional_distribution_functions.h"

namespace SLR {
    class SLR_API DiffuseEDF : public EDF {
    public:
        DiffuseEDF() : EDF(DirectionType::Emission | DirectionType::LowFreq) { }
        
        SampledSpectrum sample(const EDFQuery &query, const EDFSample &smp, EDFQueryResult* result) const override;
        SampledSpectrum evaluate(const EDFQuery &query, const Vector3D &dir) const override;
        float evaluatePDF(const EDFQuery &query, const Vector3D &dir) const override;
        float weight(const EDFQuery &query) const override { return 1.0f; };
    };
    
    
    
    class SLR_API IdealDirectionalEDF : public EDF {
        Vector3D m_direction;
    public:
        IdealDirectionalEDF(const Vector3D &dir) : EDF(DirectionType::Emission | DirectionType::Delta0D), m_direction(dir) {
            SLRAssert(m_direction.z > 0, "Z component of the direction must be positive.");
        }
        
        SampledSpectrum sample(const EDFQuery &query, const EDFSample &smp, EDFQueryResult* result) const override;
        SampledSpectrum evaluate(const EDFQuery &query, const Vector3D &dir) const override {
            return SampledSpectrum::Zero;
        }
        float evaluatePDF(const EDFQuery &query, const Vector3D &dir) const override {
            return 0.0f;
        }
        float weight(const EDFQuery &query) const override { return 1.0f; }
    };
}

#endif /* __SLR_basic_edfs__ */

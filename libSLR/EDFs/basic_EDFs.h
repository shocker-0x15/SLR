//
//  basic_EDFs.h
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__basic_EDFs__
#define __SLR__basic_EDFs__

#include "../defines.h"
#include "../references.h"
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
    public:
        IdealDirectionalEDF() : EDF(DirectionType::Emission | DirectionType::Delta0D) { }
        
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

#endif /* defined(__SLR__basic_EDFs__) */

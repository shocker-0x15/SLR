//
//  HomogeneousMedium.cpp
//
//  Created by 渡部 心 on 2016/12/31.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "HomogeneousMedium.h"

namespace SLR {
    bool HomogeneousMedium::interact(const Ray &ray, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                                     MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const {
        SLRAssert_NotImplemented();
        return false;
    }
    
    void HomogeneousMedium::getMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const {
        SLRAssert_NotImplemented();
    }
    
    void HomogeneousMedium::sample(float u0, float u1, float u2, MediumPoint* medPt, float* volumePDF) const {
        SLRAssert_NotImplemented();
    }
    
    float HomogeneousMedium::evaluateVolumePDF(const MediumPoint &medPt) const {
        return 1.0f / volume();
    }
}

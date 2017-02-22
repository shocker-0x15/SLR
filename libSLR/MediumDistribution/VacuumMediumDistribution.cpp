//
//  VacuumMediumDistribution.cpp
//
//  Created by 渡部 心 on 2017/02/22.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "VacuumMediumDistribution.h"

#include "../Core/light_path_sampler.h"

namespace SLR {    
    void VacuumMediumDistribution::calculateMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const {
        SLRAssert_ShouldNotBeCalled();
    }
    
    SampledSpectrum VacuumMediumDistribution::evaluateExtinctionCoefficient(const Point3D &param, const WavelengthSamples &wls) const {
        return SampledSpectrum::Zero;
    }
    
    SampledSpectrum VacuumMediumDistribution::evaluateAlbedo(const Point3D &param, const WavelengthSamples &wls) const {
        return SampledSpectrum::Zero;
    }
    
    void VacuumMediumDistribution::sample(float u0, float u1, float u2, MediumPoint *medPt, float *volumePDF) const {
        SLRAssert_ShouldNotBeCalled();
    }
    
    float VacuumMediumDistribution::evaluateVolumePDF(const MediumPoint &medPt) const {
        return 1.0f / volume();
    }
}

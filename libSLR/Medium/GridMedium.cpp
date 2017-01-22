//
//  GridMedium.cpp
//
//  Created by 渡部 心 on 2017/01/04.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "GridMedium.h"

namespace SLR {
    bool GridMedium::subdivide(Allocator* mem, Medium** fragments, uint32_t* numFragments) const {
        SLRAssert_NotImplemented();
        return true;
    }
    
    bool GridMedium::interact(const Ray &ray, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                              MediumInteraction *mi, SampledSpectrum *medThroughput, bool* singleWavelength) const {
        SLRAssert_NotImplemented();
        return false;
    }
    
    SampledSpectrum GridMedium::extinctionCoefficient(const Point3D &p, const WavelengthSamples &wls) const {
        SLRAssert_NotImplemented();
        return SampledSpectrum::Zero;
    }
    
    void GridMedium::getMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const {
        SLRAssert_NotImplemented();
    }
    
    void GridMedium::queryCoefficients(const Point3D &p, const WavelengthSamples &wls, SampledSpectrum* sigma_s, SampledSpectrum* sigma_e) const {
        SLRAssert_NotImplemented();
    }
    
    void GridMedium::sample(float u0, float u1, float u2, MediumPoint *medPt, float *volumePDF) const {
        SLRAssert_NotImplemented();
    }
    
    float GridMedium::evaluateVolumePDF(const MediumPoint &medPt) const {
        return 1.0f / volume();
    }
    
    
    
    bool SubGridMedium::interact(const Ray &ray, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                                 MediumInteraction *mi, SampledSpectrum *medThroughput, bool* singleWavelength) const {
        SLRAssert_NotImplemented();
        return false;
    }
    
    SampledSpectrum SubGridMedium::extinctionCoefficient(const Point3D &p, const WavelengthSamples &wls) const {
        SLRAssert_NotImplemented();
        return SampledSpectrum::Zero;
    }
    
    void SubGridMedium::getMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const {
        SLRAssert_NotImplemented();
    }
    
    void SubGridMedium::queryCoefficients(const Point3D &p, const WavelengthSamples &wls, SampledSpectrum* sigma_s, SampledSpectrum* sigma_e) const {
        SLRAssert_NotImplemented();
    }
    
    void SubGridMedium::sample(float u0, float u1, float u2, MediumPoint *medPt, float *volumePDF) const {
        SLRAssert_NotImplemented();
    }
    
    float SubGridMedium::evaluateVolumePDF(const MediumPoint &medPt) const {
        return 1.0f / volume();
    }
}

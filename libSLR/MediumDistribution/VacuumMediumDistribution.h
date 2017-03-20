//
//  VacuumMediumDistribution.h
//
//  Created by 渡部 心 on 2017/02/22.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef __SLR_VacuumMediumDistribution__
#define __SLR_VacuumMediumDistribution__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/geometry.h"

namespace SLR {
    class VacuumMediumDistribution : public MediumDistribution {
        BoundingBox3D m_region;
        
    public:
        VacuumMediumDistribution(const BoundingBox3D &region) :
        m_region(region) { }
        
        float majorantExtinctionCoefficientAtWavelength(float wl) const override {
            return 0.0f;
        }
        
        bool subdivide(Allocator* mem, MediumDistribution** fragments, uint32_t* numFragments) const override { return false; }
        
        BoundingBox3D bounds() const override { return m_region; }
        bool contains(const Point3D &p) const override { return m_region.contains(p); }
        bool intersectBoundary(const Ray &ray, const RaySegment &segment, float* distToBoundary, bool* enter) const override {
            return m_region.intersectBoundary(ray, segment, distToBoundary, enter);
        }
        bool interact(const Ray &ray, const RaySegment &segment, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                      MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const override {
            *medThroughput = SampledSpectrum::One;
            *singleWavelength = false;
            return false;
        }
        SampledSpectrum evaluateTransmittance(const Ray &ray, const RaySegment &segment, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                                              bool* singleWavelength) const override {
            *singleWavelength = false;
            return SampledSpectrum::One;
        }
        void calculateMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const override;
        SampledSpectrum evaluateExtinctionCoefficient(const Point3D &param, const WavelengthSamples &wls) const override;
        SampledSpectrum evaluateAlbedo(const Point3D &param, const WavelengthSamples &wls) const override;
        float volume() const override { return m_region.volume(); }
        void sample(float u0, float u1, float u2, MediumPoint* medPt, float* volumePDF) const override;
        float evaluateVolumePDF(const MediumPoint& medPt) const override;
    };
}

#endif /* __SLR_VacuumMediumDistribution__ */

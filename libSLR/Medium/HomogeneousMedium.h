//
//  HomogeneousMedium.h
//
//  Created by 渡部 心 on 2016/12/31.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef __SLR_HomogeneousMedium__
#define __SLR_HomogeneousMedium__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/geometry.h"

namespace SLR {
    class HomogeneousMediumDistribution : public MediumDistribution {
        BoundingBox3D m_region;
        const AssetSpectrum* m_sigma_s;
        const AssetSpectrum* m_sigma_e;
    public:
        HomogeneousMediumDistribution(const BoundingBox3D &region, const AssetSpectrum* sigma_s, const AssetSpectrum* sigma_e) :
        MediumDistribution(sigma_e->calcBounds()), m_region(region), m_sigma_s(sigma_s), m_sigma_e(sigma_e) { }
        
        bool subdivide(Allocator* mem, MediumDistribution** fragments, uint32_t* numFragments) const override { return false; }
        
        BoundingBox3D bounds() const override { return m_region; }
        bool contains(const Point3D &p) const override { return m_region.contains(p); }
        bool intersectBoundary(const Ray &ray, float* distToBoundary, bool* enter) const override {
            return m_region.intersectBoundary(ray, distToBoundary, enter);
        }
        SampledSpectrum extinctionCoefficient(const Point3D &p, const WavelengthSamples &wls) const override {
            if (m_region.contains(p))
                return m_sigma_e->evaluate(wls);
            return SampledSpectrum::Zero;
        }
        bool interact(const Ray &ray, float distanceLimit, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                      MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const override;
        SampledSpectrum evaluateTransmittance(Ray &ray, float distanceLimit, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                                              bool* singleWavelength) const override;
        void getMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const override;
        float volume() const override { return m_region.volume(); }
        void sample(float u0, float u1, float u2, MediumPoint* medPt, float* volumePDF) const override;
        float evaluateVolumePDF(const MediumPoint& medPt) const override;
    };
}

#endif /* __SLR_HomogeneousMedium__ */

//
//  HomogeneousMediumDistribution.h
//
//  Created by 渡部 心 on 2016/12/31.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef __SLR_HomogeneousMediumDistribution__
#define __SLR_HomogeneousMediumDistribution__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/geometry.h"

namespace SLR {
    class HomogeneousMediumDistribution : public MediumDistribution {
        std::array<float, NumStrataForStorage> m_majorantExtinctionCoefficient;
        BoundingBox3D m_region;
        const AssetSpectrum* m_sigma_s;
        const AssetSpectrum* m_sigma_e;
        
    public:
        HomogeneousMediumDistribution(const BoundingBox3D &region, const AssetSpectrum* sigma_s, const AssetSpectrum* sigma_e) :
        m_region(region), m_sigma_s(sigma_s), m_sigma_e(sigma_e) {
            m_sigma_e->calcBounds(NumStrataForStorage, m_majorantExtinctionCoefficient.data());
        }
        
        float majorantExtinctionCoefficientAtWavelength(float wl) const override {
            int index = (wl - WavelengthLowBound) / (WavelengthHighBound - WavelengthLowBound) * NumStrataForStorage;
            index = std::clamp(index, 0, (int)NumStrataForStorage - 1);
            return m_majorantExtinctionCoefficient[index];
        }
        
        bool subdivide(Allocator* mem, MediumDistribution** fragments, uint32_t* numFragments) const override { return false; }
        
        BoundingBox3D bounds() const override { return m_region; }
        bool contains(const Point3D &p) const override { return m_region.contains(p); }
        bool intersectBoundary(const Ray &ray, const RaySegment &segment, float* distToBoundary, bool* enter) const override {
            return m_region.intersectBoundary(ray, segment, distToBoundary, enter);
        }
        bool interact(const Ray &ray, const RaySegment &segment, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                      MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const override;
        SampledSpectrum evaluateTransmittance(const Ray &ray, const RaySegment &segment, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                                              bool* singleWavelength) const override;
        void calculateMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const override;
        SampledSpectrum evaluateExtinctionCoefficient(const Point3D &param, const WavelengthSamples &wls) const override;
        SampledSpectrum evaluateAlbedo(const Point3D &param, const WavelengthSamples &wls) const override;
        float volume() const override { return m_region.volume(); }
        void sample(float u0, float u1, float u2, MediumPoint* medPt, float* volumePDF) const override;
        float evaluateVolumePDF(const MediumPoint& medPt) const override;
    };
}

#endif /* __SLR_HomogeneousMediumDistribution__ */

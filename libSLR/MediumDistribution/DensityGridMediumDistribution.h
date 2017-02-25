//
//  DensityGridMediumDistribution.h
//
//  Created by 渡部 心 on 2017/02/14.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef __SLR_DensityGridMediumDistribution__
#define __SLR_DensityGridMediumDistribution__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/geometry.h"

namespace SLR {
    class DensityGridMediumDistribution : public MediumDistribution {
        BoundingBox3D m_region;
        const AssetSpectrum* m_base_sigma_s;
        const AssetSpectrum* m_base_sigma_e;
        const float* m_density_grid;
        uint32_t m_numX, m_numY, m_numZ;
        
        float calcDensity(const Point3D &param) const;
    public:
        DensityGridMediumDistribution(const BoundingBox3D &region, const AssetSpectrum* base_sigma_s, const AssetSpectrum* base_sigma_e, const float* density_grid,
                                      uint32_t numX, uint32_t numY, uint32_t numZ, float maxExtinctionCoefficient) :
        MediumDistribution(maxExtinctionCoefficient),
        m_region(region), m_base_sigma_s(base_sigma_s), m_base_sigma_e(base_sigma_e), m_density_grid(density_grid), 
        m_numX(numX), m_numY(numY), m_numZ(numZ) { }
        
        bool subdivide(Allocator* mem, MediumDistribution** fragments, uint32_t* numFragments) const override;
        
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

#endif /* __SLR_DensityGridMediumDistribution__ */

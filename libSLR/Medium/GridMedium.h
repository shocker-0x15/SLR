//
//  GridMedium.h
//
//  Created by 渡部 心 on 2017/01/04.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef __SLR_GridMedium__
#define __SLR_GridMedium__

#include "../defines.h"
#include "../references.h"
#include "../Core/geometry.h"

namespace SLR {
    class GridMediumDistribution : public MediumDistribution {
        BoundingBox3D m_region;
        const AssetSpectrum** m_sigma_s_grid;
        const AssetSpectrum** m_sigma_e_grid;
        uint32_t m_numX, m_numY, m_numZ;
        
        friend class SubGridMedium;
    public:
        GridMediumDistribution(const BoundingBox3D &region, const AssetSpectrum** sigma_s_grid, const AssetSpectrum** sigma_e_grid,
                               uint32_t numX, uint32_t numY, uint32_t numZ, float maxExtinctionCoefficient) :
        MediumDistribution(maxExtinctionCoefficient),
        m_region(region), m_sigma_s_grid(sigma_s_grid), m_sigma_e_grid(sigma_e_grid),
        m_numX(numX), m_numY(numY), m_numZ(numZ) { }
        
        bool subdivide(Allocator* mem, MediumDistribution** fragments, uint32_t* numFragments) const override;
        
        BoundingBox3D bounds() const override { return m_region; }
        bool contains(const Point3D &p) const override { return m_region.contains(p); }
        bool intersectBoundary(const Ray &ray, float* distToBoundary, bool* enter) const override {
            return m_region.intersectBoundary(ray, distToBoundary, enter);
        }
        SampledSpectrum extinctionCoefficient(const Point3D &p, const WavelengthSamples &wls) const override;
        bool interact(const Ray &ray, float distanceLimit, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                      MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const override;
        SampledSpectrum evaluateTransmittance(Ray &ray, float distanceLimit, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                                              bool* singleWavelength) const override;
        void getMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const override;
        float volume() const override { return m_region.volume(); }
        void sample(float u0, float u1, float u2, MediumPoint* medPt, float* volumePDF) const override;
        float evaluateVolumePDF(const MediumPoint& medPt) const override;
    };
    
    
    
    class SubGridMediumDistribution : public MediumDistribution {
        BoundingBox3D m_region;
        const GridMediumDistribution* m_entity;
    public:
        SubGridMediumDistribution(const BoundingBox3D &region, const GridMediumDistribution* entity, float maxExtinctionCoefficient) :
        MediumDistribution(maxExtinctionCoefficient), m_region(region), m_entity(entity) { }
        
        bool subdivide(Allocator* mem, MediumDistribution** fragments, uint32_t* numFragments) const override { return false; }
        
        BoundingBox3D bounds() const override { return m_region; }
        bool contains(const Point3D &p) const override { return m_region.contains(p); }
        bool intersectBoundary(const Ray &ray, float* distToBoundary, bool* enter) const override {
            return m_region.intersectBoundary(ray, distToBoundary, enter);
        }
        SampledSpectrum extinctionCoefficient(const Point3D &p, const WavelengthSamples &wls) const override;
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

#endif /* __SLR_GridMedium__ */

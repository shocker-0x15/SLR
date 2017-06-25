//
//  CloudMediumDistribution.h
//
//  Created by 渡部 心 on 2017/03/24.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef __SLR_CloudMediumDistribution__
#define __SLR_CloudMediumDistribution__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/geometry.h"
#include "../Core/distributions.h"

namespace SLR {
    class CloudMediumDistribution : public MediumDistribution {
        std::array<float, NumStrataForStorage> m_majorantExtinctionCoefficient;
        BoundingBox3D m_region;
        const AssetSpectrum* m_base_sigma_e;
        const AssetSpectrum* m_albedo;
        MultiOctaveImprovedPerlinNoise3DGenerator<float> m_distributionGenerator;
        MultiOctaveImprovedPerlinNoise3DGenerator<float> m_densityGenerator;
        
        void saveToFile(const char* fileName, uint32_t resX, uint32_t resY, uint32_t resZ) const;
        float calcDensity(const Point3D &param) const;
    public:
        CloudMediumDistribution(const BoundingBox3D &region, float featureScale, float density, uint32_t rngSeed) : 
        m_region(region),
        m_distributionGenerator(10, 1.0f / featureScale, 1.0f, true, 2.0f, 0.5f, -1),
        m_densityGenerator(3, 1.0f / featureScale, 10.0f * density / featureScale, false, 2.0f, 0.5f, -1) {
            float sigma_e_values[] = {0.1f, 0.1f};
            m_base_sigma_e = new RegularContinuousSpectrum(WavelengthLowBound, WavelengthHighBound, sigma_e_values, 2);
            float albedo_values[] = {0.9f, 0.9f};
            m_albedo = new RegularContinuousSpectrum(WavelengthLowBound, WavelengthHighBound, albedo_values, 2);
            
            m_base_sigma_e->calcBounds(NumStrataForStorage, m_majorantExtinctionCoefficient.data());
            for (int i = 0; i < NumStrataForStorage; ++i)
                m_majorantExtinctionCoefficient[i] *= 200.0f;//m_densityGenerator.getSupValue();
//            saveToFile("Cloud0_16x16x16.vdg", 16, 16, 16);
//            saveToFile("Cloud0_64x64x64.vdg", 64, 64, 64);
//            saveToFile("Cloud0_256x256x256.vdg", 256, 256, 256);
        }
        ~CloudMediumDistribution() {
            delete m_albedo;
            delete m_base_sigma_e;
        }
        
        float majorantExtinctionCoefficientAtWavelength(float wl) const override {
            int index = (wl - WavelengthLowBound) / (WavelengthHighBound - WavelengthLowBound) * NumStrataForStorage;
            index = std::clamp(index, 0, (int)NumStrataForStorage - 1);
            return m_majorantExtinctionCoefficient[index];
        }
        
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

#endif /* __SLR_CloudMediumDistribution__ */

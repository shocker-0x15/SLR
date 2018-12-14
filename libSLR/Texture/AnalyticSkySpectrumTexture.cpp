//
//  AnalyticSkySpectrumTexture.cpp
//
//  Created by 渡部 心 on 2017/04/24.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "AnalyticSkySpectrumTexture.h"

#include "../Core/distributions.h"

namespace SLR {
    class AnalyticSkySpectrumTexture::SunDiscContinuousDistribution2D : public ContinuousDistribution2D {
        ReferenceFrame m_sunDirFrame;
        float m_cosSunRadius;
    public:
        SunDiscContinuousDistribution2D(const Vector3D &sunDirection, float sunRadius) : 
        m_sunDirFrame(sunDirection), m_cosSunRadius(std::cos(sunRadius)) { }
        
        void sample(float u0, float u1, float* d0, float* d1, float* PDF) const override {
            SLRAssert(u0 >= 0 && u0 < 1, "\"u0\" must be in range [0, 1).");
            SLRAssert(u1 >= 0 && u1 < 1, "\"u1\" must be in range [0, 1).");
            Vector3D dir = uniformSampleCone(u0, u1, m_cosSunRadius);
            dir = m_sunDirFrame.fromLocal(dir);
            float theta, phi;
            dir.toPolarYUp(&theta, &phi);
            
            *d0 = phi / (2 * M_PI);
            *d1 = theta / M_PI;
            
            // convert PDF w.r.t solid angle to texture space.
            *PDF = 1.0f / (2 * M_PI * (1 - m_cosSunRadius)) * (2 * M_PI * M_PI * std::sin(theta));
        }
        
        float evaluatePDF(float d0, float d1) const override {
            float phi = d0 * 2 * M_PI;
            float theta = d1 * M_PI;
            Vector3D dir = Vector3D::fromPolarYUp(phi, theta);
            dir = m_sunDirFrame.toLocal(dir);
            
            if (dir.z >= m_cosSunRadius)
                return 1.0f / (2 * M_PI * (1 - m_cosSunRadius)) * (2 * M_PI * M_PI * std::sin(theta));
            else
                return 0.0f;
        }
    };
    
#ifdef SLR_Use_Spectral_Representation
    const uint32_t AnalyticSkySpectrumTexture::NumChannels = 11;
    const float AnalyticSkySpectrumTexture::SampledWavelengths[11] = {
        320.0f, 360.0f, 400.0f, 440.0f, 
        480.0f, 520.0f, 560.0f, 600.0f, 
        640.0f, 680.0f, 720.0f
    };
#else
    const uint32_t AnalyticSkySpectrumTexture::NumChannels = 3;
    // Scale factor to match the result of the spectral version, adjusted by hand. No theoretical basis.
    const float AnalyticSkySpectrumTexture::RadianceScale = 0.010828553542f;
#endif
    
    
    
    AnalyticSkySpectrumTexture::AnalyticSkySpectrumTexture(float solarRadius, float solarElevation, float turbidity, const AssetSpectrum* groundAlbedo, float extAngleOfHorizon, 
                                                           const Texture2DMapping* mapping) :
    m_solarRadius(std::min(solarRadius, (float)M_PI / 2)), m_solarElevation(solarElevation), m_turbidity(turbidity), m_groundAlbedo(groundAlbedo), m_extAngleOfHorizon(extAngleOfHorizon), 
    m_mapping(mapping), 
    m_distribution(nullptr) {
#ifdef SLR_Use_Spectral_Representation
        float albedo[NumChannels];
        groundAlbedo->evaluate(SampledWavelengths, NumChannels, albedo);
        for (int i = 0; i < NumChannels; ++i) {
            m_skyModelStates[i] = arhosekskymodelstate_alloc_init(solarElevation, turbidity, albedo[i]);
            m_skyModelStates[i]->solar_radius = m_solarRadius;
        }
#else
        float albedo[NumChannels];
        groundAlbedo->getRGB(albedo);
        for (int i = 0; i < NumChannels; ++i) {
            m_skyModelStates[i] = arhosek_rgb_skymodelstate_alloc_init(turbidity, albedo[i], solarElevation);
            m_skyModelStates[i]->solar_radius = m_solarRadius;
        }
#endif
        m_sunDirection = Vector3D::fromPolarYUp(M_PI, M_PI / 2 - m_solarElevation);
    }
    
    AnalyticSkySpectrumTexture::~AnalyticSkySpectrumTexture() {
        if (m_distribution) {
#ifdef SLR_Use_Spectral_Representation
            delete m_distribution;
            delete m_skyDomeDistribution;
            delete m_sunDiscDistribution;
#else
            delete m_distribution;
#endif
        }
        
        for (int i = 0; i < NumChannels; ++i)
            arhosekskymodelstate_free(m_skyModelStates[i]);
    }
    
    SampledSpectrum AnalyticSkySpectrumTexture::evaluate(const Point3D &p, const WavelengthSamples &wls) const {
        float theta = M_PI * p.y;
        float mappedTheta = theta / (1 + m_extAngleOfHorizon / (M_PI / 2));
        if (mappedTheta >= M_PI / 2)
            return SampledSpectrum::Zero;
        Vector3D viewVec = Vector3D::fromPolarYUp(2 * M_PI * p.x, theta);
        float gamma = std::acos(std::clamp(dot(viewVec, m_sunDirection), -1.0f, 1.0f));
        
#ifdef SLR_Use_Spectral_Representation
        float sampledValues[NumChannels];
        if (gamma < m_skyModelStates[0]->solar_radius) {
            for (int i = 0; i < NumChannels; ++i)
                sampledValues[i] = arhosekskymodel_solar_radiance(m_skyModelStates[i], theta, gamma, SampledWavelengths[i]);
        }
        else {
            for (int i = 0; i < NumChannels; ++i)
                sampledValues[i] = arhosekskymodel_radiance(m_skyModelStates[i], mappedTheta, gamma, SampledWavelengths[i]);
        }
        
        RegularContinuousSpectrum spectrum(SampledWavelengths[0], SampledWavelengths[NumChannels - 1], sampledValues, NumChannels);
#else
        SampledSpectrum spectrum;
        for (int i = 0; i < NumChannels; ++i)
            spectrum[i] = RadianceScale * arhosek_tristim_skymodel_radiance(m_skyModelStates[i], mappedTheta, gamma, i);
#endif
        SampledSpectrum ret = spectrum.evaluate(wls);
        ret = max(ret, 0.0f);
        SLRAssert(ret.allFinite(), "Invalid vaue.: %s", ret.toString().c_str());
        
        return ret;
    }
    
    const ContinuousDistribution2D* AnalyticSkySpectrumTexture::createIBLImportanceMap() const {
        if (m_distribution) {
            delete m_distribution;
        }
        
        // JP: 天空光の輝度分布と合計エネルギーを計算する。
        // EN: calculate the luminance distribution of the sky dome and its total energy.
        const uint32_t mapWidth = 1024;
        const uint32_t mapHeight = 512;
        FloatSum accSkyDomeEnergy = 0.0f;
        std::function<float(uint32_t, uint32_t)> pickFunc = [this, &mapWidth, &mapHeight, &accSkyDomeEnergy](uint32_t x, uint32_t y) -> float {
            float theta = M_PI * (y + 0.5f) / mapHeight;
            float mappedTheta = theta / (1 + m_extAngleOfHorizon / (M_PI / 2));
            if (mappedTheta >= M_PI / 2)
                return 0.0f;
            Vector3D viewVec = Vector3D::fromPolarYUp(2 * M_PI * (x + 0.5f) / mapWidth, mappedTheta);
            float gamma = std::acos(std::clamp(dot(viewVec, m_sunDirection), -1.0f, 1.0f));
            
#ifdef SLR_Use_Spectral_Representation
            float sampledValues[NumChannels];
            for (int i = 0; i < NumChannels; ++i)
                sampledValues[i] = arhosekskymodel_radiance(m_skyModelStates[i], mappedTheta, gamma, SampledWavelengths[i]);
            
            RegularContinuousSpectrum spectrum(SampledWavelengths[0], SampledWavelengths[NumChannels - 1], sampledValues, NumChannels);
            
            SpectrumStorage yStorage;
            const uint32_t NumDetailSampling = 5;
            for (int i = 0; i < NumDetailSampling; ++i) {
                float wlPDF;
                WavelengthSamples wls = WavelengthSamples::createWithEqualOffsets(0.5f, (float)i / NumDetailSampling, &wlPDF);
                yStorage.add(wls, spectrum.evaluate(wls) / wlPDF);
            }
            
            float luminance = yStorage.getValue().result.luminance() / NumDetailSampling;
#else
            SampledSpectrum spectrum;
            for (int i = 0; i < NumChannels; ++i)
                spectrum[i] = RadianceScale * arhosek_tristim_skymodel_radiance(m_skyModelStates[i], mappedTheta, gamma, i);
            
            float luminance = spectrum.luminance();
#endif
            SLRAssert(std::isfinite(luminance), "Invalid area average value.");
            float contribution = std::sin(M_PI * (y + 0.5f) / mapHeight) * luminance;
            accSkyDomeEnergy += contribution;
            
            return contribution;
        };
     
        m_skyDomeDistribution = new RegularConstantContinuousDistribution2D(mapWidth, mapHeight, pickFunc);
//        m_skyDomeDistribution->exportBMP("distribution.bmp", true);
        
#ifdef SLR_Use_Spectral_Representation
        // JP: 太陽のディスクのエネルギーを計算する。
        // EN: calculate energy from the sun disc.
        const uint32_t NumRadiusSamples = 5;
        const uint32_t NumAngularSamples = 12;
        ReferenceFrame sunDirFrame(m_sunDirection);
        FloatSum sunDiscEnergy = 0.0f;
        for (int ir = 0; ir < NumRadiusSamples; ++ir) {
            float r = m_solarRadius * std::sqrt((ir + 0.5f) / NumRadiusSamples);
            float rLow = m_solarRadius * std::sqrt((float)ir / NumRadiusSamples);
            float rHigh = m_solarRadius * std::sqrt((float)(ir + 1) / NumRadiusSamples);
            float rWidth = rHigh - rLow;
            for (int ia = 0; ia < NumAngularSamples; ++ia) {
                float a = 2 * M_PI * (ia + 0.5f) / NumAngularSamples;
                float solidAngle = rWidth * (2 * M_PI / NumAngularSamples) * 0.5f * (rHigh + rLow);
                
                Vector3D viewVecLocal = Vector3D::fromPolarZUp(a, r); 
                Vector3D viewVec = sunDirFrame.fromLocal(viewVecLocal);
                float theta = std::acos(std::clamp(viewVec.y, -1.0f, 1.0f));
                float gamma = std::acos(std::clamp(dot(viewVec, m_sunDirection), -1.0f, 1.0f));
                SLRAssert(gamma < m_solarRadius, "Invalid.");
                
                float sampledValues[NumChannels];
                for (int i = 0; i < NumChannels; ++i)
                    sampledValues[i] = arhosekskymodel_direct_radiance(m_skyModelStates[i], theta, gamma, SampledWavelengths[i]);
                
                RegularContinuousSpectrum spectrum(SampledWavelengths[0], SampledWavelengths[NumChannels - 1], sampledValues, NumChannels);
                
                SpectrumStorage yStorage;
                const uint32_t NumWLDetailSampling = 5;
                for (int i = 0; i < NumWLDetailSampling; ++i) {
                    float wlPDF;
                    WavelengthSamples wls = WavelengthSamples::createWithEqualOffsets(0.5f, (float)i / NumWLDetailSampling, &wlPDF);
                    yStorage.add(wls, spectrum.evaluate(wls) / wlPDF);
                }
                
                float luminance = yStorage.getValue().result.luminance() / NumWLDetailSampling;
                
                sunDiscEnergy += solidAngle * luminance; 
            }
        }
        float skyDomeEnergy = accSkyDomeEnergy * (2 * M_PI / mapWidth * M_PI / mapHeight);
        m_sunDiscDistribution = new SunDiscContinuousDistribution2D(m_sunDirection, m_solarRadius);
        
        std::array<const ContinuousDistribution2D*, 2> dists{m_skyDomeDistribution, m_sunDiscDistribution};
        std::array<float, 2> importances{skyDomeEnergy, sunDiscEnergy};
        m_distribution = new MultiContinuousDistribution2D(dists.data(), importances.data(), dists.size());
#else
        slrprintf("WARNING: Sun disc is not supported in RGB rendering.");
        m_distribution = m_skyDomeDistribution;
#endif
        
        return m_distribution;
    }
}

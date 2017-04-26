//
//  AnalyticSkySpectrumTexture.cpp
//
//  Created by 渡部 心 on 2017/04/24.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "AnalyticSkySpectrumTexture.h"

#include "../Core/distributions.h"

namespace SLR {
#ifdef Use_Spectral_Representation
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
    
    AnalyticSkySpectrumTexture::AnalyticSkySpectrumTexture(float solarElevation, float turbidity, const AssetSpectrum* groundAlbedo, const Texture2DMapping* mapping) :
    m_solarElevation(solarElevation), m_turbidity(turbidity), m_groundAlbedo(groundAlbedo), m_mapping(mapping) {
#ifdef Use_Spectral_Representation
        float albedo[NumChannels];
        groundAlbedo->evaluate(SampledWavelengths, NumChannels, albedo);
        for (int i = 0; i < NumChannels; ++i)
            m_skyModelStates[i] = arhosekskymodelstate_alloc_init(solarElevation, turbidity, albedo[i]);
#else
        float albedo[NumChannels];
        groundAlbedo->getRGB(albedo);
        for (int i = 0; i < NumChannels; ++i)
            m_skyModelStates[i] = arhosek_rgb_skymodelstate_alloc_init(turbidity, albedo[i], solarElevation);
#endif
    }
    
    AnalyticSkySpectrumTexture::~AnalyticSkySpectrumTexture() {
        for (int i = 0; i < NumChannels; ++i)
            arhosekskymodelstate_free(m_skyModelStates[i]);
    }
    
    SampledSpectrum AnalyticSkySpectrumTexture::evaluate(const Point3D &p, const WavelengthSamples &wls) const {
        float theta = M_PI * p.y;
        if (theta >= M_PI / 2)
            return SampledSpectrum::Zero;
        Vector3D viewVec = Vector3D::fromPolarYUp(2 * M_PI * p.x, theta);
        Vector3D solarVec = Vector3D::fromPolarYUp(M_PI, M_PI / 2 - m_solarElevation);
        float gamma = std::acos(std::clamp(dot(viewVec, solarVec), -1.0f, 1.0f));
        
#ifdef Use_Spectral_Representation
        float sampledValues[NumChannels];
        if (gamma < m_skyModelStates[0]->solar_radius) {
            for (int i = 0; i < NumChannels; ++i)
                sampledValues[i] = arhosekskymodel_solar_radiance(m_skyModelStates[i], theta, gamma, SampledWavelengths[i]);
        }
        else {
            for (int i = 0; i < NumChannels; ++i)
                sampledValues[i] = arhosekskymodel_radiance(m_skyModelStates[i], theta, gamma, SampledWavelengths[i]);
        }
        
        RegularContinuousSpectrum spectrum(SampledWavelengths[0], SampledWavelengths[NumChannels - 1], sampledValues, NumChannels);
#else
        SampledSpectrum spectrum;
        for (int i = 0; i < NumChannels; ++i)
            spectrum[i] = RadianceScale * arhosek_tristim_skymodel_radiance(m_skyModelStates[i], theta, gamma, i);
#endif
        SampledSpectrum ret = spectrum.evaluate(wls);
        SLRAssert(!ret.hasMinus() && ret.allFinite(), "Invalid vaue.");
        
        return ret;
    }
    
    ContinuousDistribution2D* AnalyticSkySpectrumTexture::createIBLImportanceMap() const {        
        const uint32_t mapWidth = 1024;
        const uint32_t mapHeight = 512;
        std::function<float(uint32_t, uint32_t)> pickFunc = [this, &mapWidth, &mapHeight](uint32_t x, uint32_t y) -> float {
            float theta = M_PI * (y + 0.5f) / mapHeight;
            if (theta >= M_PI / 2)
                return 0.0f;
            Vector3D viewVec = Vector3D::fromPolarYUp(2 * M_PI * (x + 0.5f) / mapWidth, theta);
            Vector3D solarVec = Vector3D::fromPolarYUp(M_PI, M_PI / 2 - m_solarElevation);
            float gamma = std::acos(std::clamp(dot(viewVec, solarVec), -1.0f, 1.0f));
            
#ifdef Use_Spectral_Representation
            float sampledValues[NumChannels];
            if (gamma < m_skyModelStates[0]->solar_radius) {
                for (int i = 0; i < NumChannels; ++i)
                    sampledValues[i] = arhosekskymodel_solar_radiance(m_skyModelStates[i], theta, gamma, SampledWavelengths[i]);
            }
            else {
                for (int i = 0; i < NumChannels; ++i)
                    sampledValues[i] = arhosekskymodel_radiance(m_skyModelStates[i], theta, gamma, SampledWavelengths[i]);
            }
            
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
                spectrum[i] = RadianceScale * arhosek_tristim_skymodel_radiance(m_skyModelStates[i], theta, gamma, i);
            
            float luminance = spectrum.luminance();
#endif
            SLRAssert(std::isfinite(luminance), "Invalid area average value.");
            return std::sin(M_PI * (y + 0.5f) / mapHeight) * luminance;
        };
        
        RegularConstantContinuousDistribution2D* ret = new RegularConstantContinuousDistribution2D(mapWidth, mapHeight, pickFunc);
//        ret->exportBMP("distribution.bmp", true);
        
        return ret;
    }
}

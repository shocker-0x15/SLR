//
//  AnalyticSkySpectrumTexture.h
//
//  Created by 渡部 心 on 2017/04/24.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef __SLR_AnalyticSkySpectrumTexture__
#define __SLR_AnalyticSkySpectrumTexture__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/textures.h"
extern "C" {
#include "../External/AnalyticSkyDomeRadiance/ArHosekSkyModel.h"
}

namespace SLR {
    // References
    // An Analytic Model for Full Spectral Sky-Dome Radiance
    class SLR_API AnalyticSkySpectrumTexture : public SpectrumTexture {
        class SunDiscContinuousDistribution2D;
        
#ifdef SLR_Use_Spectral_Representation
        static const uint32_t NumChannels;
        static const float SampledWavelengths[11];
#else
        static const uint32_t NumChannels;
        static const float RadianceScale;
#endif
        
        float m_solarRadius;
        float m_solarElevation;
        float m_turbidity;
        const AssetSpectrum* m_groundAlbedo;
#ifdef SLR_Use_Spectral_Representation
        ArHosekSkyModelState* m_skyModelStates[11];
#else
        ArHosekSkyModelState* m_skyModelStates[3];
#endif
        const Texture2DMapping* m_mapping;
        
        Vector3D m_sunDirection;
        mutable ContinuousDistribution2D* m_distribution;
        mutable SunDiscContinuousDistribution2D* m_sunDiscDistribution;
        mutable RegularConstantContinuousDistribution2D* m_skyDomeDistribution;
    public:
        AnalyticSkySpectrumTexture(float solarRadius, float solarElevation, float turbidity, const AssetSpectrum* groundAlbedo, const Texture2DMapping* mapping);
        ~AnalyticSkySpectrumTexture();
        
        SampledSpectrum evaluate(const Point3D &p, const WavelengthSamples &wls) const;
        SampledSpectrum evaluate(const SurfacePoint &surfPt, const WavelengthSamples &wls) const override {
            return evaluate(m_mapping->map(surfPt), wls);
        }
        SampledSpectrum evaluate(const MediumPoint &medPt, const WavelengthSamples &wls) const override {
            return evaluate(m_mapping->map(medPt), wls);
        }
        const ContinuousDistribution2D* createIBLImportanceMap() const override;
    };
}

#endif /* __SLR_AnalyticSkySpectrumTexture__ */

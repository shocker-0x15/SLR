//
//  MicrofacetSurfaceMaterial.h
//
//  Created by 渡部 心 on 2016/05/04.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef __MicrofacetSurfaceMaterial__
#define __MicrofacetSurfaceMaterial__

#include "../defines.h"
#include "../references.h"
#include "../Core/surface_material.h"

namespace SLR {
    class SLR_API MicrofacetReflection : public SurfaceMaterial {
        const SpectrumTexture* m_eta;
        const SpectrumTexture* m_k;
        const SVMicrofacetDistribution* m_D;
    public:
        MicrofacetReflection(const SpectrumTexture* eta, const SpectrumTexture* k, const SVMicrofacetDistribution* D) :
        m_eta(eta), m_k(k), m_D(D) {}
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
    
    class SLR_API MicrofacetScattering : public SurfaceMaterial {
        const SpectrumTexture* m_etaExt;
        const SpectrumTexture* m_etaInt;
        const SVMicrofacetDistribution* m_D;
    public:
        MicrofacetScattering(const SpectrumTexture* etaExt, const SpectrumTexture* etaInt, const SVMicrofacetDistribution* D) :
        m_etaExt(etaExt), m_etaInt(etaInt), m_D(D) {}
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
}

#endif /* __MicrofacetSurfaceMaterial__ */

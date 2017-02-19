//
//  basic_surface_materials.h
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_basic_surface_materials__
#define __SLR_basic_surface_materials__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/surface_material.h"

namespace SLR {
    class SLR_API DiffuseReflectionSurfaceMaterial : public SurfaceMaterial {
        const SpectrumTexture* m_reflectance;
        const FloatTexture* m_sigma;
    public:
        DiffuseReflectionSurfaceMaterial(const SpectrumTexture* reflectance, const FloatTexture* sigma) :
        m_reflectance(reflectance), m_sigma(sigma) {}
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
    
    
    
    class SLR_API SpecularReflectionSurfaceMaterial : public SurfaceMaterial {
        const SpectrumTexture* m_coeffR;
        const SpectrumTexture* m_eta;
        const SpectrumTexture* m_k;
    public:
        SpecularReflectionSurfaceMaterial(const SpectrumTexture* coeffR, const SpectrumTexture* eta, const SpectrumTexture* k) :
        m_coeffR(coeffR), m_eta(eta), m_k(k) { }
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
    
    
    
    class SLR_API SpecularScatteringSurfaceMaterial : public SurfaceMaterial {
        const SpectrumTexture* m_coeff;
        const SpectrumTexture* m_etaExt;
        const SpectrumTexture* m_etaInt;
    public:
        SpecularScatteringSurfaceMaterial(const SpectrumTexture* coeff, const SpectrumTexture* etaExt, const SpectrumTexture* etaInt) :
        m_coeff(coeff), m_etaExt(etaExt), m_etaInt(etaInt) { }
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
    
    
    
    class SLR_API FlippedSurfaceMaterial : public SurfaceMaterial {
        const SurfaceMaterial* m_baseMat;
    public:
        FlippedSurfaceMaterial(const SurfaceMaterial* baseMat) :
        m_baseMat(baseMat) { }
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
}

#endif /* __SLR_basic_surface_materials__ */

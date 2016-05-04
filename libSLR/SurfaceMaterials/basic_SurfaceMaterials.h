//
//  basic_SurfaceMaterials.h
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__basic_SurfaceMaterials__
#define __SLR__basic_SurfaceMaterials__

#include "../defines.h"
#include "../references.h"
#include "../Core/surface_material.h"

namespace SLR {
    class SLR_API DiffuseReflection : public SurfaceMaterial {
        const SpectrumTexture* m_reflectance;
        const FloatTexture* m_sigma;
    public:
        DiffuseReflection(const SpectrumTexture* reflectance, const FloatTexture* sigma) :
        m_reflectance(reflectance), m_sigma(sigma) {};
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
    
    
    
    class SLR_API SpecularReflection : public SurfaceMaterial {
        const SpectrumTexture* m_coeffR;
        const SVFresnel* m_fresnel;
    public:
        SpecularReflection(const SpectrumTexture* coeffR, const SVFresnel* fresnel) :
        m_coeffR(coeffR), m_fresnel(fresnel) { };
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
    
    
    
    class SLR_API SpecularTransmission : public SurfaceMaterial {
        const SpectrumTexture* m_coeffT;
        const SpectrumTexture* m_etaExt;
        const SpectrumTexture* m_etaInt;
    public:
        SpecularTransmission(const SpectrumTexture* coeffT, const SpectrumTexture* etaExt, const SpectrumTexture* etaInt) :
        m_coeffT(coeffT), m_etaExt(etaExt), m_etaInt(etaInt) { };
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
    
    
    
    class SLR_API InverseSurfaceMaterial : public SurfaceMaterial {
        const SurfaceMaterial* m_baseMat;
    public:
        InverseSurfaceMaterial(const SurfaceMaterial* baseMat) :
        m_baseMat(baseMat) { };
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
}

#endif /* defined(__SLR__basic_SurfaceMaterials__) */

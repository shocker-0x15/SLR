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
    class DiffuseReflection : public SurfaceMaterial {
        SpectrumTextureRef m_reflectance;
        FloatTextureRef m_sigma;
    public:
        DiffuseReflection(const SpectrumTextureRef &reflectance, const FloatTextureRef &sigma) :
        m_reflectance(reflectance), m_sigma(sigma) {};
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
    
    class SpecularReflection : public SurfaceMaterial {
        SpectrumTextureRef m_coeffR;
        SpatialFresnelRef m_fresnel;
    public:
        SpecularReflection(const SpectrumTextureRef &coeffR, const SpatialFresnelRef &fresnel) :
        m_coeffR(coeffR), m_fresnel(fresnel) { };
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
    
    class SpecularTransmission : public SurfaceMaterial {
        SpectrumTextureRef m_coeffT;
        SpectrumTextureRef m_etaExt;
        SpectrumTextureRef m_etaInt;
    public:
        SpecularTransmission(const SpectrumTextureRef &coeffT, const SpectrumTextureRef &etaExt, const SpectrumTextureRef &etaInt) :
        m_coeffT(coeffT), m_etaExt(etaExt), m_etaInt(etaInt) { };
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };    
}

#endif /* defined(__SLR__basic_SurfaceMaterials__) */

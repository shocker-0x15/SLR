//
//  AshikhminReflection.h
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__AshikhminReflection__
#define __SLR__AshikhminReflection__

#include "../defines.h"
#include "../references.h"
#include "../Core/surface_material.h"

namespace SLR {
    class AshikhminSpecularReflection : public SurfaceMaterial {
        const SpectrumTexture* m_Rs;
        const FloatTexture* m_nu;
        const FloatTexture* m_nv;
    public:
        AshikhminSpecularReflection(const SpectrumTexture* Rs, const FloatTexture* nu, const FloatTexture* nv) :
        m_Rs(Rs), m_nu(nu), m_nv(nv) { };
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
    
    class AshikhminDiffuseReflection : public SurfaceMaterial {
        const SpectrumTexture* m_Rs;
        const SpectrumTexture* m_Rd;
    public:
        AshikhminDiffuseReflection(const SpectrumTexture* Rs, const SpectrumTexture* Rd) :
        m_Rs(Rs), m_Rd(Rd) { };
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
}

#endif /* defined(__SLR__AshikhminReflection__) */

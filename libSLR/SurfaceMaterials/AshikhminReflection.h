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
        SpectrumTextureRef m_Rs;
        FloatTextureRef m_nu, m_nv;
    public:
        AshikhminSpecularReflection(const SpectrumTextureRef &Rs, const FloatTextureRef &nu, const FloatTextureRef &nv) :
        m_Rs(Rs), m_nu(nu), m_nv(nv) { };
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
    
    class AshikhminDiffuseReflection : public SurfaceMaterial {
        SpectrumTextureRef m_Rs;
        SpectrumTextureRef m_Rd;
    public:
        AshikhminDiffuseReflection(const SpectrumTextureRef &Rs, const SpectrumTextureRef &Rd) :
        m_Rs(Rs), m_Rd(Rd) { };
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };    
}

#endif /* defined(__SLR__AshikhminReflection__) */

//
//  MixedSurfaceMaterial.h
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__MixedSurfaceMaterial__
#define __SLR__MixedSurfaceMaterial__

#include "../defines.h"
#include "../references.h"
#include "../Core/surface_material.h"

namespace SLR {
    class MixedSurfaceMaterial : public SurfaceMaterial {
        SurfaceMaterialRef m_mat0;
        SurfaceMaterialRef m_mat1;
        FloatTextureRef m_factor;
    public:
        MixedSurfaceMaterial(const SurfaceMaterialRef &m0, const SurfaceMaterialRef &m1, const FloatTextureRef &factor) :
        m_mat0(m0), m_mat1(m1), m_factor(factor) {};
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };    
}

#endif /* defined(__SLR__MixedSurfaceMaterial__) */

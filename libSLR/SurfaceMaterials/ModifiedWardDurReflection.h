//
//  ModifiedWardDurReflection.h
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_ModifiedWardDurReflection__
#define __SLR_ModifiedWardDurReflection__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/surface_material.h"

namespace SLR {
    class SLR_API ModifiedWardDurReflectionSurfaceMaterial : public SurfaceMaterial {
        const SpectrumTexture* m_reflectance;
        const FloatTexture* m_anisoX;
        const FloatTexture* m_anisoY;
    public:
        ModifiedWardDurReflectionSurfaceMaterial(const SpectrumTexture* reflectance, const FloatTexture* anisoX, const FloatTexture* anisoY) :
        m_reflectance(reflectance), m_anisoX(anisoX), m_anisoY(anisoY) { };
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
}

#endif /* __SLR_ModifiedWardDurReflection__ */

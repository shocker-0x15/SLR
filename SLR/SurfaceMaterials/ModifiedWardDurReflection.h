//
//  ModifiedWardDurReflection.h
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__ModifiedWardDurReflection__
#define __SLR__ModifiedWardDurReflection__

#include "../defines.h"
#include "../references.h"
#include "../Core/surface_material.h"

class ModifiedWardDurReflection : public SurfaceMaterial {
    SpectrumTextureRef m_reflectance;
    FloatTextureRef m_anisoX, m_anisoY;
public:
    ModifiedWardDurReflection(const SpectrumTextureRef &reflectance, const FloatTextureRef &anisoX, const FloatTextureRef &anisoY) :
    m_reflectance(reflectance), m_anisoX(anisoX), m_anisoY(anisoY) { };
    
    BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
};

#endif /* defined(__SLR__ModifiedWardDurReflection__) */

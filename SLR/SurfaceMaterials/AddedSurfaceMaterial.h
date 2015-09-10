//
//  AddedSurfaceMaterial.h
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__AddedSurfaceMaterial__
#define __SLR__AddedSurfaceMaterial__

#include "../defines.h"
#include "../references.h"
#include "../Core/surface_material.h"

class AddedSurfaceMaterial : public SurfaceMaterial {
    SurfaceMaterialRef m_mat0;
    SurfaceMaterialRef m_mat1;
public:
    AddedSurfaceMaterial(const SurfaceMaterialRef &m0, const SurfaceMaterialRef &m1) :
    m_mat0(m0), m_mat1(m1) {};
    
    BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
};

#endif /* defined(__SLR__AddedSurfaceMaterial__) */

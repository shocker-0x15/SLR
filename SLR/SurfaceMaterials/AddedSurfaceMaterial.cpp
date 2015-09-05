//
//  AddedSurfaceMaterial.cpp
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "AddedSurfaceMaterial.h"
#include "../BSDFs/MultiBSDF.h"
#include "../Memory/ArenaAllocator.h"

BSDF* AddedSurfaceMaterial::getBSDF(const SurfacePoint &surfPt, ArenaAllocator &mem, float scale) const {
    BSDF* bsdf0 = m_mat0->getBSDF(surfPt, mem, scale);
    BSDF* bsdf1 = m_mat1->getBSDF(surfPt, mem, scale);
    MultiBSDF* bsdf = mem.create<MultiBSDF>();
    bsdf->add(bsdf0);
    bsdf->add(bsdf1);
    return bsdf;
}

//
//  SummedSurfaceMaterial.cpp
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "SummedSurfaceMaterial.h"

#include "../MemoryAllocators/ArenaAllocator.h"
#include "../BSDF/MultiBSDF.h"

namespace SLR {
    BSDF* SummedSurfaceMaterial::getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        BSDF* bsdf0 = m_mat0->getBSDF(surfPt, wls, mem, scale);
        BSDF* bsdf1 = m_mat1->getBSDF(surfPt, wls, mem, scale);
        MultiBSDF* bsdf = mem.create<MultiBSDF>();
        bsdf->add(bsdf0);
        bsdf->add(bsdf1);
        return bsdf;
    }
}

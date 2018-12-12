//
//  SummedSurfaceMaterial.h
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_SummedSurfaceMaterial__
#define __SLR_SummedSurfaceMaterial__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/surface_material.h"

namespace SLR {
    class SLR_API SummedSurfaceMaterial : public SurfaceMaterial {
        const SurfaceMaterial* m_mat0;
        const SurfaceMaterial* m_mat1;
    public:
        SummedSurfaceMaterial(const SurfaceMaterial* m0, const SurfaceMaterial* m1) :
        m_mat0(m0), m_mat1(m1) {};
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
}

#endif /* __SLR_SummedSurfaceMaterial__ */

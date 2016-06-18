//
//  AshikhminShirleyReflection.cpp
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "AshikhminReflection.h"
#include "../BSDFs/AshikhminBRDF.h"
#include "../Memory/ArenaAllocator.h"
#include "../Core/textures.h"

namespace SLR {
    BSDF* AshikhminShirleyReflection::getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        SampledSpectrum Rs = m_Rs->evaluate(surfPt.texCoord, wls);
        SampledSpectrum Rd = m_Rd->evaluate(surfPt.texCoord, wls);
        float nu = m_nu->evaluate(surfPt.texCoord);
        float nv = m_nv->evaluate(surfPt.texCoord);
        return mem.create<AshikhminShirleyBRDF>(scale * Rs, scale * Rd, nu, nv);
    }
}

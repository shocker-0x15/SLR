//
//  DiffuseEmission.cpp
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "DiffuseEmission.h"
#include "../EDFs/basic_EDFs.h"
#include "../Memory/ArenaAllocator.h"
#include "../Core/textures.h"

namespace SLR {
    SampledSpectrum DiffuseEmitterSurfaceProperty::emittance(const SurfacePoint &surfPt, const WavelengthSamples &wls) const {
        return m_emittance->evaluate(surfPt, wls);
    }
    
    EDF* DiffuseEmitterSurfaceProperty::getEDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        return mem.create<DiffuseEDF>();
    }    
}

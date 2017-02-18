//
//  IdealDirectionalEmission.cpp
//
//  Created by 渡部 心 on 2017/02/17.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "IdealDirectionalEmission.h"
#include "../EDFs/basic_EDFs.h"
#include "../MemoryAllocators/ArenaAllocator.h"
#include "../Core/textures.h"

namespace SLR {
    SampledSpectrum IdealDirectionalEmitterSurfaceProperty::emittance(const SurfacePoint &surfPt, const WavelengthSamples &wls) const {
        return m_emittance->evaluate(surfPt, wls);
    }
    
    EDF* IdealDirectionalEmitterSurfaceProperty::getEDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        return mem.create<IdealDirectionalEDF>();
    }
}

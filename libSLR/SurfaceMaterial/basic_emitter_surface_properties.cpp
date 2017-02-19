//
//  basic_emitter_surface_properties.cpp
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "basic_emitter_surface_properties.h"

#include "../MemoryAllocators/ArenaAllocator.h"
#include "../Core/textures.h"
#include "../EDF/basic_edfs.h"

namespace SLR {
    SampledSpectrum DiffuseEmitterSurfaceProperty::emittance(const SurfacePoint &surfPt, const WavelengthSamples &wls) const {
        return m_emittance->evaluate(surfPt, wls);
    }
    
    EDF* DiffuseEmitterSurfaceProperty::getEDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        return mem.create<DiffuseEDF>();
    }
    
    
    
    SampledSpectrum IdealDirectionalEmitterSurfaceProperty::emittance(const SurfacePoint &surfPt, const WavelengthSamples &wls) const {
        return m_emittance->evaluate(surfPt, wls);
    }
    
    EDF* IdealDirectionalEmitterSurfaceProperty::getEDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        return mem.create<IdealDirectionalEDF>();
    }
}

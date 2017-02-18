//
//  basic_MediumMaterials.cpp
//
//  Created by 渡部 心 on 2017/02/04.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "basic_MediumMaterials.h"
#include "../MemoryAllocators/ArenaAllocator.h"
#include "../Core/textures.h"
#include "../PhaseFunctions/basic_phase_functions.h"

namespace SLR {
    PhaseFunction* IsotropicScatteringMediumMaterial::getPhaseFunction(const MediumPoint &medPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        return mem.create<IsotropicPhaseFunction>();
    }
    
    
    
    PhaseFunction* HenyeyGreensteinScatteringMediumMaterial::getPhaseFunction(const MediumPoint &medPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        return mem.create<HenyeyGreensteinPhaseFunction>(m_g->evaluate(medPt));
    }
    
    
    
    PhaseFunction* SchlickScatteringMediumMaterial::getPhaseFunction(const MediumPoint &medPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        return mem.create<SchlickPhaseFunction>(m_g->evaluate(medPt));
    }
}

//
//  basic_MediumMaterials.h
//
//  Created by 渡部 心 on 2017/02/04.
//  Copyright © 2017年 渡部 心. All rights reserved.
//

#ifndef __SLR_basic_MediumMaterials__
#define __SLR_basic_MediumMaterials__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/medium_material.h"

namespace SLR {
    class SLR_API IsotropicScatteringMediumMaterial : public MediumMaterial {
    public:
        IsotropicScatteringMediumMaterial() {}
        
        PhaseFunction* getPhaseFunction(const MediumPoint &medPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const;
    };
    
    
    
    class SLR_API HenyeyGreensteinScatteringMediumMaterial : public MediumMaterial {
        const FloatTexture* m_g;
    public:
        HenyeyGreensteinScatteringMediumMaterial(const FloatTexture* g) :
        m_g(g) {}
        
        PhaseFunction* getPhaseFunction(const MediumPoint &medPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const;
    };
    
    
    
    class SLR_API SchlickScatteringMediumMaterial : public MediumMaterial {
        const FloatTexture* m_g;
    public:
        SchlickScatteringMediumMaterial(const FloatTexture* g) :
        m_g(g) {}
        
        PhaseFunction* getPhaseFunction(const MediumPoint &medPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const;
    };
}

#endif /* __SLR_basic_MediumMaterials__ */

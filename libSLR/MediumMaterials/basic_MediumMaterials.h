//
//  basic_MediumMaterials.h
//
//  Created by 渡部 心 on 2017/02/04.
//  Copyright © 2017年 渡部 心. All rights reserved.
//

#ifndef __SLR_basic_MediumMaterials__
#define __SLR_basic_MediumMaterials__

#include "../defines.h"
#include "../references.h"
#include "../Core/medium_material.h"

namespace SLR {
    class SLR_API IsotropicScattering : public MediumMaterial {
    public:
        IsotropicScattering() {}
        
        PhaseFunction* getPhaseFunction(const MediumPoint &medPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const;
    };
    
    
    
    class SLR_API HenyeyGreensteinScattering : public MediumMaterial {
        const FloatTexture* m_g;
    public:
        HenyeyGreensteinScattering(const FloatTexture* g) :
        m_g(g) {}
        
        PhaseFunction* getPhaseFunction(const MediumPoint &medPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const;
    };
}

#endif /* basic_MediumMaterials_hpp */

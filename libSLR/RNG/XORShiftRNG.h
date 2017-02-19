//
//  XORShiftRNG.h
//
//  Created by 渡部 心 on 11/07/21.
//  Copyright (c) 2014年 渡部 心. All rights reserved.
//

#ifndef __SLR_XORShiftRNG__
#define __SLR_XORShiftRNG__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/random_number_generator.h"

namespace SLR {
    template <typename TypeSet>
    class SLR_API XORShiftRNGTemplate : public RandomNumberGeneratorTemplate<TypeSet> {
        typename TypeSet::UInt m_state[4];
    public:
        XORShiftRNGTemplate();
        XORShiftRNGTemplate(typename TypeSet::Int seed);
        XORShiftRNGTemplate(const XORShiftRNGTemplate &rng) {
            m_state[0] = rng.m_state[0];
            m_state[1] = rng.m_state[1];
            m_state[2] = rng.m_state[2];
            m_state[3] = rng.m_state[3];
        };
        
        typename TypeSet::UInt getUInt() override;
    };
    
    template <> XORShiftRNGTemplate<Types32bit>::XORShiftRNGTemplate();
    template <> XORShiftRNGTemplate<Types32bit>::XORShiftRNGTemplate(Types32bit::Int seed);
    template <> Types32bit::UInt XORShiftRNGTemplate<Types32bit>::getUInt();    
}

#endif /* __SLR_XORShiftRNG__ */

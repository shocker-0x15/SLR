//
//  XORShift.h
//
//  Created by 渡部 心 on 11/07/21.
//  Copyright (c) 2014年 渡部 心. All rights reserved.
//

#ifndef __SLR__XORShift__
#define __SLR__XORShift__

#include "../defines.h"
#include "../references.h"
#include "RandomNumberGenerator.h"

namespace SLR {
    template <typename TypeSet>
    class XORShiftTemplate : public RandomNumberGeneratorTemplate<TypeSet> {
        typename TypeSet::UInt m_state[4];
    public:
        XORShiftTemplate();
        XORShiftTemplate(typename TypeSet::Int seed);
        XORShiftTemplate(const XORShiftTemplate &xorShift) {
            m_state[0] = xorShift.m_state[0];
            m_state[1] = xorShift.m_state[1];
            m_state[2] = xorShift.m_state[2];
            m_state[3] = xorShift.m_state[3];
        };
        
        typename TypeSet::UInt getUInt() override;
    };
    
    template <> XORShiftTemplate<Types32bit>::XORShiftTemplate();
    template <> XORShiftTemplate<Types32bit>::XORShiftTemplate(Types32bit::Int seed);
    template <> Types32bit::UInt XORShiftTemplate<Types32bit>::getUInt();    
}

#endif

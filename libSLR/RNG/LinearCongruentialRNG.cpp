//
//  LinearCongruentialRNG.cpp
//
//  Created by 渡部 心 on 2016/06/17.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "LinearCongruentialRNG.h"

namespace SLR {
    template <>
    SLR_API Types32bit::UInt LinearCongruentialRNGTemplate<Types32bit>::getUInt() {
        return (m_seed = (m_seed * 1103515245 + 12345) & 0xFFFFFFFF);
    }
    
    template <>
    SLR_API Types64bit::UInt LinearCongruentialRNGTemplate<Types64bit>::getUInt() {
        return (m_seed = (m_seed * 2862933555777941757 + 3037000493) & 0xFFFFFFFFFFFFFFFF);
    }
    
    template class SLR_API LinearCongruentialRNGTemplate<Types32bit>;
    template class SLR_API LinearCongruentialRNGTemplate<Types64bit>;
}

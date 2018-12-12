//
//  random_number_generator.cpp
//
//  Created by 渡部 心 on 2015/07/07.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "random_number_generator.h"

namespace SLR {
    template <>
    SLR_API Types32bit::Float RandomNumberGeneratorTemplate<Types32bit>::getFloat0cTo1o() {
        Types32bit::UInt fractionBits = (getUInt() >> 9) | 0x3f800000;
        return *(Types32bit::Float*)&fractionBits - 1.0f;
    }
    
    template <>
    SLR_API Types64bit::Float RandomNumberGeneratorTemplate<Types64bit>::getFloat0cTo1o() {
        Types64bit::UInt fractionBits = (getUInt() >> 12) | 0x3ff0000000000000;
        return *(Types64bit::Float*)&fractionBits - 1.0;
    }
    
    template class SLR_API RandomNumberGeneratorTemplate<Types32bit>;
    template class SLR_API RandomNumberGeneratorTemplate<Types64bit>;
}

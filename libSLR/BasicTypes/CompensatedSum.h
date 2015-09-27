//
//  CompensatedSum.hpp
//
//  Created by 渡部 心 on 2015/07/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef SLR_CompensatedSum_h
#define SLR_CompensatedSum_h

#include "../defines.h"
#include "../references.h"

namespace SLR {
    template <typename RealType>
    struct CompensatedSum {
        RealType result;
        RealType comp;
        CompensatedSum(const RealType &value) : result(value), comp(0.0) { };
        CompensatedSum &operator=(const RealType &value) {
            result = value;
            comp = 0;
        };
        CompensatedSum &operator+=(const RealType &value) {
            RealType cInput = value - comp;
            RealType sumTemp = result + cInput;
            comp = (sumTemp - result) - cInput;
            result = sumTemp;
            return *this;
        };
        operator RealType() const { return result; };
    };    
}

#endif

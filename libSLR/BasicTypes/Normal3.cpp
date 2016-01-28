//
//  Normal3.cpp
//
//  Created by 渡部 心 on 2016/01/26.
//  Copyright c 2016年 渡部 心. All rights reserved.
//

#include "Normal3.h"

namespace SLR {
    template struct SLR_API Normal3Template<float>;
//    template SLR_API const Normal3Template<float> Normal3Template<float>::Zero;
//    template SLR_API const Normal3Template<float> Normal3Template<float>::One;
//    template SLR_API const Normal3Template<float> Normal3Template<float>::Inf;
//    template SLR_API const Normal3Template<float> Normal3Template<float>::NaN;

    template struct SLR_API Normal3Template<double>;
//    template SLR_API const Normal3Template<double> Normal3Template<double>::Zero;
//    template SLR_API const Normal3Template<double> Normal3Template<double>::One;
//    template SLR_API const Normal3Template<double> Normal3Template<double>::Inf;
//    template SLR_API const Normal3Template<double> Normal3Template<double>::NaN;
}

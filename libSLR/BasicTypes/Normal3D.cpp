//
//  Normal3.cpp
//
//  Created by 渡部 心 on 2016/01/26.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "Normal3D.h"

namespace SLR {
    template struct SLR_API Normal3DTemplate<float>;
//    template SLR_API const Normal3DTemplate<float> Normal3DTemplate<float>::Zero;
//    template SLR_API const Normal3DTemplate<float> Normal3DTemplate<float>::One;
//    template SLR_API const Normal3DTemplate<float> Normal3DTemplate<float>::Inf;
//    template SLR_API const Normal3DTemplate<float> Normal3DTemplate<float>::NaN;

    template struct SLR_API Normal3DTemplate<double>;
//    template SLR_API const Normal3DTemplate<double> Normal3DTemplate<double>::Zero;
//    template SLR_API const Normal3DTemplate<double> Normal3DTemplate<double>::One;
//    template SLR_API const Normal3DTemplate<double> Normal3DTemplate<double>::Inf;
//    template SLR_API const Normal3DTemplate<double> Normal3DTemplate<double>::NaN;
}

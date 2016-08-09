//
//  Vector3.cpp
//
//  Created by 渡部 心 on 2016/01/26.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "Vector3.h"

namespace SLR {
    template struct SLR_API Vector3Template<float>;
//    template SLR_API const Vector3Template<float> Vector3Template<float>::Zero;
//    template SLR_API const Vector3Template<float> Vector3Template<float>::One;
//    template SLR_API const Vector3Template<float> Vector3Template<float>::Inf;
//    template SLR_API const Vector3Template<float> Vector3Template<float>::NaN;
//    template SLR_API const Vector3Template<float> Vector3Template<float>::Ex;
//    template SLR_API const Vector3Template<float> Vector3Template<float>::Ey;
//    template SLR_API const Vector3Template<float> Vector3Template<float>::Ez;

    template struct SLR_API Vector3Template<double>;
//    template SLR_API const Vector3Template<double> Vector3Template<double>::Zero;
//    template SLR_API const Vector3Template<double> Vector3Template<double>::One;
//    template SLR_API const Vector3Template<double> Vector3Template<double>::Inf;
//    template SLR_API const Vector3Template<double> Vector3Template<double>::NaN;
//    template SLR_API const Vector3Template<double> Vector3Template<double>::Ex;
//    template SLR_API const Vector3Template<double> Vector3Template<double>::Ey;
//    template SLR_API const Vector3Template<double> Vector3Template<double>::Ez;
}

//
//  Vector4.cpp
//
//  Created by 渡部 心 on 2016/01/26.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "Vector4.h"

namespace SLR {
    template struct SLR_API Vector4Template<float>;
//    template SLR_API const Vector4Template<float> Vector4Template<float>::Zero;
//    template SLR_API const Vector4Template<float> Vector4Template<float>::One;
//    template SLR_API const Vector4Template<float> Vector4Template<float>::Inf;
//    template SLR_API const Vector4Template<float> Vector4Template<float>::NaN;
//    template SLR_API const Vector4Template<float> Vector4Template<float>::Ex;
//    template SLR_API const Vector4Template<float> Vector4Template<float>::Ey;
//    template SLR_API const Vector4Template<float> Vector4Template<float>::Ez;
//    template SLR_API const Vector4Template<float> Vector4Template<float>::Ew;

    template struct SLR_API Vector4Template<double>;
//    template SLR_API const Vector4Template<double> Vector4Template<float>::Zero;
//    template SLR_API const Vector4Template<double> Vector4Template<float>::One;
//    template SLR_API const Vector4Template<double> Vector4Template<float>::Inf;
//    template SLR_API const Vector4Template<double> Vector4Template<float>::NaN;
//    template SLR_API const Vector4Template<double> Vector4Template<float>::Ex;
//    template SLR_API const Vector4Template<double> Vector4Template<float>::Ey;
//    template SLR_API const Vector4Template<double> Vector4Template<float>::Ez;
//    template SLR_API const Vector4Template<double> Vector4Template<float>::Ew;
}

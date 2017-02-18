//
//  Vector4.cpp
//
//  Created by 渡部 心 on 2016/01/26.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "Vector4.h"

namespace SLR {
    template struct SLR_API Vector4DTemplate<float>;
//    template SLR_API const Vector4DTemplate<float> Vector4DTemplate<float>::Zero;
//    template SLR_API const Vector4DTemplate<float> Vector4DTemplate<float>::One;
//    template SLR_API const Vector4DTemplate<float> Vector4DTemplate<float>::Inf;
//    template SLR_API const Vector4DTemplate<float> Vector4DTemplate<float>::NaN;
//    template SLR_API const Vector4DTemplate<float> Vector4DTemplate<float>::Ex;
//    template SLR_API const Vector4DTemplate<float> Vector4DTemplate<float>::Ey;
//    template SLR_API const Vector4DTemplate<float> Vector4DTemplate<float>::Ez;
//    template SLR_API const Vector4DTemplate<float> Vector4DTemplate<float>::Ew;

    template struct SLR_API Vector4DTemplate<double>;
//    template SLR_API const Vector4DTemplate<double> Vector4DTemplate<float>::Zero;
//    template SLR_API const Vector4DTemplate<double> Vector4DTemplate<float>::One;
//    template SLR_API const Vector4DTemplate<double> Vector4DTemplate<float>::Inf;
//    template SLR_API const Vector4DTemplate<double> Vector4DTemplate<float>::NaN;
//    template SLR_API const Vector4DTemplate<double> Vector4DTemplate<float>::Ex;
//    template SLR_API const Vector4DTemplate<double> Vector4DTemplate<float>::Ey;
//    template SLR_API const Vector4DTemplate<double> Vector4DTemplate<float>::Ez;
//    template SLR_API const Vector4DTemplate<double> Vector4DTemplate<float>::Ew;
}

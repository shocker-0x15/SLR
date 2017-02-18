//
//  Vector3.cpp
//
//  Created by 渡部 心 on 2016/01/26.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "Vector3D.h"

namespace SLR {
    template struct SLR_API Vector3DTemplate<float>;
//    template SLR_API const Vector3DTemplate<float> Vector3DTemplate<float>::Zero;
//    template SLR_API const Vector3DTemplate<float> Vector3DTemplate<float>::One;
//    template SLR_API const Vector3DTemplate<float> Vector3DTemplate<float>::Inf;
//    template SLR_API const Vector3DTemplate<float> Vector3DTemplate<float>::NaN;
//    template SLR_API const Vector3DTemplate<float> Vector3DTemplate<float>::Ex;
//    template SLR_API const Vector3DTemplate<float> Vector3DTemplate<float>::Ey;
//    template SLR_API const Vector3DTemplate<float> Vector3DTemplate<float>::Ez;

    template struct SLR_API Vector3DTemplate<double>;
//    template SLR_API const Vector3DTemplate<double> Vector3DTemplate<double>::Zero;
//    template SLR_API const Vector3DTemplate<double> Vector3DTemplate<double>::One;
//    template SLR_API const Vector3DTemplate<double> Vector3DTemplate<double>::Inf;
//    template SLR_API const Vector3DTemplate<double> Vector3DTemplate<double>::NaN;
//    template SLR_API const Vector3DTemplate<double> Vector3DTemplate<double>::Ex;
//    template SLR_API const Vector3DTemplate<double> Vector3DTemplate<double>::Ey;
//    template SLR_API const Vector3DTemplate<double> Vector3DTemplate<double>::Ez;
}

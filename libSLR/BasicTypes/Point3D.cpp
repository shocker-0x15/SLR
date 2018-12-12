//
//  Point3.cpp
//
//  Created by 渡部 心 on 2016/01/26.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "Point3D.h"

namespace SLR {
    template struct SLR_API Point3DTemplate<float>;
//    template SLR_API const Point3DTemplate<float> Point3DTemplate<float>::Zero;
//    template SLR_API const Point3DTemplate<float> Point3DTemplate<float>::One;
//    template SLR_API const Point3DTemplate<float> Point3DTemplate<float>::Inf;
//    template SLR_API const Point3DTemplate<float> Point3DTemplate<float>::NaN;

    template struct SLR_API Point3DTemplate<double>;
//    template SLR_API const Point3DTemplate<double> Point3DTemplate<double>::Zero;
//    template SLR_API const Point3DTemplate<double> Point3DTemplate<double>::One;
//    template SLR_API const Point3DTemplate<double> Point3DTemplate<double>::Inf;
//    template SLR_API const Point3DTemplate<double> Point3DTemplate<double>::NaN;
}

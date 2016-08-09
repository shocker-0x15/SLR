//
//  Point3.cpp
//
//  Created by 渡部 心 on 2016/01/26.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "Point3.h"

namespace SLR {
    template struct SLR_API Point3Template<float>;
//    template SLR_API const Point3Template<float> Point3Template<float>::Zero;
//    template SLR_API const Point3Template<float> Point3Template<float>::One;
//    template SLR_API const Point3Template<float> Point3Template<float>::Inf;
//    template SLR_API const Point3Template<float> Point3Template<float>::NaN;

    template struct SLR_API Point3Template<double>;
//    template SLR_API const Point3Template<double> Point3Template<double>::Zero;
//    template SLR_API const Point3Template<double> Point3Template<double>::One;
//    template SLR_API const Point3Template<double> Point3Template<double>::Inf;
//    template SLR_API const Point3Template<double> Point3Template<double>::NaN;
}

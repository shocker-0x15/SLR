//
//  Ray.h
//
//  Created by 渡部 心 on 2017/01/23.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef __SLR_Ray__
#define __SLR_Ray__

#include "Point3.h"
#include "Vector3.h"

namespace SLR {
    template <typename RealType>
    struct SLR_API RayTemplate {
        static const RealType Epsilon;
        
        Point3DTemplate<RealType> org;
        Vector3DTemplate<RealType> dir;
        RealType distMin, distMax;
        RealType time;
        
        RayTemplate() { }
        RayTemplate(const Point3DTemplate<RealType> &o, const Vector3DTemplate<RealType> &d, RealType t, RealType dMin = 0.0f, RealType dMax = INFINITY) :
        org(o), dir(d), distMin(dMin), distMax(dMax), time(t) { }
    };
    
    template <typename RealType>
    const RealType RayTemplate<RealType>::Epsilon = 0.0001f;
}

#endif /* __SLR_Ray__ */

//
//  Ray.h
//
//  Created by 渡部 心 on 2017/01/23.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef __SLR_Ray__
#define __SLR_Ray__

#include "Point3D.h"
#include "Vector3D.h"

namespace SLR {
    template <typename RealType>
    struct SLR_API RayTemplate {
        static const RealType Epsilon;
        
        Point3DTemplate<RealType> org;
        Vector3DTemplate<RealType> dir;
        RealType time;
        
        RayTemplate() { }
        RayTemplate(const Point3DTemplate<RealType> &o, const Vector3DTemplate<RealType> &d, RealType t) :
        org(o), dir(d), time(t) { }
    };
    
    template <typename RealType>
    const RealType RayTemplate<RealType>::Epsilon = 0.0001f;
    
    
    
    template <typename RealType>
    struct SLR_API RaySegmentTemplate {
        RealType distMin, distMax;
        RaySegmentTemplate(RealType dMin = 0.0f, RealType dMax = INFINITY) : 
        distMin(dMin), distMax(dMax) { }
    };
}

#endif /* __SLR_Ray__ */

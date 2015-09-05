//
//  InfiniteSphere.h
//
//  Created by 渡部 心 on 2015/08/11.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__InfiniteSphere__
#define __SLR__InfiniteSphere__

#include "../defines.h"
#include "../references.h"
#include "../Core/geometry.h"

class InfiniteSphere : public Surface {
public:
    InfiniteSphere() { };
    ~InfiniteSphere() { };
    
    BoundingBox3D bounds() const override;
    bool preTransformed() const override;
    bool intersect(const Ray &ray, Intersection* isect) const override;
    void getSurfacePoint(const Intersection &isect, SurfacePoint* surfPt) const override;
    float area() const override;
    void sample(float u0, float u1, SurfacePoint* surfPt, float* areaPDF) const override;
    float evaluateAreaPDF(const SurfacePoint& surfPt) const override;
};

#endif

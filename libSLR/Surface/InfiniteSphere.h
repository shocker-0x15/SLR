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

namespace SLR {
    class SLR_API InfiniteSphere : public Surface {
    public:
        InfiniteSphere() { };
        ~InfiniteSphere() { };
        
        float costForIntersect() const override { return 1.0f; }
        BoundingBox3D bounds() const override;
        BoundingBox3D choppedBounds(BoundingBox3D::Axis chopAxis, float minChopPos, float maxChopPos) const override;
        void splitBounds(BoundingBox3D::Axis chopAxis, float splitPos, BoundingBox3D* bbox0, BoundingBox3D* bbox1) const override;
        bool preTransformed() const override;
        bool intersect(const Ray &ray, SurfaceInteraction* si) const override;
        void getSurfacePoint(const SurfaceInteraction &si, SurfacePoint* surfPt) const override;
        float area() const override;
        void sample(float u0, float u1, SurfacePoint* surfPt, float* areaPDF) const override;
        float evaluateAreaPDF(const SurfacePoint& surfPt) const override;
    };    
}

#endif

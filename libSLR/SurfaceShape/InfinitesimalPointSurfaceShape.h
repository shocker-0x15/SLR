//
//  InfinitesimalPointSurfaceShape.h
//
//  Created by 渡部 心 on 2017/02/17.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef __SLR_InfinitesimalPointSurfaceShape__
#define __SLR_InfinitesimalPointSurfaceShape__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/geometry.h"

namespace SLR {
    class SLR_API InfinitesimalPointSurfaceShape : public SurfaceShape {
        Point3D m_position;
        Vector3D m_direction;
    public:
        InfinitesimalPointSurfaceShape(const Point3D &p, const Vector3D &d) : m_position(p), m_direction(d) { }
        ~InfinitesimalPointSurfaceShape() { }
        
        float costForIntersect() const override { return 0.0f; }
        BoundingBox3D bounds() const override { return m_position; }
        bool preTransformed() const override { return true; }
        bool intersect(const Ray &ray, const RaySegment &segment, SurfaceInteraction* si) const override { return false; }
        void calculateSurfacePoint(const SurfaceInteraction &si, SurfacePoint* surfPt) const override;
        float area() const override { return 0.0f; }
        void sample(float u0, float u1, SurfacePoint* surfPt, float* areaPDF, DirectionType* posType) const override;
        float evaluateAreaPDF(const SurfacePoint& surfPt) const override { return 1.0f; /* delta distribution: \delta(\vx) */ }
    };
}

#endif /* __SLR_InfinitesimalPointSurfaceShape__ */

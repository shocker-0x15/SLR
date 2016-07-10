
//  InfiniteSphere.cpp
//
//  Created by 渡部 心 on 2015/08/11.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "InfiniteSphere.h"
#include "../Core/SurfaceObject.h"
#include "../Core/distributions.h"
#include "../Core/textures.h"

namespace SLR {
    BoundingBox3D InfiniteSphere::bounds() const {
        SLRAssert(false, "InfiniteSphere::bounds() should not be called.");
        return BoundingBox3D(-Point3D::Inf, Point3D::Inf);
    }
    
    BoundingBox3D InfiniteSphere::choppedBounds(BoundingBox3D::Axis chopAxis, float minChopPos, float maxChopPos) const {
        SLRAssert(false, "InfiniteSphere::choppedBounds() should not be called.");
        return BoundingBox3D(-Point3D::Inf, Point3D::Inf);
    }
    
    void InfiniteSphere::splitBounds(BoundingBox3D::Axis chopAxis, float splitPos, BoundingBox3D* bbox0, BoundingBox3D* bbox1) const {
        SLRAssert(false, "InfiniteSphere::splitBounds() should not be called.");
        *bbox0 = BoundingBox3D(-Point3D::Inf, Point3D::Inf);
        *bbox1 = BoundingBox3D(-Point3D::Inf, Point3D::Inf);
    }
    
    bool InfiniteSphere::preTransformed() const {
        return true;
    }
    
    bool InfiniteSphere::intersect(const Ray &ray, Intersection *isect) const {
        if (!std::isinf(ray.distMax))
            return false;
        float theta, phi;
        ray.dir.toPolarYUp(&theta, &phi);
        isect->dist = INFINITY;
        isect->p = Point3D(ray.dir);
        isect->gNormal = -ray.dir;
        isect->u = phi;
        isect->v = theta;
        isect->texCoord = TexCoord2D(phi / (2 * M_PI), theta / M_PI);
        return true;
    }
    
    void InfiniteSphere::getSurfacePoint(const Intersection &isect, SurfacePoint *surfPt) const {
        surfPt->p = isect.p;
        surfPt->atInfinity = true;
        surfPt->gNormal = isect.gNormal;
        surfPt->u = isect.u;
        surfPt->v = isect.v;
        surfPt->texCoord = isect.texCoord;
        surfPt->texCoord0Dir = Vector3D(-std::cos(surfPt->u), 0.0f, -std::sin(surfPt->u));
        surfPt->shadingFrame.x = surfPt->texCoord0Dir;
        surfPt->shadingFrame.z = surfPt->gNormal;
        surfPt->shadingFrame.y = cross(surfPt->shadingFrame.z, surfPt->shadingFrame.x);
    }
    
    float InfiniteSphere::area() const {
        return 4 * M_PI;// * Inf * Inf;
    }
    
    void InfiniteSphere::sample(float u0, float u1, SurfacePoint *surfPt, float *areaPDF) const {
        SLRAssert_NotImplemented();
    }
    
    float InfiniteSphere::evaluateAreaPDF(const SurfacePoint &surfPt) const {
        SLRAssert_NotImplemented();
        return 0.0f;
    }    
}

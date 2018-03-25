
//  InfiniteSphereSurfaceShape.cpp
//
//  Created by 渡部 心 on 2015/08/11.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "InfiniteSphereSurfaceShape.h"

#include "../Core/distributions.h"
#include "../Core/surface_object.h"
#include "../Core/textures.h"

namespace SLR {
    BoundingBox3D InfiniteSphereSurfaceShape::bounds() const {
        SLRAssert(false, "InfiniteSphereSurfaceShape::bounds() should not be called.");
        return BoundingBox3D(-Point3D::Inf, Point3D::Inf);
    }
    
    BoundingBox3D InfiniteSphereSurfaceShape::choppedBounds(BoundingBox3D::Axis chopAxis, float minChopPos, float maxChopPos) const {
        SLRAssert(false, "InfiniteSphereSurfaceShape::choppedBounds() should not be called.");
        return BoundingBox3D(-Point3D::Inf, Point3D::Inf);
    }
    
    void InfiniteSphereSurfaceShape::splitBounds(BoundingBox3D::Axis chopAxis, float splitPos, BoundingBox3D* bbox0, BoundingBox3D* bbox1) const {
        SLRAssert(false, "InfiniteSphereSurfaceShape::splitBounds() should not be called.");
        *bbox0 = BoundingBox3D(-Point3D::Inf, Point3D::Inf);
        *bbox1 = BoundingBox3D(-Point3D::Inf, Point3D::Inf);
    }
    
    bool InfiniteSphereSurfaceShape::preTransformed() const {
        return true;
    }
    
    bool InfiniteSphereSurfaceShape::intersectWithoutAlpha(const Ray &ray, const RaySegment &segment, SurfaceInteraction* si) const {
        if (!std::isinf(segment.distMax))
            return false;
        float theta, phi;
        ray.dir.toPolarYUp(&theta, &phi);
        *si = SurfaceInteraction(ray.time, // ---------------------------------- time
                                 INFINITY, // ---------------------------------- distance
                                 Point3D(ray.dir), // -------------------------- position in world coordinate
                                 -ray.dir, // ---------------------------------- geometric normal in world coordinate
                                 phi, theta, // -------------------------------- surface parameters
                                 TexCoord2D(phi / (2 * M_PI), theta / M_PI) // - texture coordinate
                                 );
        return true;
    }
    
    bool InfiniteSphereSurfaceShape::intersect(const Ray &ray, const RaySegment &segment, LightPathSampler &pathSampler, SurfaceInteraction *si) const {
        if (!std::isinf(segment.distMax))
            return false;
        float theta, phi;
        ray.dir.toPolarYUp(&theta, &phi);
        *si = SurfaceInteraction(ray.time, // ---------------------------------- time
                                 INFINITY, // ---------------------------------- distance
                                 Point3D(ray.dir), // -------------------------- position in world coordinate
                                 -ray.dir, // ---------------------------------- geometric normal in world coordinate
                                 phi, theta, // -------------------------------- surface parameters
                                 TexCoord2D(phi / (2 * M_PI), theta / M_PI) // - texture coordinate
                                 );
        return true;
    }
    
    float InfiniteSphereSurfaceShape::testVisibility(const Ray &ray, const RaySegment &segment) const {
        return std::isinf(segment.distMax) ? 0.0f : 1.0f;
    }
    
    void InfiniteSphereSurfaceShape::calculateSurfacePoint(const SurfaceInteraction &si, SurfacePoint *surfPt) const {
        float u, v;
        si.getSurfaceParameter(&u, &v);
        
        Vector3D texCoord0Dir = Vector3D(-std::cos(u), 0.0f, -std::sin(u));
        ReferenceFrame shadingFrame;
        shadingFrame.x = texCoord0Dir;
        shadingFrame.z = si.getGeometricNormal();
        shadingFrame.y = cross(shadingFrame.z, shadingFrame.x);
        
        *surfPt = SurfacePoint(si, true, shadingFrame, texCoord0Dir);
    }
    
    float InfiniteSphereSurfaceShape::area() const {
        return 4 * M_PI;// * Inf * Inf;
    }
    
    void InfiniteSphereSurfaceShape::sample(float u0, float u1, SurfacePoint *surfPt, float *areaPDF, DirectionType* posType) const {
        SLRAssert_NotImplemented();
    }
    
    float InfiniteSphereSurfaceShape::evaluateAreaPDF(const SurfacePoint &surfPt) const {
        SLRAssert_NotImplemented();
        return 0.0f;
    }    
}

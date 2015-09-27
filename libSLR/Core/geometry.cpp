//
//  geometry.cpp
//
//  Created by 渡部 心 on 2015/08/01.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "geometry.h"

#include "SurfaceObject.h"
#include "../Memory/ArenaAllocator.h"
#include "Transform.h"

namespace SLR {
    const float Ray::Epsilon = 0.0001f;
    
    void Intersection::getSurfacePoint(SurfacePoint* surfPt) {
        obj.top()->getSurfacePoint(*this, surfPt);
    }
    
    Vector3D SurfacePoint::getShadowDirection(const SurfacePoint &shadingPoint, float* dist2) {
        SLRAssert(shadingPoint.atInfinity == false, "Shading point must be in finite region.");
        if (atInfinity) {
            *dist2 = 1.0f;
            return normalize(p - Point3D::Zero);
        }
        else {
            Vector3D ret(p - shadingPoint.p);
            *dist2 = ret.sqLength();
            return ret / std::sqrt(*dist2);
        }
    }
    
    bool SurfacePoint::isEmitting() const {
        return obj->isEmitting();
    }
    
    SampledSpectrum SurfacePoint::emittance(const WavelengthSamples &wls) const {
        return obj->emittance(*this, wls);
    }
    
    float SurfacePoint::evaluateAreaPDF() const {
        return obj->evaluateAreaPDF(*this);
    }
    
    BSDF* SurfacePoint::createBSDF(const WavelengthSamples &wls, ArenaAllocator &mem) const {
        return obj->createBSDF(*this, wls, mem);
    }
    
    EDF* SurfacePoint::createEDF(const WavelengthSamples &wls, ArenaAllocator &mem) const {
        return obj->createEDF(*this, wls, mem);
    }
    
    SurfacePoint operator*(const StaticTransform &transform, const SurfacePoint &surfPt) {
        SurfacePoint ret;
        ret.p = transform * surfPt.p;
        ret.atInfinity = surfPt.atInfinity;
        ret.gNormal = normalize(transform * surfPt.gNormal);
        ret.u = surfPt.u;
        ret.v = surfPt.v;
        ret.texCoord = surfPt.texCoord;
        ret.texCoord0Dir = normalize(transform * surfPt.texCoord0Dir);
        ret.shadingFrame.x = normalize(transform * surfPt.shadingFrame.x);
        ret.shadingFrame.y = normalize(transform * surfPt.shadingFrame.y);
        ret.shadingFrame.z = normalize(transform * surfPt.shadingFrame.z);
        ret.obj = surfPt.obj;
        
        return ret;
    }    
}

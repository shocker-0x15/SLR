//
//  accelerator.h
//
//  Created by 渡部 心 on 2016/07/11.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef __SLR_accelerator__
#define __SLR_accelerator__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/geometry.h"

namespace SLR {
    class SLR_API Accelerator {
    public:
        virtual ~Accelerator() {}
        
        virtual float costForIntersect() const = 0;
        
        virtual BoundingBox3D bounds() const = 0;
        
        virtual bool intersectWithoutAlpha(const Ray &ray, const RaySegment &segment, SurfaceInteraction* si) const = 0;
        virtual bool intersect(const Ray &ray, const RaySegment &segment, LightPathSampler &pathSampler, SurfaceInteraction* si, const SurfaceObject** closestObject) const = 0;
        virtual float testVisibility(const Ray &ray, const RaySegment &segment) const = 0; 
        
        static bool traceTraverse;
        static std::string traceTraversePrefix;
    };
}

#endif /* __SLR_accelerator__ */

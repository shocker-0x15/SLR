//
//  Accelerator.h
//
//  Created by 渡部 心 on 2016/07/11.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef Accelerator_h
#define Accelerator_h

#include "../defines.h"
#include "../references.h"
#include "../Core/geometry.h"
#include "../Core/SurfaceObject.h"

namespace SLR {
    class SLR_API Accelerator {
    public:
        virtual ~Accelerator() {}
        
        virtual float costForIntersect() const = 0;
        
        virtual BoundingBox3D bounds() const = 0;
        
        virtual bool intersect(Ray &ray, SurfaceInteraction* si) const = 0;
        
        bool intersect(Ray &ray, SurfacePoint* surfPt) const {
            SurfaceInteraction si;
            if (!intersect(ray, &si))
                return false;
            si.getSurfacePoint(surfPt);
            return true;
        }
        
        static bool traceTraverse;
        static std::string traceTraversePrefix;
    };
}

#endif /* Accelerator_h */

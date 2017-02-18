//
//  Accelerator.h
//
//  Created by 渡部 心 on 2016/07/11.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef __SLR_Accelerator__
#define __SLR_Accelerator__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/geometry.h"
#include "../Core/SurfaceObject.h"

namespace SLR {
    class SLR_API Accelerator {
    public:
        virtual ~Accelerator() {}
        
        virtual float costForIntersect() const = 0;
        
        virtual BoundingBox3D bounds() const = 0;
        
        virtual bool intersect(Ray &ray, SurfaceInteraction* si, uint32_t* closestIndex) const = 0;
        
        static bool traceTraverse;
        static std::string traceTraversePrefix;
    };
}

#endif /* __SLR_Accelerator__ */

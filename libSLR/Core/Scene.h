//
//  Scene.h
//
//  Created by 渡部 心 on 2016/09/19.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef __SLR_Scene_h__
#define __SLR_Scene_h__

#include "../defines.h"
#include "../references.h"
#include "SurfaceObject.h"
#include "MediumObject.h"

namespace SLR {
    class SLR_API Scene {
        const SurfaceObjectAggregate* m_aggregate;
        const InfiniteSphereSurfaceObject* m_envSphere;
        Point3D m_worldCenter;
        float m_worldRadius;
        float m_worldDiscArea;
        const Camera* m_camera;
    public:
        Scene() { }
        
        void build(const SurfaceObjectAggregate* aggregate, const InfiniteSphereSurfaceObject* envSphere, const Camera* camera);
        
        const Camera* getCamera() const { return m_camera; }
        Point3D getWorldCenter() const { return m_worldCenter; }
        float getWorldRadius() const { return m_worldRadius; }
        float getWorldDiscArea() const { return m_worldDiscArea; }
        
        bool intersect(Ray &ray, Intersection* isect) const;
        bool testVisibility(const SurfacePoint &shdP, const SurfacePoint &lightP, float time) const;
        void selectLight(float u, Light* light, float* prob) const;
        float evaluateProb(const Light &light) const;
    };
}

#endif /* Scene_h */

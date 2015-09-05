//
//  Scene.h
//
//  Created by 渡部 心 on 2015/07/11.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__Scene__
#define __SLR__Scene__

#include "../defines.h"
#include "../references.h"
#include "../Core/geometry.h"

class Scene {
    InternalNodeRef m_rootNode;
    InfiniteSphereNodeRef m_envNode;
    Camera* m_mainCamera;
    SurfaceObjectAggregate* m_sceneObj;
    InfiniteSphereSurfaceObject* m_envSphere;
    ArenaAllocator* m_mem;
    
    Point3D m_worldCenter;
    float m_worldRadius;
    float m_worldDiscArea;
public:
    Scene();
    ~Scene();
    
    InternalNodeRef &rootNode() { return m_rootNode; };
    void setEnvNode(const InfiniteSphereNodeRef &node) { m_envNode = node; };
    
    void build();
    
    Camera* getMainCamera() const { return m_mainCamera; };
    const SurfaceObject* getSceneObject() const { return (const SurfaceObject*)m_sceneObj; };
    Point3D getWorldCenter() const { return m_worldCenter; };
    float getWorldRadius() const { return m_worldRadius; };
    float getWorldDiscArea() const { return m_worldDiscArea; };
    
    bool intersect(Ray &ray, Intersection* isect) const;
    bool testVisiblility(const SurfacePoint &shdP, const SurfacePoint &lightP, float time) const;
    void selectLight(float u, Light* light, float* prob) const;
    float evaluateProb(const Light &light) const;
};

#endif

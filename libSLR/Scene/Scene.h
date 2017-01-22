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
#include "nodes.h"
#include "SurfaceObject.h"
#include "MediumObject.h"

namespace SLR {
    class SLR_API Scene {
        Node* m_rootNode;
        InfiniteSphereNode* m_envNode;
        
        Allocator* m_sceneMem;
        SurfaceObjectAggregate* m_surfaceAggregate;
        MediumObjectAggregate* m_mediumAggregate;
        InfiniteSphereSurfaceObject* m_envSphere;
        Point3D m_worldCenter;
        float m_worldRadius;
        float m_worldDiscArea;
        Camera* m_camera;
    public:
        Scene(Node* rootNode) : m_rootNode(rootNode), m_envNode(nullptr) { }
        
        void setEnvironmentNode(InfiniteSphereNode* envNode) { m_envNode = envNode; }
        
        void build(Allocator* sceneMem);
        void destory();
        
        const Camera* getCamera() const { return m_camera; }
        Point3D getWorldCenter() const { return m_worldCenter; }
        float getWorldRadius() const { return m_worldRadius; }
        float getWorldDiscArea() const { return m_worldDiscArea; }
        
        bool intersect(Ray &ray, SurfaceInteraction* si) const;
        bool interact(Ray &ray, const WavelengthSamples &wls, LightPathSampler &pathSampler, ArenaAllocator &mem,
                      Interaction** interact, SampledSpectrum* medThroughput, bool* singleWavelength) const;
        bool testVisibility(const SurfacePoint &shdP, const SurfacePoint &lightP, float time) const;
        bool testVisibility(const InteractionPoint* shdP, const InteractionPoint* lightP, float time,
                            const WavelengthSamples &wls, LightPathSampler &pathSampler, SampledSpectrum* fractionalVisibility, bool* singleWavelength) const;
        void selectSurfaceLight(float u, SurfaceLight* light, float* prob) const;
        float evaluateSurfaceLightProb(const SurfaceLight &light) const;
        void selectLight(float u, ArenaAllocator &mem, Light** light, float* prob) const;
        float evaluateLightProb(const Light* light) const;
    };
}

#endif /* Scene_h */

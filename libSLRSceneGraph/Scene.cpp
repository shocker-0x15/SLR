//
//  Scene.cpp
//
//  Created by 渡部 心 on 2015/07/11.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "Scene.h"
#include "nodes.h"
#include <libSLR/Memory/ArenaAllocator.h>
#include <libSLR/Core/SurfaceObject.h>
#include <libSLR/Core/cameras.h>
#include <libSLR/Core/Transform.h>
#include "InfiniteSphereNode.h"

#include <libSLR/Core/Renderer.h>

namespace SLRSceneGraph {
    Scene::Scene() : m_envNode(nullptr) {
        m_rootNode = createShared<InternalNode>();
        m_rootNode->setName("root");
        m_rootNode->setTransform(createShared<SLR::StaticTransform>());
    };
    
    Scene::~Scene() {}
    
    void Scene::build(SLR::Scene** scene, SLR::ArenaAllocator &mem) {
        RenderingData renderingData;
        m_rootNode->getRenderingData(mem, nullptr, &renderingData);
        SLRAssert(renderingData.camera != nullptr, "Camera is not set.");
        
        SLR::SurfaceObjectAggregate* aggregate = mem.create<SLR::SurfaceObjectAggregate>(renderingData.surfObjs);
        SLR::InfiniteSphereSurfaceObject* envSphere = m_envNode ? m_envNode->getSurfaceObject() : nullptr;
        
        SLR::Camera* camera = renderingData.camera;
        camera->setTransform(renderingData.camTransform);
        *scene = mem.create<SLR::Scene>(aggregate, envSphere, camera);
    }
    
    RenderingContext::RenderingContext() {
        
    }
    
    RenderingContext::~RenderingContext() {
        
    }
    
    RenderingContext& RenderingContext::operator=(RenderingContext &&ctx) {
        renderer = std::move(ctx.renderer);
        width = ctx.width;
        height = ctx.height;
        timeStart = ctx.timeStart;
        timeEnd = ctx.timeEnd;
        rngSeed = ctx.rngSeed;
        
        return *this;
    }
}

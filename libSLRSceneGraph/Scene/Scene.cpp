//
//  Scene.cpp
//
//  Created by 渡部 心 on 2015/07/11.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "Scene.h"

#include <libSLR/Core/transform.h>
#include <libSLR/Core/renderer.h>
#include <libSLR/Scene/Scene.h>
#include "node.h"

namespace SLRSceneGraph {
    Scene::Scene() : m_envNode(nullptr) {
        m_rootNode = createShared<InternalNode>(createShared<SLR::StaticTransform>());
        m_rootNode->setName("root");
        m_raw = new SLR::Scene(m_rootNode->getRaw());
    };
    
    Scene::~Scene() {
        delete m_raw;
    }
    
    void Scene::setEnvironmentNode(const InfiniteSphereNodeRef &node) {
        m_envNode = node;
        m_raw->setEnvironmentNode((SLR::InfiniteSphereNode*)node->getRaw());
    }
    
    void Scene::prepareForRendering() {
        m_rootNode->prepareForRendering();
        if (m_envNode)
            m_envNode->prepareForRendering();
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
        brightness = ctx.brightness;
        rngSeed = ctx.rngSeed;
        
        return *this;
    }
}

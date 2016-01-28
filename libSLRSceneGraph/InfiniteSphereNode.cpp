//
//  InfiniteSphereNode.cpp
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "InfiniteSphereNode.h"
#include "surface_materials.hpp"
#include <libSLR/Core/SurfaceObject.h>

namespace SLRSceneGraph {
    InfiniteSphereNode::InfiniteSphereNode(const SceneWRef &scene, const SpectrumTextureRef &IBLTex, float scale) : m_scene(scene), m_ready(false) {
        m_emitter = createShared<IBLEmission>(m_scene, IBLTex, scale);
    }
    
    InfiniteSphereNode::~InfiniteSphereNode() {
        delete m_surfObj;
    }
    
    SLR::InfiniteSphereSurfaceObject* InfiniteSphereNode::getSurfaceObject() {
        if (!m_ready) {
            m_surfObj = new SLR::InfiniteSphereSurfaceObject(&m_sphere, (SLR::IBLEmission*)m_emitter->getRaw());
            m_ready = true;
        }
        return m_surfObj;
    }    
}

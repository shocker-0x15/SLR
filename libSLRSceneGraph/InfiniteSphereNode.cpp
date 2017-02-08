//
//  InfiniteSphereNode.cpp
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "InfiniteSphereNode.h"
#include "textures.hpp"
#include <libSLR/Scene/nodes.h>
#include "Scene.h"

namespace SLRSceneGraph {
    void InfiniteSphereNode::allocateRawData() {
        m_rawData = (SLR::Node*)malloc(sizeof(SLR::InfiniteSphereNode));
    }
    
    void InfiniteSphereNode::setupRawData() {
        SceneRef scene = m_scene.lock();
        new (m_rawData) SLR::InfiniteSphereNode(scene->getRaw(), m_IBLTex->getRaw(), m_scale);
        m_setup = true;
    }
    
    void InfiniteSphereNode::terminateRawData() {
        SLR::InfiniteSphereNode &raw = *(SLR::InfiniteSphereNode*)m_rawData;
        if (m_setup)
            raw.~InfiniteSphereNode();
        m_setup = false;
    }
    
    InfiniteSphereNode::InfiniteSphereNode(const SceneWRef &scene, const SpectrumTextureRef &IBLTex, float scale) :
    m_scene(scene), m_IBLTex(IBLTex), m_scale(scale) {
        allocateRawData();
    }
    
    void InfiniteSphereNode::prepareForRendering() {
        terminateRawData();
        setupRawData();
    }
}

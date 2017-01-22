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
    void InfiniteSphereNode::setupRawData() {
        SceneRef scene = m_scene.lock();
        m_rawData = new SLR::InfiniteSphereNode(scene->getRaw(), m_IBLTex->getRaw(), m_scale);
    }
    
    InfiniteSphereNode::InfiniteSphereNode(const SceneWRef &scene, const SpectrumTextureRef &IBLTex, float scale) :
    m_scene(scene), m_IBLTex(IBLTex), m_scale(scale) {
        setupRawData();
    }
    
    void InfiniteSphereNode::prepareForRendering() {
        SLR::InfiniteSphereNode &raw = *(SLR::InfiniteSphereNode*)m_rawData;
        SceneRef scene = m_scene.lock();
        new (&raw) SLR::InfiniteSphereNode(scene->getRaw(), m_IBLTex->getRaw(), m_scale);
    }
}

//
//  InfinitesimalPointNode.cpp
//
//  Created by 渡部 心 on 2017/02/17.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "InfinitesimalPointNode.h"
#include "textures.hpp"
#include "surface_materials.hpp"
#include "medium_nodes.h"
#include <libSLR/Core/Transform.h>
#include <libSLR/Core/SurfaceObject.h>
#include <libSLR/Scene/nodes.h>

namespace SLRSceneGraph {
    void InfinitesimalPointNode::allocateRawData() {
        m_rawData = (SLR::Node*)malloc(sizeof(SLR::InfinitesimalPointNode));
    }
    
    void InfinitesimalPointNode::setupRawData() {
        new (m_rawData) SLR::InfinitesimalPointNode(m_position, m_direction, m_material->getRaw());
        m_setup = true;
    }
    
    void InfinitesimalPointNode::terminateRawData() {
        SLR::InfinitesimalPointNode &raw = *(SLR::InfinitesimalPointNode*)m_rawData;
        if (m_setup)
            raw.~InfinitesimalPointNode();
        m_setup = false;
    }
    
    InfinitesimalPointNode::InfinitesimalPointNode(const SLR::Point3D &p, const SLR::Vector3D &d, const SurfaceMaterialRef &mat) :
    m_position(p), m_direction(d), m_material(mat) {
        allocateRawData();
    }
    
    NodeRef InfinitesimalPointNode::copy() const {
        NodeRef ret = createShared<InfinitesimalPointNode>(m_position, m_direction, m_material);
        return ret;
    }
    
    void InfinitesimalPointNode::prepareForRendering() {
        terminateRawData();
        setupRawData();
    }
}

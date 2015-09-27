//
//  InfiniteSphereNode.cpp
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "InfiniteSphereNode.h"
#include <libSLR/Core/textures.h>
#include <libSLR/Core/surface_material.h>
#include <libSLR/Core/distributions.h>
#include <libSLR/Core/SurfaceObject.h>
#include <libSLR/SurfaceMaterials/IBLEmission.h>

namespace SLRSceneGraph {
    InfiniteSphereNode::InfiniteSphereNode(const SLR::Scene* scene, const SLR::SpectrumTextureRef &IBLTex, float scale) : m_scene(scene) {
        m_dist = IBLTex->createIBLImportanceMap();
        //    m_dist->exportBMP("distribution.bmp", 4);
        m_emitter = createShared<SLR::IBLEmission>(m_scene, IBLTex, scale);
        m_surfMat = createShared<SLR::EmitterSurfaceMaterial>(nullptr, m_emitter);
    }
    
    InfiniteSphereNode::~InfiniteSphereNode() {
        delete m_dist;
        delete m_surfObj;
    }
    
    SLR::InfiniteSphereSurfaceObject* InfiniteSphereNode::getSurfaceObject() {
        if (!m_ready) {
            m_surfObj = new SLR::InfiniteSphereSurfaceObject(&m_sphere, m_surfMat.get(), m_dist);
            m_ready = true;
        }
        return m_surfObj;
    }    
}

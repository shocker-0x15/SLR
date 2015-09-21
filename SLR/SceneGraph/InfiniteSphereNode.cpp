//
//  InfiniteSphereNode.cpp
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "InfiniteSphereNode.h"
#include "../Core/textures.h"
#include "../Core/surface_material.h"
#include "../Core/distributions.h"
#include "../Core/SurfaceObject.h"
#include "../SurfaceMaterials/IBLEmission.h"

InfiniteSphereNode::InfiniteSphereNode(const Scene* scene, const SpectrumTextureRef &IBLTex, const SampledSpectrum &scale) : m_scene(scene) {
    m_dist = IBLTex->createIBLImportanceMap();
//    m_dist->exportBMP("distribution.bmp", 4);
    m_emitter = createShared<IBLEmission>(m_scene, IBLTex, scale);
    m_surfMat = createShared<EmitterSurfaceMaterial>(nullptr, m_emitter);
}

InfiniteSphereNode::~InfiniteSphereNode() {
    delete m_dist;
    delete m_surfObj;
}

InfiniteSphereSurfaceObject* InfiniteSphereNode::getSurfaceObject() {
    if (!m_ready) {
        m_surfObj = new InfiniteSphereSurfaceObject(&m_sphere, m_surfMat.get(), m_dist);
        m_ready = true;
    }
    return m_surfObj;
}

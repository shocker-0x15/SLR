//
//  TriangleMeshNode.cpp
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "TriangleMeshNode.h"
#include "../Core/SurfaceObject.h"

uint64_t TriangleMeshNode::addVertex(const Vertex &v) {
    m_verticesOrg.push_back(v);
    m_vertices.push_back(v);
    return m_vertices.size() - 1;
}

void TriangleMeshNode::addTriangle(uint64_t v0, uint64_t v1, uint64_t v2, const SurfaceMaterialRef &mat, const Normal3DTextureRef &normalMap) {
    Triangle tri{this, v0, v1, v2, nullptr};
    m_triangles.push_back(tri);
    m_materials.push_back(mat);
    m_normalMaps.push_back(normalMap);
}

void TriangleMeshNode::applyTransform(const StaticTransform &t) {
    for (int i = 0; i < m_vertices.size(); ++i) {
        Vertex &vo = m_verticesOrg[i];
        Vertex &v = m_vertices[i];
        v.position = t * vo.position;
        v.normal = normalize(t * vo.normal);
        v.tangent = normalize(t * vo.tangent);
    }
}

void TriangleMeshNode::createSurfaceObjects() {
    //    for (int i = 0; i < m_vertices.size(); ++i) {
    //        Vertex &vo = m_verticesOrg[i];
    //        Vertex &v = m_vertices[i];
    //        v.position = t * vo.position;
    //        v.normal = normalize(t * vo.normal);
    //        v.tangent = normalize(t * vo.tangent);
    //    }
    
    m_numRefinedObjs = m_triangles.size();
    m_refinedObjs = new SingleSurfaceObject*[m_numRefinedObjs];
    for (int i = 0; i < m_triangles.size(); ++i) {
        if (m_normalMaps[i])
            m_refinedObjs[i] = new BumpSingleSurfaceObject(&m_triangles[i], m_materials[i].get(), m_normalMaps[i].get());
        else
            m_refinedObjs[i] = new SingleSurfaceObject(&m_triangles[i], m_materials[i].get());
    }
}

//
//  TriangleMeshNode.cpp
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "TriangleMeshNode.h"
#include "textures.hpp"
#include "surface_materials.hpp"
#include <libSLR/Core/Transform.h>
#include <libSLR/Core/SurfaceObject.h>
#include <libSLR/Surface/TriangleMesh.h>

namespace SLRSceneGraph {
    TriangleMeshNode::~TriangleMeshNode() {
        if (m_trianglesForRendering)
            delete[] m_trianglesForRendering;
    }

    uint64_t TriangleMeshNode::addVertex(const SLR::Vertex &v) {
        m_vertices.push_back(v);
        return m_vertices.size() - 1;
    }
    
    void TriangleMeshNode::addTriangle(uint64_t v0, uint64_t v1, uint64_t v2, const SurfaceMaterialRef &mat, const Normal3DTextureRef &normalMap) {
        Triangle tri{{v0, v1, v2}};
        m_triangles.push_back(tri);
        m_materials.push_back(mat);
        m_normalMaps.push_back(normalMap);
    }
    
    void TriangleMeshNode::applyTransform(const SLR::StaticTransform &t) {
        for (int i = 0; i < m_vertices.size(); ++i) {
            Vertex &v = m_vertices[i];
            v.position = t * v.position;
            v.normal = normalize(t * v.normal);
            v.tangent = normalize(t * v.tangent);
        }
    }
    
    void TriangleMeshNode::createSurfaceObjects() {
        m_numRefinedObjs = m_triangles.size();
        m_refinedObjs = new SLR::SingleSurfaceObject*[m_numRefinedObjs];
        m_trianglesForRendering = new SLR::Triangle[m_numRefinedObjs];
        for (int i = 0; i < m_triangles.size(); ++i) {
            Triangle &tri = m_triangles[i];
            SLR::Triangle &triR = m_trianglesForRendering[i];
            new (&triR) SLR::Triangle(&m_vertices[tri.vIdx[0]], &m_vertices[tri.vIdx[1]], &m_vertices[tri.vIdx[2]], nullptr);
            if (m_normalMaps[i])
                m_refinedObjs[i] = new SLR::BumpSingleSurfaceObject(&triR, m_materials[i]->getRaw(), m_normalMaps[i]->getRaw());
            else
                m_refinedObjs[i] = new SLR::SingleSurfaceObject(&triR, m_materials[i]->getRaw());
        }
    }
}

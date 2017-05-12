//
//  TriangleMeshNode.cpp
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "TriangleMeshNode.h"

#include <libSLR/Core/transform.h>
#include <libSLR/Core/surface_object.h>
#include <libSLR/SurfaceShape/TriangleSurfaceShape.h>
#include <libSLR/Scene/TriangleMeshNode.h>
#include "../textures.h"
#include "../surface_materials.h"
#include "medium_nodes.h"

namespace SLRSceneGraph {
    void TriangleMeshNode::allocateRawData() {
        m_rawData = (SLR::Node*)malloc(sizeof(SLR::TriangleMeshNode));
    }
    
    void TriangleMeshNode::setupRawData() {
        new (m_rawData) SLR::TriangleMeshNode((uint32_t)m_vertices.size(), (uint32_t)m_matGroups.size(), m_onlyForBoundary, m_axisForRadialTangent);
        
        SLR::TriangleMeshNode &raw = *(SLR::TriangleMeshNode*)getRaw();
        
        uint32_t numVertices = (uint32_t)m_vertices.size();
        SLR::Vertex* vertices = raw.getVertexArray();
        for (int i = 0; i < numVertices; ++i) {
            SLR::Vertex &vertex = vertices[i];
            const Vertex &srcVtx = m_vertices[i];
            vertex = SLR::Vertex(srcVtx.position, srcVtx.normal, srcVtx.tangent, srcVtx.texCoord);
        }
        
        uint32_t numMatGroups = (uint32_t)m_matGroups.size();
        SLR::MaterialGroupInTriangleMesh* matGroups = raw.getMaterialGroupArray();
        for (int i = 0; i < numMatGroups; ++i) {
            SLR::MaterialGroupInTriangleMesh &matGroup = matGroups[i];
            const MaterialGroup &srcMatGroup = m_matGroups[i];
            matGroup.material = srcMatGroup.material->getRaw();
            matGroup.normalMap = srcMatGroup.normalMap ? srcMatGroup.normalMap->getRaw() : nullptr;
            matGroup.alphaMap = srcMatGroup.alphaMap ? srcMatGroup.alphaMap->getRaw() : nullptr;
            
            uint32_t numTriangles = (uint32_t)srcMatGroup.triangles.size();
            SLR::Vertex** vertexReferences = new SLR::Vertex*[3 * numTriangles];
            for (int j = 0; j < numTriangles; ++j) {
                const Triangle &srcTri = srcMatGroup.triangles[j];
                vertexReferences[3 * j + 0] = &vertices[srcTri.vIdx[0]];
                vertexReferences[3 * j + 1] = &vertices[srcTri.vIdx[1]];
                vertexReferences[3 * j + 2] = &vertices[srcTri.vIdx[2]];
            }
            
            std::unique_ptr<SLR::Vertex*[]> vertexReferencesHolder(vertexReferences);
            matGroup.setTriangles(vertexReferencesHolder, numTriangles);
        }
        
        if (m_enclosedMediumNode) {
            m_enclosedMediumNode->prepareForRendering();
            raw.setInternalMedium((SLR::MediumNode*)m_enclosedMediumNode->getRaw());
        }
        
        m_setup = true;
    }
    
    void TriangleMeshNode::terminateRawData() {
        SLR::TriangleMeshNode &raw = *(SLR::TriangleMeshNode*)m_rawData;
        if (m_setup)
            raw.~TriangleMeshNode();
        m_setup = false;
    }
    
    TriangleMeshNode::TriangleMeshNode() : m_onlyForBoundary(false) {
        allocateRawData();
    }

    uint64_t TriangleMeshNode::addVertex(const SLR::Vertex &v) {
        m_vertices.push_back(v);
        return m_vertices.size() - 1;
    }
    
    void TriangleMeshNode::addMaterialGroup(const SurfaceMaterialRef &mat, const NormalTextureRef &normalMap, const FloatTextureRef &alphaMap,
                                            const std::vector<Triangle> &&triangles) {
        m_matGroups.emplace_back();
        MaterialGroup &matGroup = m_matGroups.back();
        matGroup.material = mat;
        matGroup.normalMap = normalMap;
        matGroup.alphaMap = alphaMap;
        matGroup.triangles = triangles;
        
        for (int i = 0; i < matGroup.triangles.size(); ++i) {
            const Triangle &tri = matGroup.triangles[i];
            const SLR::Vertex &vtx0 = m_vertices[tri.vIdx[0]];
            const SLR::Vertex &vtx1 = m_vertices[tri.vIdx[1]];
            const SLR::Vertex &vtx2 = m_vertices[tri.vIdx[2]];
            SLR::Vector3D gNormal = SLR::cross(vtx1.position - vtx0.position, vtx2.position - vtx0.position);
            bool faceAway = false;
            faceAway |= SLR::dot(vtx0.normal, gNormal) <= 0;
            faceAway |= SLR::dot(vtx1.normal, gNormal) <= 0;
            faceAway |= SLR::dot(vtx2.normal, gNormal) <= 0;
            if (faceAway)
                printf("Warning: triangle [%llu, %llu, %llu] has a shading normal too faced away from a geometric normal.\n", tri.vIdx[0], tri.vIdx[1], tri.vIdx[2]);
        }
    }
    
    NodeRef TriangleMeshNode::copy() const {
        TriangleMeshNodeRef ret = createShared<TriangleMeshNode>();
        ret->m_vertices = m_vertices;
        ret->m_matGroups = m_matGroups;
        return ret;
    }
    
    void TriangleMeshNode::applyTransform(const SLR::StaticTransform &t) {
        for (int i = 0; i < m_vertices.size(); ++i) {
            Vertex &v = m_vertices[i];
            v.position = t * v.position;
            v.normal = normalize(t * v.normal);
            v.tangent = normalize(t * v.tangent);
        }
    }
    
    void TriangleMeshNode::prepareForRendering() {
        terminateRawData();
        setupRawData();
    }
}

//
//  TriangleMeshNode.cpp
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "TriangleMeshNode.h"
#include "textures.hpp"
#include "surface_materials.hpp"
#include "medium_nodes.h"
#include <libSLR/Core/Transform.h>
#include <libSLR/Core/SurfaceObject.h>
#include <libSLR/Surface/TriangleMesh.h>
#include <libSLR/Scene/TriangleMeshNode.h>

namespace SLRSceneGraph {
    void TriangleMeshNode::allocateRawData() {
        m_rawData = (SLR::Node*)malloc(sizeof(SLR::TriangleMeshNode));
    }
    
    void TriangleMeshNode::setupRawData() {
        new (m_rawData) SLR::TriangleMeshNode(m_onlyForBoundary);
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
    
    void TriangleMeshNode::addTriangles(const SurfaceMaterialRef &mat, const Normal3DTextureRef &normalMap, const FloatTextureRef &alphaMap,
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
        SLR::TriangleMeshNode &raw = *(SLR::TriangleMeshNode*)getRaw();
        
        uint32_t numVertices = (uint32_t)m_vertices.size();
        SLR::Vertex* vertices = new SLR::Vertex[numVertices];
        for (int i = 0; i < numVertices; ++i) {
            SLR::Vertex &vertex = vertices[i];
            const Vertex &srcVtx = m_vertices[i];
            vertex = SLR::Vertex(srcVtx.position, srcVtx.normal, srcVtx.tangent, srcVtx.texCoord);
        }
        
        uint32_t numMatGroups = (uint32_t)m_matGroups.size();
        SLR::TriangleMeshNode::MaterialGroup* matGroups = new SLR::TriangleMeshNode::MaterialGroup[numMatGroups];
        for (int i = 0; i < numMatGroups; ++i) {
            SLR::TriangleMeshNode::MaterialGroup &matGroup = matGroups[i];
            const MaterialGroup &srcMatGroup = m_matGroups[i];
            const SLR::SurfaceMaterial* surfMat = srcMatGroup.material->getRaw();
            const SLR::Normal3DTexture* normalMap = srcMatGroup.normalMap ? srcMatGroup.normalMap->getRaw() : nullptr;
            const SLR::FloatTexture* alphaMap = srcMatGroup.alphaMap ? srcMatGroup.alphaMap->getRaw() : nullptr;
            new (&matGroup) SLR::TriangleMeshNode::MaterialGroup(surfMat, normalMap, alphaMap);
            
            uint32_t numTriangles = (uint32_t)srcMatGroup.triangles.size();
            SLR::Triangle* triangles = new SLR::Triangle[numTriangles];
            for (int j = 0; j < numTriangles; ++j) {
                SLR::Triangle &triangle = triangles[j];
                const Triangle &srcTri = srcMatGroup.triangles[j];
                triangle = SLR::Triangle(&vertices[srcTri.vIdx[0]], &vertices[srcTri.vIdx[1]], &vertices[srcTri.vIdx[2]], alphaMap);
            }
            
            std::unique_ptr<SLR::Triangle[]> trianglesHolder(triangles);
            matGroup.setTriangles(trianglesHolder, numTriangles);
        }
        std::unique_ptr<SLR::Vertex[]> verticesHolder(vertices);
        std::unique_ptr<SLR::TriangleMeshNode::MaterialGroup[]> matGroupsHolder(matGroups);
        raw.setVertices(verticesHolder, numVertices);
        raw.setMaterialGroups(matGroupsHolder, numMatGroups);
        
        if (m_enclosedMediumNode) {
            m_enclosedMediumNode->prepareForRendering();
            raw.setInternalMedium((SLR::MediumNode*)m_enclosedMediumNode->getRaw());
        }
    }
}

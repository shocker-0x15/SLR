//
//  TriangleMeshNode.cpp
//
//  Created by 渡部 心 on 2017/01/05.
//  Copyright © 2017年 渡部 心. All rights reserved.
//

#include "TriangleMeshNode.h"
#include "../MemoryAllocators/ArenaAllocator.h"
#include "../Core/SurfaceObject.h"
#include "../Core/MediumObject.h"
#include "../Core/Transform.h"
#include "../Scene/medium_nodes.h"

namespace SLR {
    void TriangleMeshNode::MaterialGroup::setTriangles(std::unique_ptr<TriangleSurfaceShape[]> &triangles, uint32_t numTriangles) {
        m_triangles = std::move(triangles);
        m_numTriangles = numTriangles;
    }
    
    
    
    void TriangleMeshNode::setVertices(std::unique_ptr<Vertex[]> &vertices, uint32_t numVertices) {
        m_vertices = std::move(vertices);
        m_numVertices = numVertices;
    }
    
    void TriangleMeshNode::setMaterialGroups(std::unique_ptr<MaterialGroup[]> &matGroups, uint32_t numMatGroups) {
        m_matGroups = std::move(matGroups);
        m_numMatGroups = numMatGroups;
    }
    
    void TriangleMeshNode::createRenderingData(Allocator* mem, const Transform* subTF, RenderingData* data) {
        uint32_t numObjects = 0;
        for (int i = 0; i < m_numMatGroups; ++i)
            numObjects += m_matGroups[i].m_numTriangles;
        size_t objBaseIdx = data->surfObjs.size();
        m_objs.resize(numObjects);
        if (!m_onlyForBoundary)
            data->surfObjs.resize(objBaseIdx + numObjects);
        
        // apply transform
        StaticTransform transform;
        if (subTF) {
            SLRAssert(subTF->isStatic(), "Transformation given to TriangleMeshNode must be static.");
            subTF->sample(0.0f, &transform);
        }
        if (!transform.isIdentity()) {
            for (int i = 0; i < m_numVertices; ++i) {
                Vertex &v = m_vertices[i];
                v.position = transform * v.position;
                v.normal = normalize(transform * v.normal);
                v.tangent = normalize(transform * v.tangent);
                v.texCoord = v.texCoord;
            }
        }
        
        // create surface objects
        uint32_t triBaseIdx = 0;
        for (int i = 0; i < m_numMatGroups; ++i) {
            const MaterialGroup &matGroup = m_matGroups[i];
            
            const SurfaceMaterial* mat = matGroup.m_material;
            const NormalTexture* normalMap = matGroup.m_normalMap;
//            const FloatTexture* alphaMap = matGroup.m_alphaMap;
            
            for (int tIdx = 0; tIdx < matGroup.m_numTriangles; ++tIdx) {
                TriangleSurfaceShape &tri = matGroup.m_triangles[tIdx];

                SurfaceObject* obj;
                if (normalMap)
                    obj = mem->create<BumpSingleSurfaceObject>(&tri, mat, normalMap);
                else
                    obj = mem->create<SingleSurfaceObject>(&tri, mat);
                m_objs[triBaseIdx + tIdx] = obj;
                if (!m_onlyForBoundary)
                    data->surfObjs[objBaseIdx + triBaseIdx + tIdx] = obj;
            }
            triBaseIdx += matGroup.m_numTriangles;
        }
        
        // create an enclosed medium object
        if (m_enclosedMediumNode) {
            // create a dedicated acceleration structure and enclosed medium object.
            if (subTF)
                m_mediumTransform = subTF->copy(mem);
            
            RenderingData subData(nullptr);
            m_enclosedMediumNode->createRenderingData(mem, nullptr, &subData);
            m_boundarySurfObj = mem->create<SurfaceObjectAggregate>(m_objs);
            m_enclosedMedObj = mem->create<EnclosedMediumObject>(subData.medObjs[0], m_boundarySurfObj,
                                                                 m_mediumTransform ? *(StaticTransform*)m_mediumTransform : StaticTransform());
            if (subTF && !transform.isIdentity()) {
                m_TFMedObj = mem->create<TransformedMediumObject>(m_enclosedMedObj, m_mediumTransform);
                data->medObjs.push_back(m_TFMedObj);
            }
            else {
                data->medObjs.push_back(m_enclosedMedObj);
            }
        }
    }
    
    void TriangleMeshNode::destroyRenderingData(Allocator* mem) {
        if (m_enclosedMediumNode) {
            if (m_TFMedObj)
                mem->destroy(m_TFMedObj);
            m_TFMedObj = nullptr;
            mem->destroy(m_enclosedMedObj);
            m_enclosedMedObj = nullptr;
            mem->destroy(m_boundarySurfObj);
            m_boundarySurfObj = nullptr;
            m_enclosedMediumNode->destroyRenderingData(mem);
            
            if (m_mediumTransform)
                mem->destroy(m_mediumTransform);
            m_mediumTransform = nullptr;
        }
        
        for (int i = (int)m_objs.size() - 1; i >= 0; --i)
            mem->destroy(m_objs[i]);
        m_objs.clear();
    }
}

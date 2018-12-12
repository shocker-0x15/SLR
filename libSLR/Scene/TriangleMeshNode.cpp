//
//  TriangleMeshNode.cpp
//
//  Created by 渡部 心 on 2017/01/05.
//  Copyright c 2017年 渡部 心. All rights reserved.
//

#include "TriangleMeshNode.h"

#include "../MemoryAllocators/ArenaAllocator.h"
#include "../Core/transform.h"
#include "../Core/surface_object.h"
#include "../Core/medium_object.h"
#include "../SurfaceShape/TriangleSurfaceShape.h"
#include "../Scene/medium_nodes.h"

namespace SLR {
    void MaterialGroupInTriangleMesh::setTriangles(std::unique_ptr<Vertex*[]> &_vertexReferences, uint32_t _numTriangles) {
        vertexReferences = std::move(_vertexReferences);
        numTriangles = _numTriangles;
        triangles = new TriangleSurfaceShape[_numTriangles];
        for (int i = 0; i < numTriangles; ++i) {
            new (triangles + i) TriangleSurfaceShape(this, 3 * i);
        }
    }
    
    MaterialGroupInTriangleMesh::~MaterialGroupInTriangleMesh() {
        if (triangles)
            delete[] triangles;
        triangles = nullptr;
    }

    
    
    TriangleMeshNode::TriangleMeshNode(uint32_t numVertices, uint32_t numMatGroups, bool onlyForBoundary, int8_t axisForRadialTangent) : 
    m_numVertices(numVertices), m_numMatGroups(numMatGroups), m_onlyForBoundary(onlyForBoundary), m_axisForRadialTangent(axisForRadialTangent) {
        m_vertices = new Vertex[m_numVertices];
        m_matGroups = new MaterialGroupInTriangleMesh[m_numMatGroups];
        for (int i = 0; i < m_numMatGroups; ++i)
            m_matGroups[i].parent = this;
    }
    
    TriangleMeshNode::~TriangleMeshNode() {
        if (m_matGroups)
            delete[] m_matGroups;
        if (m_vertices)
            delete[] m_vertices;
        m_matGroups = nullptr;
        m_vertices = nullptr;
    }
    
    void TriangleMeshNode::createRenderingData(Allocator* mem, const Transform* subTF, RenderingData* data) {
        uint32_t numObjects = 0;
        for (int i = 0; i < m_numMatGroups; ++i)
            numObjects += m_matGroups[i].numTriangles;
        size_t objBaseIdx = data->surfObjs.size();
        m_objs.resize(numObjects);
        if (!m_onlyForBoundary)
            data->surfObjs.resize(objBaseIdx + numObjects);
        
        // apply transform
        if (subTF) {
            SLRAssert(subTF->isStatic(), "Transformation given to TriangleMeshNode must be static.");
            subTF->sample(0.0f, &m_appliedTransform);
        }
        m_appliedTFIsIdentity = m_appliedTransform.isIdentity();
        if (!m_appliedTFIsIdentity) {
            for (int i = 0; i < m_numVertices; ++i) {
                Vertex &v = m_vertices[i];
                v.position = m_appliedTransform * v.position;
                v.normal = normalize(m_appliedTransform * v.normal);
                v.tangent = normalize(m_appliedTransform * v.tangent);
                v.texCoord = v.texCoord;
            }
        }
        
        // create surface objects
        uint32_t triBaseIdx = 0;
        for (int i = 0; i < m_numMatGroups; ++i) {
            const MaterialGroupInTriangleMesh &matGroup = m_matGroups[i];
            
            const SurfaceMaterial* mat = matGroup.material;
            const NormalTexture* normalMap = matGroup.normalMap;
//            const FloatTexture* alphaMap = matGroup.alphaMap;
            
            for (int tIdx = 0; tIdx < matGroup.numTriangles; ++tIdx) {
                TriangleSurfaceShape &tri = matGroup.triangles[tIdx];

                SurfaceObject* obj;
                if (normalMap)
                    obj = mem->create<BumpSingleSurfaceObject>(&tri, mat, normalMap);
                else
                    obj = mem->create<SingleSurfaceObject>(&tri, mat);
                m_objs[triBaseIdx + tIdx] = obj;
                if (!m_onlyForBoundary)
                    data->surfObjs[objBaseIdx + triBaseIdx + tIdx] = obj;
            }
            triBaseIdx += matGroup.numTriangles;
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
            if (subTF && !m_appliedTFIsIdentity) {
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

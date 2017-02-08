//
//  TriangleMeshNode.cpp
//
//  Created by 渡部 心 on 2017/01/05.
//  Copyright © 2017年 渡部 心. All rights reserved.
//

#include "TriangleMeshNode.h"
#include "../Memory/ArenaAllocator.h"
#include "../Core/SurfaceObject.h"
#include "../Core/Transform.h"

namespace SLR {
    void TriangleMeshNode::MaterialGroup::setTriangles(std::unique_ptr<Triangle[]> &triangles, uint32_t numTriangles) {
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
        data->surfObjs.resize(objBaseIdx + numObjects);
        
        StaticTransform transform;
        if (subTF) {
            SLRAssert(subTF->isStatic(), "Transformation given to TriangleMeshNode must be static.");
            subTF->sample(0.0f, &transform);
        }
        for (int i = 0; i < m_numVertices; ++i) {
            Vertex &v = m_vertices[i];
            v.position = transform * v.position;
            v.normal = normalize(transform * v.normal);
            v.tangent = normalize(transform * v.tangent);
            v.texCoord = v.texCoord;
        }
        
        uint32_t triBaseIdx = 0;
        for (int i = 0; i < m_numMatGroups; ++i) {
            const MaterialGroup &matGroup = m_matGroups[i];
            
            const SurfaceMaterial* mat = matGroup.m_material;
            const Normal3DTexture* normalMap = matGroup.m_normalMap;
//            const FloatTexture* alphaMap = matGroup.m_alphaMap;
            
            for (int tIdx = 0; tIdx < matGroup.m_numTriangles; ++tIdx) {
                Triangle &tri = matGroup.m_triangles[tIdx];

                SurfaceObject* obj;
                if (normalMap)
                    obj = mem->create<BumpSingleSurfaceObject>(&tri, mat, normalMap);
                else
                    obj = mem->create<SingleSurfaceObject>(&tri, mat);
                m_objs[triBaseIdx + tIdx] = obj;
                data->surfObjs[objBaseIdx + triBaseIdx + tIdx] = obj;
            }
            triBaseIdx += matGroup.m_numTriangles;
        }
    }
    
    void TriangleMeshNode::destroyRenderingData(Allocator* mem) {
        for (int i = (int)m_objs.size() - 1; i >= 0; --i)
            mem->destroy(m_objs[i]);
    }
}

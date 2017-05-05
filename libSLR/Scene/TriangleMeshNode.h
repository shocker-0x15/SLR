//
//  TriangleMeshNode.h
//
//  Created by 渡部 心 on 2017/01/05.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef __SLR_TriangleMeshNode__
#define __SLR_TriangleMeshNode__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/geometry.h"
#include "node.h"

namespace SLR {
    struct SLR_API MaterialGroupInTriangleMesh {
        const TriangleMeshNode* parent;
        const SurfaceMaterial* material;
        const NormalTexture* normalMap;
        const FloatTexture* alphaMap;
        std::unique_ptr<Vertex*[]> vertexReferences;
        TriangleSurfaceShape* triangles;
        uint32_t numTriangles;
        
        MaterialGroupInTriangleMesh() :
        material(nullptr), normalMap(nullptr), alphaMap(nullptr),
        numTriangles(0) {
        }
        ~MaterialGroupInTriangleMesh();
        
        void setTriangles(std::unique_ptr<Vertex*[]> &vertexReferences, uint32_t numTriangles);
    };
    
    
    
    class SLR_API TriangleMeshNode : public SurfaceNode {
        Vertex* m_vertices;
        uint32_t m_numVertices;
        MaterialGroupInTriangleMesh* m_matGroups;
        uint32_t m_numMatGroups;
        bool m_onlyForBoundary;
        
        std::vector<SurfaceObject*> m_objs;
    public:
        TriangleMeshNode(uint32_t numVertices, uint32_t numMatGroups, bool onlyForBoundary);
        ~TriangleMeshNode();
        
        Vertex* getVertexArray() {
            return m_vertices;
        }
        MaterialGroupInTriangleMesh* getMaterialGroupArray() {
            return m_matGroups;
        }
        
        bool isDirectlyTransformable() const override { return true; }
        void createRenderingData(Allocator* mem, const Transform* subTF, RenderingData *data) override;
        void destroyRenderingData(Allocator* mem) override;
    };
}

#endif /* __SLR_TriangleMeshNode__ */

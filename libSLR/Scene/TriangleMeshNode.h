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
        triangles(nullptr), numTriangles(0) {
        }
        ~MaterialGroupInTriangleMesh();
        
        void setTriangles(std::unique_ptr<Vertex*[]> &vertexReferences, uint32_t numTriangles);
    };
    
    
    
    class SLR_API TriangleMeshNode : public SurfaceNode {        
        Vertex* m_vertices;
        uint32_t m_numVertices;
        MaterialGroupInTriangleMesh* m_matGroups;
        uint32_t m_numMatGroups;
        // Should these be the member of SurfaceNode?
        bool m_onlyForBoundary;
        int8_t m_axisForRadialTangent; // 0:X, 1:Y, 2:Z
        bool m_appliedTFIsIdentity;
        StaticTransform m_appliedTransform;
        
        std::vector<SurfaceObject*> m_objs;
    public:
        TriangleMeshNode(uint32_t numVertices, uint32_t numMatGroups, bool onlyForBoundary, int8_t axisForRadialTangent);
        ~TriangleMeshNode();
        
        Vertex* getVertexArray() {
            return m_vertices;
        }
        MaterialGroupInTriangleMesh* getMaterialGroupArray() {
            return m_matGroups;
        }
        int8_t getAxisForRadialTangent() const {
            return m_axisForRadialTangent;
        }
        bool appliedTransformIsIdentity() const {
            return m_appliedTFIsIdentity;
        }
        const StaticTransform &getAppliedTransform() const {
            return m_appliedTransform;
        }
        
        bool isDirectlyTransformable() const override { return true; }
        void createRenderingData(Allocator* mem, const Transform* subTF, RenderingData *data) override;
        void destroyRenderingData(Allocator* mem) override;
    };
}

#endif /* __SLR_TriangleMeshNode__ */

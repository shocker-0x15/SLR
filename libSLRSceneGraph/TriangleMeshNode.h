//
//  TriangleMeshNode.h
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLRSceneGraph__TriangleMeshNode__
#define __SLRSceneGraph__TriangleMeshNode__

#include <libSLR/defines.h>
#include "references.h"
#include "nodes.h"
#include <libSLR/Core/geometry.h>

namespace SLRSceneGraph {
    // There is a possibility to define a unique vertex format suitable for editing.
    // Current implementation uses the vertex for rendering as is.
    typedef SLR::Vertex Vertex;
    
    struct SLR_SCENEGRAPH_API Triangle {
        uint64_t vIdx[3];
        Triangle() {}
        Triangle(uint64_t v0, uint64_t v1, uint64_t v2) : vIdx{v0, v1, v2} {}
    };
    
    class SLR_SCENEGRAPH_API TriangleMeshNode : public SurfaceObjectNode {
    public:
        struct MaterialGroup {
            SurfaceMaterialRef material;
            Normal3DTextureRef normalMap;
            FloatTextureRef alphaMap;
            std::vector<Triangle> triangles;
        };
    private:
        SLR::Triangle* m_trianglesForRendering;
        Vertex* m_verticesForRendering;
        
        std::vector<Vertex> m_vertices;
        std::vector<MaterialGroup> m_matGroups;
    public:
        TriangleMeshNode() : SurfaceObjectNode(), m_trianglesForRendering(nullptr), m_verticesForRendering(nullptr) { }
        ~TriangleMeshNode();
        
        uint64_t addVertex(const SLR::Vertex &v);
        void addTriangles(const SurfaceMaterialRef &mat, const Normal3DTextureRef &normalMap, const FloatTextureRef &alphaMap, const std::vector<Triangle> &&triangles);
        
        NodeRef copy() const override;
        
        void applyTransform(const SLR::StaticTransform &t) final;
        
        void applyTransformForRendering(const SLR::StaticTransform &tf) override;
        void createSurfaceObjects() final;
    };
}

#endif /* defined(__SLRSceneGraph__TriangleMeshNode__) */

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
    
    struct Triangle {
        uint64_t vIdx[3];
    };
    
    class TriangleMeshNode : public SurfaceObjectNode {
        SLR::Triangle* m_trianglesForRendering;
        
        std::vector<Vertex> m_vertices;
        std::vector<Triangle> m_triangles;
        std::vector<SurfaceMaterialRef> m_materials;
        std::vector<Normal3DTextureRef> m_normalMaps;
    public:
        TriangleMeshNode() : SurfaceObjectNode(), m_trianglesForRendering(nullptr) { };
        ~TriangleMeshNode();
        
        uint64_t addVertex(const SLR::Vertex &v);
        void addTriangle(uint64_t v0, uint64_t v1, uint64_t v2, const SurfaceMaterialRef &mat, const Normal3DTextureRef &normalMap = nullptr);
        
        void applyTransform(const SLR::StaticTransform &t) final;
        void createSurfaceObjects() final;
    };
}

#endif /* defined(__SLRSceneGraph__TriangleMeshNode__) */

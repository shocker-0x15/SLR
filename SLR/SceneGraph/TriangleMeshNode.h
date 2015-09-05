//
//  TriangleMeshNode.h
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__TriangleMeshNode__
#define __SLR__TriangleMeshNode__

#include "../defines.h"
#include "../references.h"
#include "nodes.h"
#include "../Core/geometry.h"
#include "../Surface/TriangleMesh.h"

class TriangleMeshNode : public SurfaceObjectNode {
    friend class Triangle;
    std::vector<Vertex> m_verticesOrg;
    std::vector<Vertex> m_vertices;
    std::vector<Triangle> m_triangles;
    std::vector<SurfaceMaterialRef> m_materials;
    std::vector<Normal3DTextureRef> m_normalMaps;
public:
    TriangleMeshNode() : SurfaceObjectNode() { };
    
    uint64_t addVertex(const Vertex &v);
    void addTriangle(uint64_t v0, uint64_t v1, uint64_t v2, const SurfaceMaterialRef &mat, const Normal3DTextureRef &normalMap = nullptr);
    
    void applyTransform(const StaticTransform &t) final;
    void createSurfaceObjects() final;
};

#endif /* defined(__SLR__TriangleMeshNode__) */

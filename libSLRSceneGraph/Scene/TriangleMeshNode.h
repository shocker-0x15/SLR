//
//  TriangleMeshNode.h
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLRSceneGraph_TriangleMeshNode__
#define __SLRSceneGraph_TriangleMeshNode__

#include <libSLR/defines.h>
#include <libSLR/Core/geometry.h>
#include "../declarations.h"
#include "node.h"

namespace SLRSceneGraph {
    // There is a possibility to define a unique vertex format suitable for editing.
    // Current implementation uses the vertex for rendering as is.
    typedef SLR::Vertex Vertex;
    
    struct SLR_SCENEGRAPH_API Triangle {
        uint64_t vIdx[3];
        Triangle() {}
        Triangle(uint64_t v0, uint64_t v1, uint64_t v2) : vIdx{v0, v1, v2} {}
    };
    
    class SLR_SCENEGRAPH_API TriangleMeshNode : public SurfaceNode {
    public:
        struct MaterialGroup {
            SurfaceMaterialRef material;
            NormalTextureRef normalMap;
            FloatTextureRef alphaMap;
            std::vector<Triangle> triangles;
        };
    private:
        std::vector<Vertex> m_vertices;
        std::vector<MaterialGroup> m_matGroups;
        bool m_onlyForBoundary;
        int8_t m_axisForRadialTangent; // -1: don't use radial tangent, 0:X, 1:Y, 2:Z
        
        void allocateRawData() override;
        void setupRawData() override;
        void terminateRawData() override;
    public:
        TriangleMeshNode();
        
        uint64_t addVertex(const SLR::Vertex &v);
        void addMaterialGroup(const SurfaceMaterialRef &mat, const NormalTextureRef &normalMap, const FloatTextureRef &alphaMap, 
                              const std::vector<Triangle> &&triangles);
        void useOnlyForBoundary(bool b) {
            m_onlyForBoundary = b;
        }
        void setAxisForRadialTangent(int8_t axisForRadialTangent) {
            m_axisForRadialTangent = axisForRadialTangent;
        }
        
        NodeRef copy() const override;
        
        void applyTransform(const SLR::StaticTransform &t) override;
        
        void prepareForRendering() override;
    };
}

#endif /* __SLRSceneGraph_TriangleMeshNode__ */

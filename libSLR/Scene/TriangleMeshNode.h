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
#include "../SurfaceShape/TriangleSurfaceShape.h"
#include "node.h"

namespace SLR {
    class SLR_API TriangleMeshNode : public SurfaceNode {
    public:
        struct SLR_API MaterialGroup {
            const SurfaceMaterial* m_material;
            const NormalTexture* m_normalMap;
            const FloatTexture* m_alphaMap;
            std::unique_ptr<TriangleSurfaceShape[]> m_triangles;
            uint32_t m_numTriangles;

            MaterialGroup() :
            m_material(nullptr), m_normalMap(nullptr), m_alphaMap(nullptr),
            m_triangles(nullptr), m_numTriangles(0) {
            }
            MaterialGroup(const SurfaceMaterial* material, const NormalTexture* normalMap, const FloatTexture* alphaMap) :
            m_material(material), m_normalMap(normalMap), m_alphaMap(alphaMap),
            m_triangles(nullptr), m_numTriangles(0) {
            }
            ~MaterialGroup() {
            }
            
            void setTriangles(std::unique_ptr<TriangleSurfaceShape[]> &triangles, uint32_t numTriangles);
        };
        
    private:
        std::unique_ptr<Vertex[]> m_vertices;
        uint32_t m_numVertices;
        std::unique_ptr<MaterialGroup[]> m_matGroups;
        uint32_t m_numMatGroups;
        bool m_onlyForBoundary;
        
        std::vector<SurfaceObject*> m_objs;
    public:
        TriangleMeshNode(bool onlyForBoundary) : m_vertices(nullptr), m_onlyForBoundary(onlyForBoundary) { }
        ~TriangleMeshNode() { }
        
        void setVertices(std::unique_ptr<Vertex[]> &vertices, uint32_t numVertices);
        void setMaterialGroups(std::unique_ptr<MaterialGroup[]> &matGroups, uint32_t numMatGroups);
        
        bool isDirectlyTransformable() const override { return true; }
        void createRenderingData(Allocator* mem, const Transform* subTF, RenderingData *data) override;
        void destroyRenderingData(Allocator* mem) override;
    };
}

#endif /* __SLR_TriangleMeshNode__ */

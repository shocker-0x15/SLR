//
//  TriangleMesh.h
//
//  Created by 渡部 心 on 2015/03/22.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_TriangleMesh__
#define __SLR_TriangleMesh__

#include "../defines.h"
#include "../references.h"
#include "../Core/geometry.h"

namespace SLR {
    struct SLR_API Vertex {
        Point3D position;
        Normal3D normal;
        Tangent3D tangent;
        TexCoord2D texCoord;
        
        Vertex() { }
        Vertex(const Point3D &pos, const Normal3D &norm, const Tangent3D &tang, const TexCoord2D &tc) : position(pos), normal(norm), tangent(tang), texCoord(tc) { }
    };
    
    class SLR_API TriangleSurfaceShape : public SurfaceShape {
        const Vertex* m_v[3];
        const FloatTexture* m_alphaTex;
    public:
        TriangleSurfaceShape() : m_v{nullptr, nullptr, nullptr}, m_alphaTex(nullptr) { }
        TriangleSurfaceShape(const Vertex* v0, const Vertex* v1, const Vertex* v2, const FloatTexture* aTex) :
        m_v{v0, v1, v2}, m_alphaTex(aTex) { }
        
        // TODO: consider a better cost value.
        float costForIntersect() const override { return 1.0f; }
        BoundingBox3D bounds() const override;
        BoundingBox3D choppedBounds(BoundingBox3D::Axis chopAxis, float minChopPos, float maxChopPos) const override;
        void splitBounds(BoundingBox3D::Axis splitAxis, float splitPos, BoundingBox3D* bbox0, BoundingBox3D* bbox1) const override;
        bool preTransformed() const override;
        bool intersect(const Ray &ray, SurfaceInteraction* si) const override;
        void getSurfacePoint(const SurfaceInteraction &si, SurfacePoint* surfPt) const override;
        float area() const override;
        void sample(float u0, float u1, SurfacePoint* surfPt, float* areaPDF, DirectionType* posType) const override;
        float evaluateAreaPDF(const SurfacePoint& surfPt) const override;
    };
}

#endif /* __SLR_TriangleMesh__ */

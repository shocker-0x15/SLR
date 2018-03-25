//
//  TriangleSurfaceShape.h
//
//  Created by 渡部 心 on 2015/03/22.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_TriangleSurfaceShape__
#define __SLR_TriangleSurfaceShape__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/geometry.h"

namespace SLR {    
    class SLR_API TriangleSurfaceShape : public SurfaceShape {
        const MaterialGroupInTriangleMesh* m_matGroup;
        uint32_t m_index;
        Vector3D m_texCoord0Dir;
    public:
        TriangleSurfaceShape() : m_matGroup(nullptr), m_index(UINT32_MAX) { }
        TriangleSurfaceShape(const MaterialGroupInTriangleMesh* matGroup, uint32_t index);
        
        // TODO: consider a better cost value.
        float costForIntersect() const override { return 1.0f; }
        BoundingBox3D bounds() const override;
        BoundingBox3D choppedBounds(BoundingBox3D::Axis chopAxis, float minChopPos, float maxChopPos) const override;
        void splitBounds(BoundingBox3D::Axis splitAxis, float splitPos, BoundingBox3D* bbox0, BoundingBox3D* bbox1) const override;
        bool preTransformed() const override;
        bool intersectWithoutAlpha(const Ray &ray, const RaySegment &segment, SurfaceInteraction* si) const;
        bool intersect(const Ray &ray, const RaySegment &segment, LightPathSampler &pathSampler, SurfaceInteraction* si) const override;
        float testVisibility(const Ray &ray, const RaySegment &segment) const override;
        void calculateSurfacePoint(const SurfaceInteraction &si, SurfacePoint* surfPt) const override;
        float area() const override;
        void sample(float u0, float u1, SurfacePoint* surfPt, float* areaPDF, DirectionType* posType) const override;
        float evaluateAreaPDF(const SurfacePoint& surfPt) const override;
    };
}

#endif /* __SLR_TriangleSurfaceShape__ */

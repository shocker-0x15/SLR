//
//  TriangleMesh.cpp
//
//  Created by 渡部 心 on 2015/03/22.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "TriangleMesh.h"
#include "../BasicTypes/CompensatedSum.h"
#include "../Core/distributions.h"
#include "../Core/SurfaceObject.h"
#include "../Core/textures.h"

namespace SLR {
    BoundingBox3D TriangleSurfaceShape::bounds() const {
        return BoundingBox3D(m_v[0]->position).unify(m_v[1]->position).unify(m_v[2]->position);
    }
    
    BoundingBox3D TriangleSurfaceShape::choppedBounds(BoundingBox3D::Axis chopAxis, float minChopPos, float maxChopPos) const {
        const float chopPos[2] = {minChopPos, maxChopPos};
        
        Point3D p[3] = {m_v[0]->position, m_v[1]->position, m_v[2]->position};
        std::sort(p, p + 3, [chopAxis](const Point3D &pa, const Point3D &pb) { return pa[chopAxis] < pb[chopAxis]; });
        float minPos = p[0][chopAxis];
        float maxPos = p[2][chopAxis];
        
        // early exits for the following 3 patterns:
        // -- |    |
        //    |    | --
        //    | -- |
        BoundingBox3D baseBBox = bounds();
        if (minPos >= maxChopPos || maxPos <= minChopPos)
            return BoundingBox3D();
        if (minPos >= minChopPos && maxPos <= maxChopPos)
            return baseBBox;
        
        // A point P is in left of an ordinate O: P < O
        //                 right                : P >= O
        // remaining 3 patterns:
        //  --|----|--
        //  --|--  |
        //    |  --|--
        uint32_t idxIsect = 0;
        Point3D isectPs[4];
        for (int from = 0; from < 3 - 1; ++from) {
            const Point3D &pa = p[from];
            for (int to = from + 1; to < 3; ++to) {
                const Point3D &pb = p[to];
                float deltaAB = pb[chopAxis] - pa[chopAxis];
                for (int cp = 0; cp < 2; ++cp) {
                    float deltaAP = chopPos[cp] - pa[chopAxis];
                    float deltaPB = chopPos[cp] - pb[chopAxis];
                    float t = deltaAP / deltaAB;
                    if (deltaAP > 0 && deltaPB <= 0) // pa[chopAxis] < chopPos[cp] && pb[chopAxis] >= chopPos[cp]
                        isectPs[idxIsect++] = (1 - t) * pa + t * pb;// pa + (pb - pa) * t = (1 - t) * pa + pb // This form is preferrable for precision?
                }
            }
        }
        SLRAssert(idxIsect == 2 || idxIsect == 4, "The number of intersection should be 2 or 4.");
        
        BoundingBox3D ret;
        if (p[1][chopAxis] >= minChopPos && p[1][chopAxis] < maxChopPos)
            ret.unify(p[1]);
        for (int i = 0; i < idxIsect; ++i)
            ret.unify(isectPs[i]);
        if (idxIsect == 2)
            ret.unify(maxPos < maxChopPos ? p[2] : p[0]);
        
        SLRAssert(realGE(ret.minP[chopAxis], minChopPos, 1e-6f) && realLE(ret.maxP[chopAxis], maxChopPos, 1e-6f) &&
                  realLE(ret.surfaceArea(), baseBBox.surfaceArea(), 1e-6f), "invalid chopped bounds.");
        
        return ret;
    }
    
    void TriangleSurfaceShape::splitBounds(BoundingBox3D::Axis splitAxis, float splitPos, BoundingBox3D* bbox0, BoundingBox3D* bbox1) const {
        Point3D p[3] = {m_v[0]->position, m_v[1]->position, m_v[2]->position};
        std::sort(p, p + 3, [splitAxis](const Point3D &pa, const Point3D &pb) { return pa[splitAxis] < pb[splitAxis]; });
        float minPos = p[0][splitAxis];
        float maxPos = p[2][splitAxis];
        
        // early exits for the following 2 patterns:
        // -- |
        //    | --
        BoundingBox3D baseBBox = bounds();
        if (splitPos <= minPos) {
            *bbox0 = BoundingBox3D();
            *bbox1 = baseBBox;
            return;
        }
        if (splitPos >= maxPos) {
            *bbox0 = baseBBox;
            *bbox1 = BoundingBox3D();
            return;
        }
        
        // A point P is in left of an ordinate O: P < O
        //                 right                : P >= O
        uint32_t idxIsect = 0;
        Point3D isectPs[2];
        for (int from = 0; from < 3 - 1; ++from) {
            const Point3D &pa = p[from];
            for (int to = from + 1; to < 3; ++to) {
                const Point3D &pb = p[to];
                float deltaAB = pb[splitAxis] - pa[splitAxis];
                float deltaAP = splitPos - pa[splitAxis];
                float deltaPB = splitPos - pb[splitAxis];
                float t = deltaAP / deltaAB;
                if (deltaAP > 0 && deltaPB <= 0) // pa[chopAxis] < chopPos[cp] && pb[chopAxis] >= chopPos[cp]
                    isectPs[idxIsect++] = (1 - t) * pa + t * pb;// pa + (pb - pa) * t = (1 - t) * pa + pb // This form is preferrable for precision?
            }
        }
        SLRAssert(idxIsect == 2, "The number of intersection should be 2.");
        
        *bbox0 = BoundingBox3D(p[0]);
        *bbox1 = BoundingBox3D(p[2]);
        if (p[1][splitAxis] < splitPos)
            bbox0->unify(p[1]);
        else
            bbox1->unify(p[1]);
        for (int i = 0; i < idxIsect; ++i) {
            bbox0->unify(isectPs[i]);
            bbox1->unify(isectPs[i]);
        }
        SLRAssert(realLE(bbox0->maxP[splitAxis], splitPos, 1e-6f) && realGE(bbox1->minP[splitAxis], splitPos, 1e-6f), "invalid split bounds.");
    }
    
    bool TriangleSurfaceShape::preTransformed() const {
        return true;
    }
    
    bool TriangleSurfaceShape::intersect(const Ray &ray, SurfaceInteraction* si) const {
        const Vertex &v0 = *m_v[0];
        const Vertex &v1 = *m_v[1];
        const Vertex &v2 = *m_v[2];
        
        Vector3D edge01 = v1.position - v0.position;
        Vector3D edge02 = v2.position - v0.position;
        
        Vector3D p = cross(ray.dir, edge02);
        float det = dot(edge01, p);
        if (det == 0.0f)
            return false;
        float invDet = 1.0f / det;
        
        Vector3D d = ray.org - v0.position;
        
        float b1 = dot(d, p) * invDet;
        if (b1 < 0.0f || b1 > 1.0f)
            return false;
        
        Vector3D q = cross(d, edge01);
        
        float b2 = dot(ray.dir, q) * invDet;
        if (b2 < 0.0f || b1 + b2 > 1.0f)
            return false;
        
        float tt = dot(edge02, q) * invDet;
        if (tt < ray.distMin || tt > ray.distMax)
            return false;
        
        float b0 = 1.0f - b1 - b2;
        TexCoord2D texCoord = b0 * v0.texCoord + b1 * v1.texCoord + b2 * v2.texCoord;
        // Check if an alpha value at the intersection point is zero or not.
        // If zero, intersection doesn't occur.
        if (m_alphaTex) {
            if (m_alphaTex->evaluate(texCoord) == 0.0f)
                return false;
        }
        
        *si = SurfaceInteraction(ray.time, // ------------------------- time
                                 tt, // ------------------------------- distance
                                 ray.org + ray.dir * tt, // ----------- position in world coordinate
                                 normalize(cross(edge01, edge02)), // - geometric normal in world coordinate
                                 b0, b1, // --------------------------- surface parameters
                                 texCoord //--------------------------- texture coordinate
                                 );
        
        return true;
    }
    
    void TriangleSurfaceShape::getSurfacePoint(const SurfaceInteraction &si, SurfacePoint* surfPt) const {
        const Vertex &v0 = *m_v[0];
        const Vertex &v1 = *m_v[1];
        const Vertex &v2 = *m_v[2];
        float b0, b1;
        si.getSurfaceParameter(&b0, &b1);
        float b2 = 1.0f - b0 - b1;
        
        // TODO: It might be preferable to explicitly calculate this direction from texture coordinates at the setup time and add an attribute to Triangle.
        Vector3D texCoord0Dir;
        Vector3D dP0 = v0.position - v2.position;
        Vector3D dP1 = v1.position - v2.position;
        TexCoord2D dTC0 = v0.texCoord - v2.texCoord;
        TexCoord2D dTC1 = v1.texCoord - v2.texCoord;
        float detTC = dTC0.u * dTC1.v - dTC0.v * dTC1.u;
        if (detTC != 0)
            texCoord0Dir = normalize((1.0f / detTC) * Vector3D(dTC1.v * dP0.x - dTC0.v * dP1.x,
                                                               dTC1.v * dP0.y - dTC0.v * dP1.y,
                                                               dTC1.v * dP0.z - dTC0.v * dP1.z));
        else
            texCoord0Dir = normalize(dP0);
        
        ReferenceFrame shadingFrame;
        shadingFrame.z = normalize(b0 * v0.normal + b1 * v1.normal + b2 * v2.normal);
        shadingFrame.x = normalize(b0 * v0.tangent + b1 * v1.tangent + b2 * v2.tangent);
        // guarantee the orthogonality between the normal and tangent.
        // Orthogonality break might be caused by barycentric interpolation?
        float dotNT = dot(shadingFrame.z, shadingFrame.x);
        if (std::fabs(dotNT) >= 0.01f)
            shadingFrame.x = normalize(shadingFrame.x - dotNT * shadingFrame.z);
        shadingFrame.y = cross(shadingFrame.z, shadingFrame.x);
        
        *surfPt = SurfacePoint(si, false, shadingFrame, texCoord0Dir);
    }
    
    float TriangleSurfaceShape::area() const {
        const Point3D &p0 = m_v[0]->position;
        const Point3D &p1 = m_v[1]->position;
        const Point3D &p2 = m_v[2]->position;
        return 0.5f * cross(p1 - p0, p2 - p0).length();
    }
    
    void TriangleSurfaceShape::sample(float u0, float u1, SurfacePoint* surfPt, float* areaPDF, DirectionType* posType) const {
        //    const SurfaceMaterial* mat = m_mat;
        float b0, b1, b2;
        uniformSampleTriangle(u0, u1, &b0, &b1);
        b2 = 1.0f - b0 - b1;
        
        const Vertex v0 = *m_v[0];
        const Vertex v1 = *m_v[1];
        const Vertex v2 = *m_v[2];

        Vector3D dP0 = v0.position - v2.position;
        Vector3D dP1 = v1.position - v2.position;
        TexCoord2D dTC0 = v0.texCoord - v2.texCoord;
        TexCoord2D dTC1 = v1.texCoord - v2.texCoord;
        float invDetTC = 1.0f / (dTC0.u * dTC1.v - dTC0.v * dTC1.u);
        Vector3D texCoord0Dir = invDetTC * Vector3D(dTC1.v * dP0.x - dTC0.v * dP1.x,
                                                    dTC1.v * dP0.y - dTC0.v * dP1.y,
                                                    dTC1.v * dP0.z - dTC0.v * dP1.z);
        
        ReferenceFrame shadingFrame;
        shadingFrame.z = normalize(b0 * v0.normal + b1 * v1.normal + b2 * v2.normal);
        shadingFrame.x = normalize(b0 * v0.tangent + b1 * v1.tangent + b2 * v2.tangent);
        shadingFrame.y = cross(shadingFrame.z, shadingFrame.x);
        
        *surfPt = SurfacePoint(b0 * v0.position + b1 * v1.position + b2 * v2.position,
                               false,
                               shadingFrame,
                               cross(v1.position - v0.position, v2.position - v0.position).normalize(),
                               b0, b1,
                               b0 * v0.texCoord + b1 * v1.texCoord + b2 * v2.texCoord,
                               texCoord0Dir
                               );
        *areaPDF = 1.0f / area();
        *posType = DirectionType::LowFreq;
    }
    
    float TriangleSurfaceShape::evaluateAreaPDF(const SurfacePoint& surfPt) const {
        float u, v;
        surfPt.getSurfaceParameter(&u, &v);
        SLRAssert(u + v <= 1.0f, "Invalid parameters for a triangle.");
        return 1.0f / area();
    }
}

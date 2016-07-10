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
    BoundingBox3D Triangle::bounds() const {
        return BoundingBox3D(m_v[0]->position).unify(m_v[1]->position).unify(m_v[2]->position);
    }
    
    // TODO: improve the logic.
    BoundingBox3D Triangle::choppedBounds(BoundingBox3D::Axis chopAxis, float minChopPos, float maxChopPos) const {
        const Point3D &p0 = m_v[0]->position;
        const Point3D &p1 = m_v[1]->position;
        const Point3D &p2 = m_v[2]->position;
        float minPos = std::min(p0[chopAxis], std::min(p1[chopAxis], p2[chopAxis]));
        float maxPos = std::max(p0[chopAxis], std::max(p1[chopAxis], p2[chopAxis]));
        
        BoundingBox3D baseBBox = bounds();
        if (maxChopPos < minPos)
            return BoundingBox3D();
        if (minChopPos > maxPos)
            return BoundingBox3D();
        if (minChopPos < minPos && maxChopPos > maxPos)
            return baseBBox;
        
        if (minChopPos >= minPos && maxChopPos <= maxPos) {
            float p0P0 = minChopPos - p0[chopAxis], p1P0 = minChopPos - p1[chopAxis], p2P0 = minChopPos - p2[chopAxis];
            float p0P1 = maxChopPos - p0[chopAxis], p1P1 = maxChopPos - p1[chopAxis], p2P1 = maxChopPos - p2[chopAxis];
            
            bool isectP0[3] = {p0P0 * p1P0 <= 0, p0P0 * p2P0 <= 0, p1P0 * p2P0 <= 0}; // 01, 02, 12
            bool isectP1[3] = {p0P1 * p1P1 <= 0, p0P1 * p2P1 <= 0, p1P1 * p2P1 <= 0}; // 01, 02, 12

            Point3D iP0[3];
            Point3D iP1[3];
            if (isectP0[0])
                iP0[0] = p0 + (p1 - p0) * (p0P0 / (p1[chopAxis] - p0[chopAxis]));
            if (isectP0[1])
                iP0[1] = p0 + (p2 - p0) * (p0P0 / (p2[chopAxis] - p0[chopAxis]));
            if (isectP0[2])
                iP0[2] = p1 + (p2 - p1) * (p1P0 / (p2[chopAxis] - p1[chopAxis]));
            
            if (isectP1[0])
                iP1[0] = p0 + (p1 - p0) * (p0P1 / (p1[chopAxis] - p0[chopAxis]));
            if (isectP1[1])
                iP1[1] = p0 + (p2 - p0) * (p0P1 / (p2[chopAxis] - p0[chopAxis]));
            if (isectP1[2])
                iP1[2] = p1 + (p2 - p1) * (p1P1 / (p2[chopAxis] - p1[chopAxis]));
            
            BoundingBox3D ret;
            for (int i = 0; i < 3; ++i) {
                if (isectP0[i])
                    ret.unify(iP0[i]);
                if (isectP1[i])
                    ret.unify(iP1[i]);
            }
            
            if (p0[chopAxis] >= minChopPos && p0[chopAxis] <= maxChopPos)
                ret.unify(p0);
            else if (p1[chopAxis] >= minChopPos && p1[chopAxis] <= maxChopPos)
                ret.unify(p1);
            else if (p2[chopAxis] >= minChopPos && p2[chopAxis] <= maxChopPos)
                ret.unify(p2);
            SLRAssert(!ret.hasNaN() && !ret.hasInf(), "Invalid chopped bounding box.");
            
            return ret;
        }
        else {
            float splitPos = maxPos < maxChopPos ? minChopPos : maxChopPos;
            
            float p0PlaneDist = splitPos - p0[chopAxis];
            float p1PlaneDist = splitPos - p1[chopAxis];
            float p2PlaneDist = splitPos - p2[chopAxis];
            BoundingBox3D bboxA, bboxB;
            if (p0PlaneDist * p1PlaneDist < 0 && p0PlaneDist * p2PlaneDist < 0) { // p0 is alone
                Point3D p01 = p0 + (p1 - p0) * (p0PlaneDist / (p1[chopAxis] - p0[chopAxis]));
                Point3D p02 = p0 + (p2 - p0) * (p0PlaneDist / (p2[chopAxis] - p0[chopAxis]));
                bboxA.unify(p0).unify(p01).unify(p02);
                bboxB.unify(p1).unify(p2).unify(p01).unify(p02);
            }
            else if (p0PlaneDist * p1PlaneDist >= 0) { // p2 is alone
                Point3D p02 = p0 + (p2 - p0) * (p0PlaneDist / (p2[chopAxis] - p0[chopAxis]));
                Point3D p12 = p1 + (p2 - p1) * (p1PlaneDist / (p2[chopAxis] - p1[chopAxis]));
                bboxA.unify(p0).unify(p1).unify(p02).unify(p12);
                bboxB.unify(p2).unify(p02).unify(p12);
            }
            else { // p1 is alone
                Point3D p01 = p0 + (p1 - p0) * (p0PlaneDist / (p1[chopAxis] - p0[chopAxis]));
                Point3D p21 = p2 + (p1 - p2) * (p2PlaneDist / (p1[chopAxis] - p2[chopAxis]));
                bboxA.unify(p0).unify(p2).unify(p01).unify(p21);
                bboxB.unify(p1).unify(p01).unify(p21);
            }
            
            if (maxPos < maxChopPos)
                return p0PlaneDist > 0 ? bboxB : bboxA;
            else
                return p0PlaneDist > 0 ? bboxA : bboxB;
        }
    }
    
    void Triangle::splitBounds(BoundingBox3D::Axis splitAxis, float splitPos, BoundingBox3D* bbox0, BoundingBox3D* bbox1) const {
        const Point3D &p0 = m_v[0]->position;
        const Point3D &p1 = m_v[1]->position;
        const Point3D &p2 = m_v[2]->position;
        float minPos = std::min(p0[splitAxis], std::min(p1[splitAxis], p2[splitAxis]));
        float maxPos = std::max(p0[splitAxis], std::max(p1[splitAxis], p2[splitAxis]));
        
        BoundingBox3D baseBBox = bounds();
        if (splitPos < minPos) {
            *bbox0 = BoundingBox3D();
            *bbox1 = baseBBox;
            return;
        }
        if (splitPos > maxPos) {
            *bbox0 = baseBBox;
            *bbox1 = BoundingBox3D();
            return;
        }
        
        float p0PlaneDist = splitPos - p0[splitAxis];
        float p1PlaneDist = splitPos - p1[splitAxis];
        float p2PlaneDist = splitPos - p2[splitAxis];
        BoundingBox3D bboxA, bboxB;
        if (p0PlaneDist * p1PlaneDist < 0 && p0PlaneDist * p2PlaneDist < 0) {
            Point3D p01 = p0 + (p1 - p0) * (p0PlaneDist / (p1[splitAxis] - p0[splitAxis]));
            Point3D p02 = p0 + (p2 - p0) * (p0PlaneDist / (p2[splitAxis] - p0[splitAxis]));
            bboxA.unify(p0).unify(p01).unify(p02);
            bboxB.unify(p1).unify(p2).unify(p01).unify(p02);
        }
        else if (p0PlaneDist * p1PlaneDist >= 0) {
            Point3D p02 = p0 + (p2 - p0) * (p0PlaneDist / (p2[splitAxis] - p0[splitAxis]));
            Point3D p12 = p1 + (p2 - p1) * (p1PlaneDist / (p2[splitAxis] - p1[splitAxis]));
            bboxA.unify(p0).unify(p1).unify(p02).unify(p12);
            bboxB.unify(p2).unify(p02).unify(p12);
        }
        else {
            Point3D p01 = p0 + (p1 - p0) * (p0PlaneDist / (p1[splitAxis] - p0[splitAxis]));
            Point3D p21 = p2 + (p1 - p2) * (p2PlaneDist / (p1[splitAxis] - p2[splitAxis]));
            bboxA.unify(p0).unify(p2).unify(p01).unify(p21);
            bboxB.unify(p1).unify(p01).unify(p21);
        }
        
        *bbox0 = p0PlaneDist > 0 ? bboxA : bboxB;
        *bbox1 = p0PlaneDist > 0 ? bboxB : bboxA;
    }
    
    bool Triangle::preTransformed() const {
        return true;
    }
    
    bool Triangle::intersect(const Ray &ray, Intersection* isect) const {
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
        
        isect->dist = tt;
        isect->p = ray.org + ray.dir * tt;
        isect->gNormal = normalize(cross(edge01, edge02));
        isect->u = b0;
        isect->v = b1;
        isect->texCoord = texCoord;
        
        return true;
    }
    
    void Triangle::getSurfacePoint(const Intersection &isect, SurfacePoint* surfPt) const {
        surfPt->p = isect.p;
        surfPt->atInfinity = false;
        surfPt->gNormal = isect.gNormal;
        surfPt->u = isect.u;
        surfPt->v = isect.v;
        surfPt->texCoord = isect.texCoord;
        
        const Vertex &v0 = *m_v[0];
        const Vertex &v1 = *m_v[1];
        const Vertex &v2 = *m_v[2];
        float b0 = isect.u, b1 = isect.v;
        float b2 = 1.0f - b0 - b1;
        
        // TODO: It might be preferable to explicitly calculate this direction from texture coordinates at the setup time and add an attribute to Triangle.
        Vector3D dP0 = v0.position - v2.position;
        Vector3D dP1 = v1.position - v2.position;
        TexCoord2D dTC0 = v0.texCoord - v2.texCoord;
        TexCoord2D dTC1 = v1.texCoord - v2.texCoord;
        float detTC = dTC0.u * dTC1.v - dTC0.v * dTC1.u;
        if (detTC != 0)
            surfPt->texCoord0Dir = normalize((1.0f / detTC) * Vector3D(dTC1.v * dP0.x - dTC0.v * dP1.x,
                                                                       dTC1.v * dP0.y - dTC0.v * dP1.y,
                                                                       dTC1.v * dP0.z - dTC0.v * dP1.z));
        else
            surfPt->texCoord0Dir = normalize(dP0);
        
        surfPt->shadingFrame.z = normalize(b0 * v0.normal + b1 * v1.normal + b2 * v2.normal);
        surfPt->shadingFrame.x = normalize(b0 * v0.tangent + b1 * v1.tangent + b2 * v2.tangent);
        // guarantee the orthogonality between the normal and tangent.
        // Orthogonality break might be caused by barycentric interpolation?
        float dotNT = dot(surfPt->shadingFrame.z, surfPt->shadingFrame.x);
        if (std::fabs(dotNT) >= 0.01f)
            surfPt->shadingFrame.x = SLR::normalize(surfPt->shadingFrame.x - dotNT * surfPt->shadingFrame.z);
        surfPt->shadingFrame.y = cross(surfPt->shadingFrame.z, surfPt->shadingFrame.x);
    }
    
    float Triangle::area() const {
        const Point3D &p0 = m_v[0]->position;
        const Point3D &p1 = m_v[1]->position;
        const Point3D &p2 = m_v[2]->position;
        return 0.5f * cross(p1 - p0, p2 - p0).length();
    }
    
    void Triangle::sample(float u0, float u1, SurfacePoint* surfPt, float* areaPDF) const {
        //    const SurfaceMaterial* mat = m_mat;
        float b0, b1, b2;
        uniformSampleTriangle(u0, u1, &b0, &b1);
        b2 = 1.0f - b0 - b1;
        
        const Vertex v0 = *m_v[0];
        const Vertex v1 = *m_v[1];
        const Vertex v2 = *m_v[2];
        
        surfPt->p = b0 * v0.position + b1 * v1.position + b2 * v2.position;
        surfPt->atInfinity = false;
        surfPt->gNormal = cross(v1.position - v0.position, v2.position - v0.position).normalize();
        surfPt->u = b0;
        surfPt->v = b1;
        surfPt->texCoord = b0 * v0.texCoord + b1 * v1.texCoord + b2 * v2.texCoord;
        
        Vector3D dP0 = v0.position - v2.position;
        Vector3D dP1 = v1.position - v2.position;
        TexCoord2D dTC0 = v0.texCoord - v2.texCoord;
        TexCoord2D dTC1 = v1.texCoord - v2.texCoord;
        float invDetTC = 1.0f / (dTC0.u * dTC1.v - dTC0.v * dTC1.u);
        surfPt->texCoord0Dir = invDetTC * Vector3D(dTC1.v * dP0.x - dTC0.v * dP1.x,
                                                   dTC1.v * dP0.y - dTC0.v * dP1.y,
                                                   dTC1.v * dP0.z - dTC0.v * dP1.z);
        
        surfPt->shadingFrame.z = normalize(b0 * v0.normal + b1 * v1.normal + b2 * v2.normal);
        surfPt->shadingFrame.x = normalize(b0 * v0.tangent + b1 * v1.tangent + b2 * v2.tangent);
        surfPt->shadingFrame.y = cross(surfPt->shadingFrame.z, surfPt->shadingFrame.x);
        
        *areaPDF = 1.0f / area();
    }
    
    float Triangle::evaluateAreaPDF(const SurfacePoint& surfPt) const {
        SLRAssert(surfPt.u + surfPt.v <= 1.0f, "Invalid parameters for a triangle.");
        return 1.0f / area();
    }    
}

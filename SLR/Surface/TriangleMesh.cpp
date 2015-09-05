//
//  TriangleMesh.cpp
//
//  Created by 渡部 心 on 2015/03/22.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "TriangleMesh.h"
#include "TriangleMeshNode.h"
#include "CompensatedSum.h"
#include "distributions.h"
#include "SurfaceObject.h"

BoundingBox3D Triangle::bounds() const {
    return BoundingBox3D(m_mesh->m_vertices[m_v[0]].position).unify(m_mesh->m_vertices[m_v[1]].position).unify(m_mesh->m_vertices[m_v[2]].position);
}

bool Triangle::preTransformed() const {
    return true;
}

bool Triangle::intersect(const Ray &ray, Intersection* isect) const {
    const Vertex &v0 = m_mesh->m_vertices[m_v[0]];
    const Vertex &v1 = m_mesh->m_vertices[m_v[1]];
    const Vertex &v2 = m_mesh->m_vertices[m_v[2]];
    
    Vector3D edge01 = v1.position - v0.position;
    Vector3D edge02 = v2.position - v0.position;
    
    Vector3D p = cross(ray.dir, edge02);
    float det = dot(edge01, p);
    if (det > -0.000001f && det < 0.000001f)
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
    isect->dist = tt;
    isect->p = ray.org + ray.dir * tt;
    isect->gNormal = normalize(cross(edge01, edge02));
    isect->u = b0;
    isect->v = b1;
    isect->texCoord = b0 * v0.texCoord + b1 * v1.texCoord + b2 * v2.texCoord;
    
//    // Check if the alpha value at the intersection point is zero or not.
//    // If zero, intersection doesn't occur.
//    if (alphaTexp) {
//        if (evaluateAlphaTexture(scene->texturesData + face->alphaTexPtr, isect->uv) == 0.0f)
//            return false;
//    }

    return true;
}

void Triangle::getSurfacePoint(const Intersection &isect, SurfacePoint* surfPt) const {
    surfPt->p = isect.p;
    surfPt->atInfinity = false;
    surfPt->gNormal = isect.gNormal;
    surfPt->u = isect.u;
    surfPt->v = isect.v;
    surfPt->texCoord = isect.texCoord;
    
    const Vertex &v0 = m_mesh->m_vertices[m_v[0]];
    const Vertex &v1 = m_mesh->m_vertices[m_v[1]];
    const Vertex &v2 = m_mesh->m_vertices[m_v[2]];
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
    surfPt->shadingFrame.y = cross(surfPt->shadingFrame.z, surfPt->shadingFrame.x);
}

float Triangle::area() const {
    const Point3D &p0 = m_mesh->m_vertices[m_v[0]].position;
    const Point3D &p1 = m_mesh->m_vertices[m_v[1]].position;
    const Point3D &p2 = m_mesh->m_vertices[m_v[2]].position;
    return 0.5f * cross(p1 - p0, p2 - p0).length();
}

void Triangle::sample(float u0, float u1, SurfacePoint* surfPt, float* areaPDF) const {
//    const SurfaceMaterial* mat = m_mat;
    float b0, b1, b2;
    uniformSampleTriangle(u0, u1, &b0, &b1);
    b2 = 1.0f - b0 - b1;
    
    const Vertex v0 = m_mesh->m_vertices[m_v[0]];
    const Vertex v1 = m_mesh->m_vertices[m_v[1]];
    const Vertex v2 = m_mesh->m_vertices[m_v[2]];
    
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

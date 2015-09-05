//
//  PerspectiveCamera.cpp
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "PerspectiveCamera.h"
#include "../Core/Transform.h"
#include "../Memory/ArenaAllocator.h"
#include "../Core/distributions.h"

Spectrum PerspectiveCamera::sample(const LensPosQuery &query, const LensPosSample &smp, LensPosQueryResult* result) const {
    StaticTransform staticTF;
    if (m_transform)
        m_transform->sample(query.time, &staticTF);
    
    float lx, ly;
    concentricSampleDisk(smp.uPos[0], smp.uPos[1], &lx, &ly);
    Point3D orgLocal = Point3D(m_lensRadius * lx, m_lensRadius * ly, 0.0f);
    SurfacePoint &surfPt = result->surfPt;
    surfPt.p = staticTF * orgLocal;
    surfPt.gNormal = staticTF * Normal3D(0, 0, 1);
    surfPt.u = lx;
    surfPt.v = ly;
    surfPt.texCoord = TexCoord2D::Zero;
    surfPt.texCoord0Dir = Vector3D::Zero;
    surfPt.shadingFrame.z = (Vector3D)surfPt.gNormal;
    surfPt.shadingFrame.x = staticTF * Vector3D(1, 0, 0);// assume the transform doesn't include scaling.
    surfPt.shadingFrame.y = cross(surfPt.shadingFrame.z, surfPt.shadingFrame.x);
    surfPt.obj = nullptr;
    result->isDeltaPos = m_lensRadius == 0.0f;
    result->areaPDF = m_lensRadius > 0.0f ? 1.0f / (M_PI * m_lensRadius * m_lensRadius) : 1.0f;
    
    return Spectrum::One;
}

IDF* PerspectiveCamera::createIDF(const SurfacePoint &surfPt, ArenaAllocator &mem) const {
    return mem.create<PerspectiveIDF>(*this, Point3D(m_lensRadius * surfPt.u, m_lensRadius * surfPt.v, 0.0f));
}

Spectrum PerspectiveIDF::sample(const IDFSample &smp, IDFQueryResult *result) const {
    Point3D pFocus = Point3D(m_cam.m_opWidth * (0.5f - smp.uDir[0]),
                             m_cam.m_opHeight * (0.5f - smp.uDir[1]),
                             m_cam.m_objPlaneDistance);
    
    Vector3D dirLocal = normalize(pFocus - m_orgLocal);
    result->dirLocal = dirLocal;
    result->dirPDF = m_cam.m_imgPlaneDistance * m_cam.m_imgPlaneDistance / ((dirLocal.z * dirLocal.z * dirLocal.z) * m_cam.m_imgPlaneArea);
    
    return Spectrum::One;
}

Spectrum PerspectiveIDF::evaluate(const Vector3D &dirIn) const {
    Point3D pFocas = dirIn * (m_cam.m_objPlaneDistance / dirIn.z) + m_orgLocal;
    bool valid = (pFocas.x >= -m_cam.m_opWidth * 0.5f && pFocas.x <= m_cam.m_opWidth * 0.5f &&
                  pFocas.y >= -m_cam.m_opHeight * 0.5f && pFocas.y <= m_cam.m_opHeight * 0.5f);
    return valid ? Spectrum::One : Spectrum::Zero;
}

float PerspectiveIDF::evaluatePDF(const Vector3D &dirIn) const {
    Point3D pFocas = dirIn * (m_cam.m_objPlaneDistance / dirIn.z) + m_orgLocal;
    bool valid = (pFocas.x >= -m_cam.m_opWidth * 0.5f && pFocas.x <= m_cam.m_opWidth * 0.5f &&
                  pFocas.y >= -m_cam.m_opHeight * 0.5f && pFocas.y <= m_cam.m_opHeight * 0.5f);
    float dirPDF = m_cam.m_imgPlaneDistance * m_cam.m_imgPlaneDistance / ((dirIn.z * dirIn.z * dirIn.z) * m_cam.m_imgPlaneArea);
    return valid ? dirPDF : 0.0f;
}

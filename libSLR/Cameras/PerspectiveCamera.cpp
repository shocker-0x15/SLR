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
#include "../Core/ImageSensor.h"

namespace SLR {
    PerspectiveCamera::PerspectiveCamera(float sensitivity,
                                         float aspect, float fovY, float lensRadius, float imgPDist, float objPDist) :
    Camera(nullptr),
    m_aspect(aspect), m_fovY(fovY), m_lensRadius(lensRadius), m_imgPlaneDistance(imgPDist), m_objPlaneDistance(objPDist) {
        m_opHeight = 2.0f * m_objPlaneDistance * std::tan(m_fovY * 0.5f);
        m_opWidth = m_opHeight * m_aspect;
        m_imgPlaneArea = m_opWidth * m_opHeight * std::pow(m_imgPlaneDistance / m_objPlaneDistance, 2);
        
        m_sensor = createUnique<ImageSensor>(sensitivity > 0 ? sensitivity : (1.0f / (M_PI * lensRadius * lensRadius)));
    }
    
    PerspectiveCamera::~PerspectiveCamera() {
    }
    
    ImageSensor* PerspectiveCamera::getSensor() const {
        return m_sensor.get();
    }
    
    SampledSpectrum PerspectiveCamera::sample(const LensPosQuery &query, const LensPosSample &smp, LensPosQueryResult* result) const {
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
        
        return SampledSpectrum::One;
    }
    
    IDF* PerspectiveCamera::createIDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const {
        return mem.create<PerspectiveIDF>(*this, Point3D(m_lensRadius * surfPt.u, m_lensRadius * surfPt.v, 0.0f));
    }
    
    SampledSpectrum PerspectiveIDF::sample(const IDFSample &smp, IDFQueryResult *result) const {
        Point3D pFocus = Point3D(m_cam.m_opWidth * (0.5f - smp.uDir[0]),
                                 m_cam.m_opHeight * (0.5f - smp.uDir[1]),
                                 m_cam.m_objPlaneDistance);
        
        Vector3D dirLocal = normalize(pFocus - m_orgLocal);
        result->dirLocal = dirLocal;
        result->dirPDF = m_cam.m_imgPlaneDistance * m_cam.m_imgPlaneDistance / ((dirLocal.z * dirLocal.z * dirLocal.z) * m_cam.m_imgPlaneArea);
        
        return SampledSpectrum::One;
    }
    
    SampledSpectrum PerspectiveIDF::evaluate(const Vector3D &dirIn) const {
        Point3D pFocas = dirIn * (m_cam.m_objPlaneDistance / dirIn.z) + m_orgLocal;
        bool valid = (pFocas.x >= -m_cam.m_opWidth * 0.5f && pFocas.x <= m_cam.m_opWidth * 0.5f &&
                      pFocas.y >= -m_cam.m_opHeight * 0.5f && pFocas.y <= m_cam.m_opHeight * 0.5f);
        return valid ? SampledSpectrum::One : SampledSpectrum::Zero;
    }
    
    float PerspectiveIDF::evaluatePDF(const Vector3D &dirIn) const {
        Point3D pFocas = dirIn * (m_cam.m_objPlaneDistance / dirIn.z) + m_orgLocal;
        bool valid = (pFocas.x >= -m_cam.m_opWidth * 0.5f && pFocas.x <= m_cam.m_opWidth * 0.5f &&
                      pFocas.y >= -m_cam.m_opHeight * 0.5f && pFocas.y <= m_cam.m_opHeight * 0.5f);
        float dirPDF = m_cam.m_imgPlaneDistance * m_cam.m_imgPlaneDistance / ((dirIn.z * dirIn.z * dirIn.z) * m_cam.m_imgPlaneArea);
        return valid ? dirPDF : 0.0f;
    }    
}

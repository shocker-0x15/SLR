//
//  PerspectiveCamera.cpp
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "PerspectiveCamera.h"
#include "../Core/Transform.h"
#include "../MemoryAllocators/ArenaAllocator.h"
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
        
        Normal3D geometricNormal = staticTF * Normal3D(0, 0, 1);
        
        ReferenceFrame shadingFrame;
        shadingFrame.z = (Vector3D)geometricNormal;
        shadingFrame.x = staticTF * Vector3D(1, 0, 0);// assume the transform doesn't include scaling.
        shadingFrame.y = cross(shadingFrame.z, shadingFrame.x);
        
        result->surfPt = SurfacePoint(staticTF * orgLocal, // - position in world coordinate
                                      false, // --------------- atInfinity
                                      shadingFrame, // -------- shading frame
                                      geometricNormal, // ----- geometric normal in world coordinate
                                      lx, ly, // -------------- surface parameter
                                      TexCoord2D::Zero, // ---- texture coordinate
                                      Vector3D::Zero // ------- direction of texture coordinate 0
                                      );
        result->surfPt.setObject(nullptr);
        
        result->posType = m_lensRadius > 0.0f ? DirectionType::LowFreq : DirectionType::Delta0D;
        result->areaPDF = m_lensRadius > 0.0f ? 1.0f / (M_PI * m_lensRadius * m_lensRadius) : 1.0f;
        
        return SampledSpectrum::One;
    }
    
    IDF* PerspectiveCamera::createIDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const {
        float u, v;
        surfPt.getSurfaceParameter(&u, &v);
        return mem.create<PerspectiveIDF>(*this, Point3D(m_lensRadius * u, m_lensRadius * v, 0.0f));
    }
    
    SampledSpectrum PerspectiveIDF::sample(const IDFSample &smp, IDFQueryResult *result) const {
        Point3D pFocus = Point3D(m_cam.m_opWidth * (0.5f - smp.uDir[0]),
                                 m_cam.m_opHeight * (0.5f - smp.uDir[1]),
                                 m_cam.m_objPlaneDistance);
        
        Vector3D dirLocal = normalize(pFocus - m_orgLocal);
        result->dirLocal = dirLocal;
        result->dirPDF = m_cam.m_imgPlaneDistance * m_cam.m_imgPlaneDistance / ((dirLocal.z * dirLocal.z * dirLocal.z) * m_cam.m_imgPlaneArea);
        result->dirType = m_type;
        
        return SampledSpectrum::One;
    }
    
    SampledSpectrum PerspectiveIDF::evaluate(const Vector3D &dirIn) const {
        Point3D pFocas = dirIn * (m_cam.m_objPlaneDistance / dirIn.z) + m_orgLocal;
        bool valid = (pFocas.x >= -m_cam.m_opWidth * 0.5f && pFocas.x <= m_cam.m_opWidth * 0.5f &&
                      pFocas.y >= -m_cam.m_opHeight * 0.5f && pFocas.y <= m_cam.m_opHeight * 0.5f &&
                      dirIn.z >= 0);
        return valid ? SampledSpectrum::One : SampledSpectrum::Zero;
    }
    
    float PerspectiveIDF::evaluatePDF(const Vector3D &dirIn) const {
        Point3D pFocas = dirIn * (m_cam.m_objPlaneDistance / dirIn.z) + m_orgLocal;
        bool valid = (pFocas.x >= -m_cam.m_opWidth * 0.5f && pFocas.x <= m_cam.m_opWidth * 0.5f &&
                      pFocas.y >= -m_cam.m_opHeight * 0.5f && pFocas.y <= m_cam.m_opHeight * 0.5f &&
                      dirIn.z >= 0);
        float dirPDF = m_cam.m_imgPlaneDistance * m_cam.m_imgPlaneDistance / ((dirIn.z * dirIn.z * dirIn.z) * m_cam.m_imgPlaneArea);
        return valid ? dirPDF : 0.0f;
    }
    
    void PerspectiveIDF::calculatePixel(const Vector3D &dirIn, float* hitPx, float* hitPy) const {
        Point3D pFocas = dirIn * (m_cam.m_objPlaneDistance / dirIn.z) + m_orgLocal;
        float smpX = 0.5f - pFocas.x / m_cam.m_opWidth;
        float smpY = 0.5f - pFocas.y / m_cam.m_opHeight;
        *hitPx = m_cam.m_sensor->width() * smpX;
        *hitPy = m_cam.m_sensor->height() * smpY;
    }
}

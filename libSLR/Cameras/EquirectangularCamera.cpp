//
//  EquirectangularCamera.cpp
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "EquirectangularCamera.h"
#include "../Core/Transform.h"
#include "../Memory/ArenaAllocator.h"
#include "../Core/ImageSensor.h"

namespace SLR {
    EquirectangularCamera::EquirectangularCamera(float sensitivity,
                                                 float phiAngle, float thetaAngle) :
    Camera(nullptr),
    m_phiAngle(phiAngle), m_thetaAngle(thetaAngle) {
        m_sensor = createUnique<ImageSensor>(sensitivity > 0 ? sensitivity : 1.0f);
    }
    
    EquirectangularCamera::~EquirectangularCamera() {
    }
    
    ImageSensor* EquirectangularCamera::getSensor() const {
        return m_sensor.get();
    }
    
    SampledSpectrum EquirectangularCamera::sample(const LensPosQuery &query, const LensPosSample &smp, LensPosQueryResult* result) const {
        StaticTransform staticTF;
        if (m_transform)
            m_transform->sample(query.time, &staticTF);
        
        SurfacePoint &surfPt = result->surfPt;
        surfPt.p = staticTF * Point3D::Zero;
        surfPt.gNormal = staticTF * Normal3D(0, 0, 1);
        surfPt.u = 0;
        surfPt.v = 0;
        surfPt.texCoord = TexCoord2D::Zero;
        surfPt.texCoord0Dir = Vector3D::Zero;
        surfPt.shadingFrame.z = (Vector3D)surfPt.gNormal;
        surfPt.shadingFrame.x = staticTF * Vector3D(1, 0, 0);// assume the transform doesn't include scaling.
        surfPt.shadingFrame.y = cross(surfPt.shadingFrame.z, surfPt.shadingFrame.x);
        surfPt.obj = nullptr;
        result->posType = DirectionType::Delta0D;
        result->areaPDF = 1.0f;
        
        return SampledSpectrum::One;
    }
    
    IDF* EquirectangularCamera::createIDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const {
        return mem.create<EquirectangularIDF>(*this);
    }
    
    SampledSpectrum EquirectangularIDF::sample(const IDFSample &smp, IDFQueryResult *result) const {
        float phi = m_cam.m_phiAngle * smp.uDir[0];
        float theta = m_cam.m_thetaAngle * smp.uDir[1];
        result->dirLocal = Vector3D::fromPolarYUp(phi, theta);
        float sinTheta = (1.0f - result->dirLocal.y * result->dirLocal.y);
        result->dirPDF = 1.0f / (m_cam.m_phiAngle * m_cam.m_thetaAngle * sinTheta);
        result->dirType = m_type;
        
        return SampledSpectrum::One;
    }
    
    SampledSpectrum EquirectangularIDF::evaluate(const Vector3D &dirIn) const {
        float phi, theta;
        dirIn.toPolarYUp(&theta, &phi);
        bool valid = (phi >= -m_cam.m_phiAngle * 0.5f && phi < m_cam.m_phiAngle * 0.5f &&
                      theta >= -m_cam.m_thetaAngle * 0.5f && theta < m_cam.m_thetaAngle * 0.5f);
        return valid ? SampledSpectrum::One : SampledSpectrum::Zero;
    }
    
    float EquirectangularIDF::evaluatePDF(const Vector3D &dirIn) const {
        float phi, theta;
        dirIn.toPolarYUp(&theta, &phi);
        bool valid = (phi >= -m_cam.m_phiAngle * 0.5f && phi < m_cam.m_phiAngle * 0.5f &&
                      theta >= -m_cam.m_thetaAngle * 0.5f && theta < m_cam.m_thetaAngle * 0.5f);
        float sinTheta = (1.0f - dirIn.y * dirIn.y);
        float dirPDF = 1.0f / (m_cam.m_phiAngle * m_cam.m_thetaAngle * sinTheta);
        return valid ? dirPDF : 0.0f;
    }
    
    void EquirectangularIDF::calculatePixel(const Vector3D &dirIn, float* hitPx, float* hitPy) const {
        float phi, theta;
        dirIn.toPolarYUp(&theta, &phi);
        float smpX = phi / m_cam.m_phiAngle;
        float smpY = theta / m_cam.m_thetaAngle;
        *hitPx = smpX * m_cam.m_sensor->width();
        *hitPy = smpY * m_cam.m_sensor->height();
    }
}

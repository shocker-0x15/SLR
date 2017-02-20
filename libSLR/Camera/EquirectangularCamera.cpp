//
//  EquirectangularCamera.cpp
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "EquirectangularCamera.h"

#include "../MemoryAllocators/ArenaAllocator.h"
#include "../Core/transform.h"
#include "../Core/ImageSensor.h"

namespace SLR {
    EquirectangularCamera::EquirectangularCamera(float sensitivity, float phiAngle, float thetaAngle) :
    Camera(nullptr),
    m_phiAngle(phiAngle), m_thetaAngle(thetaAngle) {
        m_sensor = new ImageSensor(sensitivity > 0 ? sensitivity : 1.0f);
    }
    
    EquirectangularCamera::~EquirectangularCamera() {
        delete m_sensor;
    }
    
    SampledSpectrum EquirectangularCamera::sample(const LensPosQuery &query, const LensPosSample &smp, LensPosQueryResult* result) const {
        StaticTransform staticTF;
        if (m_transform)
            m_transform->sample(query.time, &staticTF);
        
        Normal3D geometricNormal = staticTF * Normal3D(0, 0, 1);
        
        ReferenceFrame shadingFrame;
        shadingFrame.z = (Vector3D)geometricNormal;
        shadingFrame.x = staticTF * Vector3D(1, 0, 0);// assume the transform doesn't include scaling.
        shadingFrame.y = cross(shadingFrame.z, shadingFrame.x);
        
        result->surfPt = SurfacePoint(staticTF * Point3D::Zero, // - position in world coordinate
                                      false, // -------------------- atInifnity
                                      shadingFrame, // ------------- shading frame
                                      geometricNormal, // ---------- geometric normal in world coordinate
                                      0, 0, // --------------------- surface parameter
                                      TexCoord2D::Zero, // --------- texture coordinate
                                      Vector3D::Zero // ------------ direction of texture coordinate 0
                                      );
        result->surfPt.setObject(nullptr);
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

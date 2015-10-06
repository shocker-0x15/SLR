//
//  PerspectiveCamera.h
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__PerspectiveCamera__
#define __SLR__PerspectiveCamera__

#include "../defines.h"
#include "../references.h"
#include "../Core/cameras.h"
#include "../Core/directional_distribution_functions.h"

namespace SLR {
    class PerspectiveCamera : public Camera {
        float m_aspect;
        float m_fovY;
        float m_objPlaneDistance;
        float m_imgPlaneDistance;
        float m_lensRadius;
        
        float m_opWidth;
        float m_opHeight;
        float m_imgPlaneArea;
        
        friend class PerspectiveIDF;
    public:
        PerspectiveCamera(float aspect, float fovY, float lensRadius, float imgPDist, float objPDist) :
        Camera(nullptr),
        m_aspect(aspect), m_fovY(fovY), m_lensRadius(lensRadius), m_imgPlaneDistance(imgPDist), m_objPlaneDistance(objPDist) {
            m_opHeight = 2.0f * m_objPlaneDistance * std::tan(m_fovY * 0.5f);
            m_opWidth = m_opHeight * m_aspect;
            m_imgPlaneArea = m_opWidth * m_opHeight * std::pow(m_imgPlaneDistance / m_objPlaneDistance, 2);
        };
        
        SampledSpectrum sample(const LensPosQuery &query, const LensPosSample &smp, LensPosQueryResult* result) const override;
        IDF* createIDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const override;
    };
    
    class PerspectiveIDF : public IDF {
        const PerspectiveCamera &m_cam;
        Point3D m_orgLocal;
    public:
        PerspectiveIDF(const PerspectiveCamera &cam, const Point3D &orgLocal) : IDF(DirectionType::WholeSphere | DirectionType::LowFreq), m_cam(cam), m_orgLocal(orgLocal) { };
        
        SampledSpectrum sample(const IDFSample &smp, IDFQueryResult* result) const override;
        SampledSpectrum evaluate(const Vector3D &dirIn) const override;
        float evaluatePDF(const Vector3D &dirIn) const override;
    };    
}

#endif /* defined(__SLR__PerspectiveCamera__) */
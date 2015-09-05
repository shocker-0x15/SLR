//
//  EquirectangularCamera.h
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__EquirectangularCamera__
#define __SLR__EquirectangularCamera__

#include "../defines.h"
#include "../references.h"
#include "../Core/cameras.h"
#include "../Core/directional_distribution_functions.h"

class EquirectangularCamera : public Camera {
    float m_phiAngle;
    float m_thetaAngle;
    friend class EquirectangularIDF;
public:
    EquirectangularCamera(float phiAngle, float thetaAngle) :
    Camera(nullptr),
    m_phiAngle(phiAngle), m_thetaAngle(thetaAngle) { };
    
    Spectrum sample(const LensPosQuery &query, const LensPosSample &smp, LensPosQueryResult* result) const override;
    IDF* createIDF(const SurfacePoint &surfPt, ArenaAllocator &mem) const override;
};

class EquirectangularIDF : public IDF {
    const EquirectangularCamera &m_cam;
public:
    EquirectangularIDF(const EquirectangularCamera &cam) : IDF(DirectionType::WholeSphere | DirectionType::LowFreq), m_cam(cam) { };
    
    Spectrum sample(const IDFSample &smp, IDFQueryResult* result) const override;
    Spectrum evaluate(const Vector3D &dirIn) const override;
    float evaluatePDF(const Vector3D &dirIn) const override;
};

#endif /* defined(__SLR__EquirectangularCamera__) */
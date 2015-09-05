//
//  camera_nodes.h
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__camera_nodes__
#define __SLR__camera_nodes__

#include "../defines.h"
#include "../references.h"
#include "nodes.h"

class PerspectiveCameraNode : public CameraNode {
    float m_aspect;
    float m_fovY;
    float m_objPlaneDistance;
    float m_imgPlaneDistance;
    float m_lensRadius;
public:
    PerspectiveCameraNode(float aspect, float fovY, float lensRadius, float imgPDist, float objPDist) :
    m_aspect(aspect), m_fovY(fovY), m_lensRadius(lensRadius), m_imgPlaneDistance(imgPDist), m_objPlaneDistance(objPDist) { };
    
    void applyTransform(const StaticTransform &t) final { };
    void createCamera() final;
};

class EquirectangularCameraNode : public CameraNode {
    float m_phiAngle;
    float m_thetaAngle;
public:
    EquirectangularCameraNode(float phiAngle, float thetaAngle) :
    m_phiAngle(phiAngle), m_thetaAngle(thetaAngle) { };
    
    void applyTransform(const StaticTransform &t) final { };
    void createCamera() final;
};

#endif /* defined(__SLR__camera_nodes__) */
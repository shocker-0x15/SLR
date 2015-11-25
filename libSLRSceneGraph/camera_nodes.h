//
//  camera_nodes.h
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLRSceneGraph__camera_nodes__
#define __SLRSceneGraph__camera_nodes__

#include <libSLR/defines.h>
#include "references.h"
#include "nodes.h"

namespace SLRSceneGraph {
    class PerspectiveCameraNode : public CameraNode {
        float m_sensitivity;
        float m_aspect;
        float m_fovY;
        float m_objPlaneDistance;
        float m_imgPlaneDistance;
        float m_lensRadius;
    public:
        PerspectiveCameraNode(float sensitivity, float aspect, float fovY, float lensRadius, float imgPDist, float objPDist) :
        m_sensitivity(sensitivity), m_aspect(aspect), m_fovY(fovY), m_lensRadius(lensRadius), m_imgPlaneDistance(imgPDist), m_objPlaneDistance(objPDist) { };
        
        void applyTransform(const SLR::StaticTransform &t) final { };
        void createCamera() final;
    };
    
    class EquirectangularCameraNode : public CameraNode {
        float m_sensitivity;
        float m_phiAngle;
        float m_thetaAngle;
    public:
        EquirectangularCameraNode(float sensitivity, float phiAngle, float thetaAngle) :
        m_sensitivity(sensitivity), m_phiAngle(phiAngle), m_thetaAngle(thetaAngle) { };
        
        void applyTransform(const SLR::StaticTransform &t) final { };
        void createCamera() final;
    };    
}

#endif /* defined(__SLRSceneGraph__camera_nodes__) */

//
//  camera_nodes.cpp
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "camera_nodes.h"
#include <libSLR/Scene/camera_nodes.h>

namespace SLRSceneGraph {
    void PerspectiveCameraNode::setupRawData() {
        m_rawData = new SLR::PerspectiveCameraNode(m_sensitivity, m_aspect, m_fovY, m_lensRadius, m_imgPlaneDistance, m_objPlaneDistance);
    }
    
    PerspectiveCameraNode::PerspectiveCameraNode(float sensitivity, float aspect, float fovY, float lensRadius, float imgPDist, float objPDist) :
    m_sensitivity(sensitivity), m_aspect(aspect), m_fovY(fovY), m_lensRadius(lensRadius), m_imgPlaneDistance(imgPDist), m_objPlaneDistance(objPDist) {
        setupRawData();
    }
    
    NodeRef PerspectiveCameraNode::copy() const {
        PerspectiveCameraNodeRef ret = createShared<PerspectiveCameraNode>(m_sensitivity, m_aspect, m_fovY, m_lensRadius, m_imgPlaneDistance, m_objPlaneDistance);
        return ret;
    }
    
    void PerspectiveCameraNode::prepareForRendering() {
        new (m_rawData) SLR::PerspectiveCameraNode(m_sensitivity, m_aspect, m_fovY, m_lensRadius, m_imgPlaneDistance, m_objPlaneDistance);
    }
    
    
    
    void EquirectangularCameraNode::setupRawData() {
        m_rawData = new SLR::EquirectangularCameraNode(m_sensitivity, m_phiAngle, m_thetaAngle);
    }
    
    EquirectangularCameraNode::EquirectangularCameraNode(float sensitivity, float phiAngle, float thetaAngle) :
    m_sensitivity(sensitivity), m_phiAngle(phiAngle), m_thetaAngle(thetaAngle) {
        setupRawData();
    }
    
    NodeRef EquirectangularCameraNode::copy() const {
        EquirectangularCameraNodeRef ret = createShared<EquirectangularCameraNode>(m_sensitivity, m_phiAngle, m_thetaAngle);
        return ret;
    }
    
    void EquirectangularCameraNode::prepareForRendering() {
        new (m_rawData) SLR::EquirectangularCameraNode(m_sensitivity, m_phiAngle, m_thetaAngle);
    }
}

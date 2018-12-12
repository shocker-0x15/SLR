//
//  camera_nodes.cpp
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "camera_nodes.h"

#include <libSLR/Scene/camera_nodes.h>

namespace SLRSceneGraph {
    void PerspectiveCameraNode::allocateRawData() {
        m_rawData = (SLR::Node*)malloc(sizeof(SLR::PerspectiveCameraNode));
    }
    
    void PerspectiveCameraNode::setupRawData() {
        new (m_rawData) SLR::PerspectiveCameraNode(m_sensitivity, m_aspect, m_fovY, m_lensRadius, m_imgPlaneDistance, m_objPlaneDistance);
        m_setup = true;
    }
    
    void PerspectiveCameraNode::terminateRawData() {
        SLR::PerspectiveCameraNode &raw = *(SLR::PerspectiveCameraNode*)m_rawData;
        if (m_setup)
            raw.~PerspectiveCameraNode();
        m_setup = false;
    }
    
    PerspectiveCameraNode::PerspectiveCameraNode(float sensitivity, float aspect, float fovY, float lensRadius, float imgPDist, float objPDist) :
    m_sensitivity(sensitivity), m_aspect(aspect), m_fovY(fovY), m_lensRadius(lensRadius), m_imgPlaneDistance(imgPDist), m_objPlaneDistance(objPDist) {
        allocateRawData();
    }
    
    NodeRef PerspectiveCameraNode::copy() const {
        return createShared<PerspectiveCameraNode>(m_sensitivity, m_aspect, m_fovY, m_lensRadius, m_imgPlaneDistance, m_objPlaneDistance);
    }
    
    void PerspectiveCameraNode::prepareForRendering() {
        terminateRawData();
        setupRawData();
    }
    
    
    
    void EquirectangularCameraNode::allocateRawData() {
        m_rawData = (SLR::Node*)malloc(sizeof(SLR::EquirectangularCameraNode));
    }
    
    void EquirectangularCameraNode::setupRawData() {
        new (m_rawData) SLR::EquirectangularCameraNode(m_sensitivity, m_phiAngle, m_thetaAngle);
        m_setup = true;
    }
    
    void EquirectangularCameraNode::terminateRawData() {
        SLR::EquirectangularCameraNode &raw = *(SLR::EquirectangularCameraNode*)m_rawData;
        if (m_setup)
            raw.~EquirectangularCameraNode();
        m_setup = false;
    }
    
    EquirectangularCameraNode::EquirectangularCameraNode(float sensitivity, float phiAngle, float thetaAngle) :
    m_sensitivity(sensitivity), m_phiAngle(phiAngle), m_thetaAngle(thetaAngle) {
        allocateRawData();
    }
    
    NodeRef EquirectangularCameraNode::copy() const {
        return createShared<EquirectangularCameraNode>(m_sensitivity, m_phiAngle, m_thetaAngle);
    }
    
    void EquirectangularCameraNode::prepareForRendering() {
        terminateRawData();
        setupRawData();
    }
}

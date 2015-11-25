//
//  camera_nodes.cpp
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "camera_nodes.h"
#include <libSLR/Cameras/PerspectiveCamera.h>
#include <libSLR/Cameras/EquirectangularCamera.h>

namespace SLRSceneGraph {
    void PerspectiveCameraNode::createCamera() {
        m_camera = new SLR::PerspectiveCamera(m_sensitivity, m_aspect, m_fovY, m_lensRadius, m_imgPlaneDistance, m_objPlaneDistance);
    }
    
    void EquirectangularCameraNode::createCamera() {
        m_camera = new SLR::EquirectangularCamera(m_sensitivity, m_phiAngle, m_thetaAngle);
    }    
}

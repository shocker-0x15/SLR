//
//  camera_nodes.cpp
//
//  Created by 渡部 心 on 2017/01/06.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "camera_nodes.h"

#include "../MemoryAllocators/ArenaAllocator.h"
#include "../Camera/PerspectiveCamera.h"
#include "../Camera/EquirectangularCamera.h"

namespace SLR {
    void PerspectiveCameraNode::createRenderingData(Allocator* mem, const Transform* subTF, RenderingData* data) {
        m_camera = mem->create<PerspectiveCamera>(m_sensitivity, m_aspect, m_fovY, m_lensRadius, m_imgPlaneDistance, m_objPlaneDistance);
        data->camera = m_camera;
    }
    
    void PerspectiveCameraNode::destroyRenderingData(Allocator* mem) {
        mem->destroy(m_camera);
    }
    
    
    
    void EquirectangularCameraNode::createRenderingData(Allocator* mem, const Transform* subTF, RenderingData* data) {
        m_camera = mem->create<EquirectangularCamera>(m_sensitivity, m_phiAngle, m_thetaAngle);
        data->camera = m_camera;
    }
    
    void EquirectangularCameraNode::destroyRenderingData(Allocator* mem) {
        mem->destroy(m_camera);
    }
}

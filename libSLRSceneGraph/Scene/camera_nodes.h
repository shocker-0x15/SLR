//
//  camera_nodes.h
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLRSceneGraph_camera_nodes__
#define __SLRSceneGraph_camera_nodes__

#include <libSLR/defines.h>
#include "../declarations.h"
#include "node.h"

namespace SLRSceneGraph {
    class SLR_SCENEGRAPH_API PerspectiveCameraNode : public Node {
        float m_sensitivity;
        float m_aspect;
        float m_fovY;
        float m_lensRadius;
        float m_imgPlaneDistance;
        float m_objPlaneDistance;
        
        void allocateRawData() override;
        void setupRawData() override;
        void terminateRawData() override;
    public:
        PerspectiveCameraNode(float sensitivity, float aspect, float fovY, float lensRadius, float imgPDist, float objPDist);
        
        NodeRef copy() const override;
        
        void prepareForRendering() override;
    };
    
    
    
    class SLR_SCENEGRAPH_API EquirectangularCameraNode : public Node {
        float m_sensitivity;
        float m_phiAngle;
        float m_thetaAngle;
        
        void allocateRawData() override;
        void setupRawData() override;
        void terminateRawData() override;
    public:
        EquirectangularCameraNode(float sensitivity, float phiAngle, float thetaAngle);
        
        NodeRef copy() const override;
        
        void prepareForRendering() override;
    };    
}

#endif /* __SLRSceneGraph_camera_nodes__ */

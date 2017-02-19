//
//  camera_nodes.h
//
//  Created by 渡部 心 on 2017/01/06.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef __SLR_camera_nodes__
#define __SLR_camera_nodes__

#include "../defines.h"
#include "../declarations.h"
#include "node.h"

namespace SLR {
    class SLR_API PerspectiveCameraNode : public Node {
        float m_sensitivity;
        float m_aspect;
        float m_fovY;
        float m_lensRadius;
        float m_imgPlaneDistance;
        float m_objPlaneDistance;
        
        PerspectiveCamera* m_camera;
    public:
        PerspectiveCameraNode(float sensitivity, float aspect, float fovY, float lensRadius, float imgPDist, float objPDist) :
        m_sensitivity(sensitivity), m_aspect(aspect), m_fovY(fovY), m_lensRadius(lensRadius), m_imgPlaneDistance(imgPDist), m_objPlaneDistance(objPDist),
        m_camera(nullptr) { }
        
        bool isDirectlyTransformable() const override { return false; }
        void createRenderingData(Allocator* mem, const Transform* subTF, RenderingData *data) override;
        void destroyRenderingData(Allocator* mem) override;
    };
    
    
    
    class SLR_API EquirectangularCameraNode : public Node {
        float m_sensitivity;
        float m_phiAngle;
        float m_thetaAngle;
        
        EquirectangularCamera* m_camera;
    public:
        EquirectangularCameraNode(float sensitivity, float phiAngle, float thetaAngle) :
        m_sensitivity(sensitivity), m_phiAngle(phiAngle), m_thetaAngle(thetaAngle),
        m_camera(nullptr) { }
        
        bool isDirectlyTransformable() const override { return false; }
        void createRenderingData(Allocator* mem, const Transform* subTF, RenderingData *data) override;
        void destroyRenderingData(Allocator* mem) override;
    };
}

#endif /* __SLR_camera_nodes__ */

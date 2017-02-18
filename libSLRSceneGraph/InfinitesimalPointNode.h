//
//  InfinitesimalPointNode.h
//
//  Created by 渡部 心 on 2017/02/17.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef __SLRSceneGraph_InfinitesimalPointNode__
#define __SLRSceneGraph_InfinitesimalPointNode__

#include <libSLR/defines.h>
#include "declarations.h"
#include "nodes.h"
#include <libSLR/Core/geometry.h>

namespace SLRSceneGraph {
    class SLR_SCENEGRAPH_API InfinitesimalPointNode : public SurfaceNode {
        SLR::Point3D m_position;
        SLR::Vector3D m_direction;
        SurfaceMaterialRef m_material;
        
        void allocateRawData() override;
        void setupRawData() override;
        void terminateRawData() override;
    public:
        InfinitesimalPointNode(const SLR::Point3D &p, const SLR::Vector3D &d, const SurfaceMaterialRef &mat);
        
        NodeRef copy() const override;
        
        void prepareForRendering() override;
    };
}

#endif /* __SLRSceneGraph_InfinitesimalPointNode__ */

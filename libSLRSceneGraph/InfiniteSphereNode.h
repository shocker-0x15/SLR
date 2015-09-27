//
//  InfiniteSphereNode.h
//  SLR
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLRSceneGraph__InfiniteSphereNode__
#define __SLRSceneGraph__InfiniteSphereNode__

#include <libSLR/defines.h>
#include <libSLR/references.h>
#include "references.h"
#include "nodes.h"
#include <libSLR/Surface/InfiniteSphere.h>

namespace SLRSceneGraph {
    class InfiniteSphereNode : public Node {
        const SLR::Scene* m_scene;
        SLR::InfiniteSphere m_sphere;
        std::shared_ptr<SLR::IBLEmission> m_emitter;
        SLR::SurfaceMaterialRef m_surfMat;
        SLR::InfiniteSphereSurfaceObject* m_surfObj;
        SLR::RegularConstantContinuous2D* m_dist;
        bool m_ready;
    public:
        InfiniteSphereNode(const SLR::Scene* scene, const SLR::SpectrumTextureRef &IBLTex, float scale);
        ~InfiniteSphereNode();
        
        void getRenderingData(SLR::ArenaAllocator &mem, const SLR::Transform* subTF, RenderingData *data) {
            SLRAssert(false, "InfiniteSphereNode::getRenderingData() should not be called.");
        };
        SLR::InfiniteSphereSurfaceObject* getSurfaceObject();
    };    
}

#endif /* defined(__SLRSceneGraph__InfiniteSphereNode__) */

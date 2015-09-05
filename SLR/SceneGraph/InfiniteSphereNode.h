//
//  InfiniteSphereNode.h
//  SLR
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__InfiniteSphereNode__
#define __SLR__InfiniteSphereNode__

#include "../defines.h"
#include "../references.h"
#include "nodes.h"
#include "../Surface/InfiniteSphere.h"

class InfiniteSphereNode : public Node {
    const Scene* m_scene;
    InfiniteSphere m_sphere;
    std::shared_ptr<IBLEmission> m_emitter;
    SurfaceMaterialRef m_surfMat;
    InfiniteSphereSurfaceObject* m_surfObj;
    RegularConstantContinuous2D* m_dist;
    bool m_ready;
public:
    InfiniteSphereNode(const Scene* scene, const SpectrumTextureRef &IBLTex, const Spectrum &scale);
    ~InfiniteSphereNode();
    
    void getRenderingData(ArenaAllocator &mem, const Transform* subTF, RenderingData *data) {
        SLRAssert(false, "InfiniteSphereNode::getRenderingData() should not be called.");
    };
    InfiniteSphereSurfaceObject* getSurfaceObject();
};

#endif /* defined(__SLR__InfiniteSphereNode__) */

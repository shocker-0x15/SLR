//
//  InfiniteSphereNode.h
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLRSceneGraph_InfiniteSphereNode__
#define __SLRSceneGraph_InfiniteSphereNode__

#include <libSLR/defines.h>
#include "references.h"
#include "nodes.h"

namespace SLRSceneGraph {
    class SLR_SCENEGRAPH_API InfiniteSphereNode : public SurfaceNode {
        SceneWRef m_scene;
        SpectrumTextureRef m_IBLTex;
        float m_scale;
        
        void allocateRawData() override;
        void setupRawData() override;
        void terminateRawData() override;
    public:
        InfiniteSphereNode(const SceneWRef &scene, const SpectrumTextureRef &IBLTex, float scale);
        
        NodeRef copy() const override { SLRAssert_ShouldNotBeCalled(); return nullptr; }
        
        void prepareForRendering() override;
    };
}

#endif /* __SLRSceneGraph_InfiniteSphereNode__ */

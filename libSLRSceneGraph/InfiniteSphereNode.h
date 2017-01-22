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
#include <libSLR/Surface/InfiniteSphere.h>

namespace SLRSceneGraph {
    class SLR_SCENEGRAPH_API InfiniteSphereNode : public Node {
        SceneWRef m_scene;
        SpectrumTextureRef m_IBLTex;
        float m_scale;
        
        void setupRawData() override;
    public:
        InfiniteSphereNode(const SceneWRef &scene, const SpectrumTextureRef &IBLTex, float scale);
        
        NodeRef copy() const override { SLRAssert_ShouldNotBeCalled(); return nullptr; }
        
        void prepareForRendering() override;
    };
}

#endif /* defined(__SLRSceneGraph__InfiniteSphereNode__) */

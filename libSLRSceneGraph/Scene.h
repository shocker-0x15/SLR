//
//  Scene.h
//
//  Created by 渡部 心 on 2015/07/11.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLRSceneGraph__Scene__
#define __SLRSceneGraph__Scene__

#include <libSLR/defines.h>
#include "references.h"

namespace SLRSceneGraph {
    class Scene {
        InternalNodeRef m_rootNode;
        InfiniteSphereNodeRef m_envNode;
    public:
        Scene();
        ~Scene();
        
        InternalNodeRef &rootNode() { return m_rootNode; };
        void setEnvNode(const InfiniteSphereNodeRef &node) { m_envNode = node; };
        
        void build(SLR::Scene** scene, SLR::ArenaAllocator &mem);
    };
    
    struct RenderingContext {
        std::unique_ptr<SLR::Renderer> renderer;
        int32_t width;
        int32_t height;
        float timeStart;
        float timeEnd;
        float brightness;
        int32_t rngSeed;
        
        RenderingContext();
        ~RenderingContext();
        
        RenderingContext &operator=(RenderingContext &&ctx);
    };
}

#endif

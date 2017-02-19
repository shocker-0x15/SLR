//
//  Scene.h
//
//  Created by 渡部 心 on 2015/07/11.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLRSceneGraph_Scene__
#define __SLRSceneGraph_Scene__

#include <libSLR/defines.h>
#include "../declarations.h"

namespace SLRSceneGraph {
    class SLR_SCENEGRAPH_API Scene {
        SLR::Scene* m_raw;
        InternalNodeRef m_rootNode;
        InfiniteSphereNodeRef m_envNode;
    public:
        Scene();
        ~Scene();
        
        InternalNodeRef &rootNode() {
            return m_rootNode;
        }
        void setEnvironmentNode(const InfiniteSphereNodeRef &node);
        
        void prepareForRendering();
        SLR::Scene* getRaw() const { return m_raw; }
    };
    
    struct SLR_SCENEGRAPH_API RenderingContext {
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

#endif /* __SLRSceneGraph_Scene__ */

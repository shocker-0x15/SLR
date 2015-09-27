//
//  Renderer.h
//
//  Created by 渡部 心 on 2015/07/11.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef SLR_Renderer_h
#define SLR_Renderer_h

#include "../defines.h"
#include "../references.h"

namespace SLR {
    class Renderer {
    public:
        virtual void render(const Scene &scene, const RenderSettings &settings) const = 0;
    };
}

#endif

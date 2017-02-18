//
//  Renderer.h
//
//  Created by 渡部 心 on 2015/07/11.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_Renderer__
#define __SLR_Renderer__

#include "../defines.h"
#include "../declarations.h"

namespace SLR {
    class SLR_API Renderer {
    public:
        virtual ~Renderer() { };
        virtual void render(const Scene &scene, const RenderSettings &settings) const = 0;
    };
}

#endif /* __SLR_Renderer__ */

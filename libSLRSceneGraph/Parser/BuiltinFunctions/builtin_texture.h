//
//  builtin_texture.h
//
//  Created by 渡部 心 on 2016/08/20.
//  Copyright c 2016年 渡部 心. All rights reserved.
//

#ifndef SLRSceneGraph_builtin_texture_hpp
#define SLRSceneGraph_builtin_texture_hpp

#include <libSLR/defines.h>
#include "../../declarations.h"

#include "../SceneParser.h"

namespace SLRSceneGraph {
    namespace BuiltinFunctions {
        namespace Texture {
            extern const Element Texture2DMapping;
            extern const Element Texture3DMapping;
            extern const Element SpectrumTexture;
            extern const Element NormalTexture;
            extern const Element FloatTexture;
        }
    }
}

#endif /* builtin_texture_hpp */

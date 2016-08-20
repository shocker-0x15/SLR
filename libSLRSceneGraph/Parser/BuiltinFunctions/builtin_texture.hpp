//
//  builtin_texture.hpp
//
//  Created by 渡部 心 on 2016/08/20.
//  Copyright © 2016年 渡部 心. All rights reserved.
//

#ifndef SLRSceneGraph_builtin_texture_hpp
#define SLRSceneGraph_builtin_texture_hpp

#include <libSLR/defines.h>
#include "../../references.h"

#include "../SceneParser.hpp"

namespace SLRSceneGraph {
    namespace BuiltinFunctions {
        namespace Texture {
            extern const Function Texture2DMapping;
            extern const Function Texture3DMapping;
            extern const Function SpectrumTexture;
            extern const Function NormalTexture;
            extern const Function FloatTexture;
        }
    }
}

#endif /* builtin_texture_hpp */

//
//  builtin_transform.hpp
//
//  Created by 渡部 心 on 2016/08/20.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef SLRSceneGraph_builtin_transform_hpp
#define SLRSceneGraph_builtin_transform_hpp

#include <libSLR/defines.h>
#include "../../references.h"

#include "../SceneParser.hpp"

namespace SLRSceneGraph {
    namespace BuiltinFunctions {
        namespace Transform {
            extern const Function translate;
            extern const Function rotate;
            extern const Function rotateX;
            extern const Function rotateY;
            extern const Function rotateZ;
            extern const Function scale;
            extern const Function lookAt;
            extern const Function AnimatedTransform;
        }
    }
}

#endif /* transform_hpp */

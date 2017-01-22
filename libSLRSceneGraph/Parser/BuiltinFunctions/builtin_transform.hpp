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
            extern const Element translate;
            extern const Element rotate;
            extern const Element rotateX;
            extern const Element rotateY;
            extern const Element rotateZ;
            extern const Element scale;
            extern const Element lookAt;
            extern const Element AnimatedTransform;
        }
    }
}

#endif /* transform_hpp */

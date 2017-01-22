//
//  builtin_math.hpp
//
//  Created by 渡部 心 on 2016/08/20.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef SLRSceneGraph_builtin_math_hpp
#define SLRSceneGraph_builtin_math_hpp

#include <libSLR/defines.h>
#include "../../references.h"

#include "../SceneParser.hpp"

namespace SLRSceneGraph {
    namespace BuiltinFunctions {
        namespace Math {
            extern const Element min;
            extern const Element max;
            extern const Element clamp;
            extern const Element sqrt;
            extern const Element pow;
            extern const Element sin;
            extern const Element cos;
            extern const Element tan;
            extern const Element asin;
            extern const Element acos;
            extern const Element atan;
            extern const Element dot;
            extern const Element cross;
            extern const Element distance;
        }
    }
}

#endif /* math_hpp */

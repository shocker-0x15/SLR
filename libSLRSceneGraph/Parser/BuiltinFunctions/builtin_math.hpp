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
            extern const Function min;
            extern const Function max;
            extern const Function clamp;
            extern const Function sqrt;
            extern const Function pow;
            extern const Function sin;
            extern const Function cos;
            extern const Function tan;
            extern const Function asin;
            extern const Function acos;
            extern const Function atan;
            extern const Function dot;
            extern const Function cross;
        }
    }
}

#endif /* math_hpp */

//
//  TexCoord2.cpp
//
//  Created by 渡部 心 on 2016/01/26.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "TexCoord2.h"

namespace SLR {
    template struct SLR_API TexCoord2Template<float>;
//    template SLR_API const TexCoord2Template<float> TexCoord2Template<float>::Zero;
//    template SLR_API const TexCoord2Template<float> TexCoord2Template<float>::One;
//    template SLR_API const TexCoord2Template<float> TexCoord2Template<float>::Inf;
//    template SLR_API const TexCoord2Template<float> TexCoord2Template<float>::NaN;

    template struct SLR_API TexCoord2Template<double>;
//    template SLR_API const TexCoord2Template<double> TexCoord2Template<double>::Zero;
//    template SLR_API const TexCoord2Template<double> TexCoord2Template<double>::One;
//    template SLR_API const TexCoord2Template<double> TexCoord2Template<double>::Inf;
//    template SLR_API const TexCoord2Template<double> TexCoord2Template<double>::NaN;
}

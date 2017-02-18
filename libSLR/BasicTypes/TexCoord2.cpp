//
//  TexCoord2.cpp
//
//  Created by 渡部 心 on 2016/01/26.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "TexCoord2.h"

namespace SLR {
    template struct SLR_API TexCoord2DTemplate<float>;
//    template SLR_API const TexCoord2DTemplate<float> TexCoord2DTemplate<float>::Zero;
//    template SLR_API const TexCoord2DTemplate<float> TexCoord2DTemplate<float>::One;
//    template SLR_API const TexCoord2DTemplate<float> TexCoord2DTemplate<float>::Inf;
//    template SLR_API const TexCoord2DTemplate<float> TexCoord2DTemplate<float>::NaN;

    template struct SLR_API TexCoord2DTemplate<double>;
//    template SLR_API const TexCoord2DTemplate<double> TexCoord2DTemplate<double>::Zero;
//    template SLR_API const TexCoord2DTemplate<double> TexCoord2DTemplate<double>::One;
//    template SLR_API const TexCoord2DTemplate<double> TexCoord2DTemplate<double>::Inf;
//    template SLR_API const TexCoord2DTemplate<double> TexCoord2DTemplate<double>::NaN;
}

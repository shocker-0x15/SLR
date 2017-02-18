//
//  image_loader.h
//
//  Created by 渡部 心 on 2014/05/05.
//  Copyright (c) 2014年 渡部 心. All rights reserved.
//

#ifndef __SLR_image_loader__
#define __SLR_image_loader__

#include <libSLR/defines.h>
#include <libSLR/declarations.h>
#include <libSLR/BasicTypes/spectrum_base.h>
#include "../declarations.h"
#include <string>
#include <cstdint>

enum class ColorFormat {
    RGB8x3 = 0,
    RGB_8x4,
    RGBA8x4,
    RGBA16Fx4,
    Gray8,
    Num
};

bool getImageInfo(const std::string &filePath, uint32_t* width, uint32_t* height, uint64_t* requiredSize, ColorFormat* color);
bool loadImage(const std::string &filePath, uint8_t* storage, bool gammaCorrection);

#endif /* __SLR_image_loader__ */

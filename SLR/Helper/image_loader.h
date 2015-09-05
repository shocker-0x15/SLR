//
//  image_loader.h
//
//  Created by 渡部 心 on 2014/05/05.
//  Copyright (c) 2014年 渡部 心. All rights reserved.
//

#ifndef __CLeaR__image_loader__
#define __CLeaR__image_loader__

#include <string>
#include <cstdint>
#include <half.h>

enum class ColorFormat {
    RGB8x3 = 0,
    RGB_8x4,
    RGBA8x4,
    RGBA16Fx4,
    Gray8,
    Num
};

struct RGB8x3 { uint8_t r, g, b; };
struct RGB_8x4 { uint8_t r, g, b, dummy; };
struct RGBA8x4  { uint8_t r, g, b, a; };
struct RGBA16Fx4 { half r, g, b, a; };
struct Gray8 { uint8_t v; };

extern const size_t sizesOfColorFormats[(uint32_t)ColorFormat::Num];

bool getImageInfo(const std::string &filePath, uint32_t* width, uint32_t* height, uint64_t* requiredSize, ColorFormat* color);
bool loadImage(const std::string &filePath, uint8_t* storage, bool gammaCorrection);

#endif

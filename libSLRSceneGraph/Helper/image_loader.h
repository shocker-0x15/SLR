//
//  image_loader.h
//
//  Created by 渡部 心 on 2014/05/05.
//  Copyright (c) 2014年 渡部 心. All rights reserved.
//

#ifndef __CLeaR__image_loader__
#define __CLeaR__image_loader__

#include <libSLR/defines.h>
#include <libSLR/references.h>
#include <libSLR/BasicTypes/Spectrum.h>
#include "references.h"
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

namespace SLRSceneGraph {
    extern std::map<std::string, SLR::Image2DRef> s_imageDB;
    
    std::shared_ptr<SLR::TiledImage2D> createTiledImage(const std::string &filepath, SLR::Allocator *mem, SLR::SpectrumType spType, bool gammaCorrection = false);
}

#endif

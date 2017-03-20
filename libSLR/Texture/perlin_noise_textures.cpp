//
//  perlin_noise_textures.cpp
//
//  Created by 渡部 心 on 2017/03/17.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "perlin_noise_textures.h"

#include "../Core/geometry.h"
#include "../RNG/LinearCongruentialRNG.h"

namespace SLR {
    float PerlinNoiseFloatTexture::evaluate(const Point3D &p) const {
        return m_generator.evaluate(p.x, p.y, p.z);
    }
}

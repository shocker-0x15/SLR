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
    Normal3D PerlinNoiseNormalTexture::evaluate(const Point3D &p) const {
        float phi = 2 * M_PI * m_generator[0].evaluate(p);
        float theta = m_thetaMax * m_generator[1].evaluate(p);
        
        return Normal3D::fromPolarZUp(phi, theta);
    }
    
    
    
    float PerlinNoiseFloatTexture::evaluate(const Point3D &p) const {
        return m_generator.evaluate(p);
    }
}

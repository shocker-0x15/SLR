//
//  IBLEmission.cpp
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "IBLEmission.h"
#include "../EDFs/IBLEDF.h"
#include "../Memory/ArenaAllocator.h"
#include "../Core/textures.h"
#include "../Core/SurfaceObject.h"

namespace SLR {
    SampledSpectrum IBLEmission::emittance(const SurfacePoint &surfPt, const WavelengthSamples &wls) const {
        return m_scene->getWorldDiscArea() * m_coeffM->evaluate(surfPt.texCoord, wls) * m_scale;
    }
    
    EDF* IBLEmission::getEDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        return mem.create<IBLEDF>(m_scene->getWorldDiscArea());
    }    
}

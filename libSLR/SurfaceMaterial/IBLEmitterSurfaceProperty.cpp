//
//  IBLEmitterSurfaceProperty.cpp
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "IBLEmitterSurfaceProperty.h"

#include "../MemoryAllocators/ArenaAllocator.h"
#include "../Core/textures.h"
#include "../EDF/IBLEDF.h"
#include "../Scene/Scene.h"

namespace SLR {
    SampledSpectrum IBLEmitterSurfaceProperty::emittance(const SurfacePoint &surfPt, const WavelengthSamples &wls) const {
        return M_PI * m_coeffM->evaluate(surfPt, wls) * m_scale;
    }
    
    EDF* IBLEmitterSurfaceProperty::getEDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        return mem.create<IBLEDF>(m_scene->getWorldDiscArea());
    }
    
    ContinuousDistribution2D* IBLEmitterSurfaceProperty::createIBLImportanceMap() const {
        return m_coeffM->createIBLImportanceMap();
    }
}

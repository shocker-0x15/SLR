//
//  DiffuseEmission.h
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__DiffuseEmission__
#define __SLR__DiffuseEmission__

#include "../defines.h"
#include "../references.h"
#include "../Core/surface_material.h"

class DiffuseEmission : public EmitterSurfaceProperty {
    SpectrumTextureRef m_emittance;
public:
    DiffuseEmission(const SpectrumTextureRef &emittance) :
    m_emittance(emittance) {};
    
    SampledSpectrum emittance(const SurfacePoint &surfPt, const WavelengthSamples &wls) const override;
    EDF* getEDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
};

#endif /* defined(__SLR__DiffuseEmission__) */

//
//  IBLEmission.h
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__IBLEmission__
#define __SLR__IBLEmission__

#include "../defines.h"
#include "../references.h"
#include "../Core/surface_material.h"

class IBLEmission : public EmitterSurfaceProperty {
    const Scene* m_scene;
    SpectrumTextureRef m_coeffM;
    SampledSpectrum m_scale;
public:
    IBLEmission(const Scene* scene, const SpectrumTextureRef &coeffM, const SampledSpectrum &scale) :
    m_scene(scene), m_coeffM(coeffM), m_scale(scale) {};
    
    SampledSpectrum emittance(const SurfacePoint &surfPt, const WavelengthSamples &wls) const override;
    EDF* getEDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
};

#endif /* defined(__SLR__IBLEmission__) */

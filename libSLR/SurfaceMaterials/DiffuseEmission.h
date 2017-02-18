//
//  DiffuseEmission.h
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_DiffuseEmission__
#define __SLR_DiffuseEmission__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/surface_material.h"

namespace SLR {
    class SLR_API DiffuseEmitterSurfaceProperty : public EmitterSurfaceProperty {
        const SpectrumTexture* m_emittance;
    public:
        DiffuseEmitterSurfaceProperty(const SpectrumTexture* emittance) :
        m_emittance(emittance) {};
        
        SampledSpectrum emittance(const SurfacePoint &surfPt, const WavelengthSamples &wls) const override;
        EDF* getEDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
}

#endif /* __SLR_DiffuseEmission__ */

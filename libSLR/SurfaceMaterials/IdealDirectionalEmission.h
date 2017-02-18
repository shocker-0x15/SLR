//
//  IdealDirectionalEmission.h
//
//  Created by 渡部 心 on 2017/02/17.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef __SLR_IdealDirectionalEmission__
#define __SLR_IdealDirectionalEmission__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/surface_material.h"

namespace SLR {
    class SLR_API IdealDirectionalEmitterSurfaceProperty : public EmitterSurfaceProperty {
        const SpectrumTexture* m_emittance;
    public:
        IdealDirectionalEmitterSurfaceProperty(const SpectrumTexture* emittance) :
        m_emittance(emittance) {};
        
        SampledSpectrum emittance(const SurfacePoint &surfPt, const WavelengthSamples &wls) const override;
        EDF* getEDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
}

#endif /* __SLR_IdealDirectionalEmission__ */

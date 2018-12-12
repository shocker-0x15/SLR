//
//  basic_emitter_surface_properties.h
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_basic_emitter_surface_properties__
#define __SLR_basic_emitter_surface_properties__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/geometry.h"
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
    
    
    
    class SLR_API IdealDirectionalEmitterSurfaceProperty : public EmitterSurfaceProperty {
        const SpectrumTexture* m_emittance;
        const Vector3D m_direction;
    public:
        IdealDirectionalEmitterSurfaceProperty(const SpectrumTexture* emittance, const Vector3D &dir) :
        m_emittance(emittance), m_direction(dir) {};
        
        SampledSpectrum emittance(const SurfacePoint &surfPt, const WavelengthSamples &wls) const override;
        EDF* getEDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
}

#endif /* __SLR_basic_emitter_surface_properties__ */

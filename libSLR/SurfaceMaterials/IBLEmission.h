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

namespace SLR {
    class SLR_API IBLEmission : public EmitterSurfaceProperty {
        const Scene* m_scene;
        const SpectrumTexture* m_coeffM;
        float m_scale;
    public:
        IBLEmission(const Scene* scene, const SpectrumTexture* coeffM, float scale) :
        m_scene(scene), m_coeffM(coeffM), m_scale(scale) {};
        
        SampledSpectrum emittance(const SurfacePoint &surfPt, const WavelengthSamples &wls) const override;
        EDF* getEDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
        
        RegularConstantContinuous2D* createIBLImportanceMap() const;
    };
}

#endif /* defined(__SLR__IBLEmission__) */

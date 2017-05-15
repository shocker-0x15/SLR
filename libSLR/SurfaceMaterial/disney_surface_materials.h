//
//  disney_surface_materials.h
//
//  Created by 渡部 心 on 2017/05/14.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef __SLR_disney_surface_materials__
#define __SLR_disney_surface_materials__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/surface_material.h"

namespace SLR {
    class SLR_API DisneyReflectionSurfaceMaterial : public SurfaceMaterial {
        const SpectrumTexture* m_baseColor;
        const FloatTexture* m_subsurface;
        const FloatTexture* m_metallic;
        const FloatTexture* m_specular;
        const FloatTexture* m_specularTint;
        const FloatTexture* m_roughness;
        const FloatTexture* m_anisotropic;
        const FloatTexture* m_sheen;
        const FloatTexture* m_sheenTint;
        const FloatTexture* m_clearCoat;
        const FloatTexture* m_clearCoatGloss;
    public:
        DisneyReflectionSurfaceMaterial(const SpectrumTexture* baseColor, 
                                        const FloatTexture* subsurface, const FloatTexture* metallic, const FloatTexture* specular, const FloatTexture* specularTint, 
                                        const FloatTexture* roughness, const FloatTexture* anisotropic, const FloatTexture* sheen, const FloatTexture* sheenTint, 
                                        const FloatTexture* clearCoat, const FloatTexture* clearCoatGloss) :
        m_baseColor(baseColor), 
        m_subsurface(subsurface), m_metallic(metallic), m_specular(specular), m_specularTint(specularTint), m_roughness(roughness), m_anisotropic(anisotropic), 
        m_sheen(sheen), m_sheenTint(sheenTint), m_clearCoat(clearCoat), m_clearCoatGloss(clearCoatGloss) { }
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
}

#endif /* __SLR_disney_surface_materials__ */

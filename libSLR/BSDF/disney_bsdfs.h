//
//  disney_bsdfs.h
//
//  Created by 渡部 心 on 2017/05/04.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef __SLR_disney_bsdfs__
#define __SLR_disney_bsdfs__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/directional_distribution_functions.h"

namespace SLR {
    // References
    // Physically-Based Shading at Disney
    // BRDF Explorer
    class SLR_API DisneyBRDF : public BSDF {
        SampledSpectrum m_baseColor;
        float m_subsurface;
        float m_metallic;
        float m_specular;
        float m_specularTint;
        float m_roughness;
        float m_anisotropic;
        float m_sheen;
        float m_sheenTint;
        float m_clearCoat;
        float m_clearCoatGloss;
        
        float m_baseColorLuminance;
        GGXMicrofacetDistribution* m_base_D;
        BerryMicrofacetDistribution* m_clearcoat_D;
        GGXMicrofacetDistribution* m_clearcoat_forG;
        
        SampledSpectrum sampleInternal(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult* result) const override;
        SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const override;
        float evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, float* revPDF) const override;
        float weightInternal(const BSDFQuery &query) const override;
        SampledSpectrum getBaseColorInternal(DirectionType flags) const override;
    public:
        DisneyBRDF(const SampledSpectrum &baseColor, float baseColorLuminance, 
                   float subsurface, float metallic, float specular, float specularTint, float roughness, float anisotropic, 
                   float sheen, float sheenTint, float clearCoat, float clearCoatGloss);
        ~DisneyBRDF();
    };
    
    
    
    
//    // References
//    // Extending the Disney BRDF to a BSDF with Integrated Subsurface Scattering
//    class SLR_API DisneyBSDF : public BSDF {
//    };
}


#endif /* __SLR_disney_bsdfs__ */

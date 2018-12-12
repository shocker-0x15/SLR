//
//  microfacet_bsdfs.h
//
//  Created by 渡部 心 on 2016/05/03.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef __SLR_microfacet_bsdfs__
#define __SLR_microfacet_bsdfs__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/directional_distribution_functions.h"

namespace SLR {
    // References
    // Microfacet Models for Refraction through Rough Surfaces
    
    class SLR_API MicrofacetBRDF : public BSDF {
        FresnelConductor m_F;
        const MicrofacetDistribution* m_D;
        
        SampledSpectrum sampleInternal(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult* result) const override;
        SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const override;
        float evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, float* revPDF) const override;
        float weightInternal(const BSDFQuery &query) const override;
        SampledSpectrum getBaseColorInternal(DirectionType flags) const override;
    public:
        MicrofacetBRDF(const SampledSpectrum &eta, const SampledSpectrum &k, const MicrofacetDistribution* D) :
        BSDF(DirectionType::Reflection | DirectionType::HighFreq), m_F(eta, k), m_D(D) { }
    };
    
    
    
    class SLR_API MicrofacetBSDF : public BSDF {
        FresnelDielectric m_F;
        const MicrofacetDistribution* m_D;
        
        SampledSpectrum sampleInternal(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult* result) const override;
        SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const override;
        float evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, float* revPDF) const override;
        float weightInternal(const BSDFQuery &query) const override;
        SampledSpectrum getBaseColorInternal(DirectionType flags) const override;
    public:
        MicrofacetBSDF(const SampledSpectrum &etaExt, const SampledSpectrum &etaInt, const MicrofacetDistribution* D) :
        BSDF(DirectionType::Reflection | DirectionType::Transmission | DirectionType::HighFreq),
        m_F(etaExt, etaInt), m_D(D) { }
    };
}

#endif /* __SLR_microfacet_bsdfs__ */

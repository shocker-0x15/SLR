//
//  MicrofacetBSDF.h
//  SLR
//
//  Created by 渡部 心 on 2016/05/03.
//  Copyright © 2016年 渡部 心. All rights reserved.
//

#ifndef MicrofacetBSDF_h
#define MicrofacetBSDF_h

#include "../defines.h"
#include "../references.h"
#include "../Core/directional_distribution_functions.h"

namespace SLR {    
    class SLR_API MicrofacetBRDF : public BSDF {
        const Fresnel* m_F;
        const MicrofacetDistribution* m_D;
        
        SampledSpectrum sampleInternal(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const override;
        SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const override;
        float evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, float* revPDF) const override;
        float weightInternal(const BSDFQuery &query) const override;
        SampledSpectrum getBaseColorInternal(DirectionType flags) const override;
    public:
        MicrofacetBRDF(const Fresnel* F, const MicrofacetDistribution* D) :
        BSDF(DirectionType::Reflection | DirectionType::HighFreq), m_F(F), m_D(D) { }
    };
    
    
    
    class SLR_API MicrofacetBTDF : public BSDF {
        const FresnelDielectric m_F;
        const MicrofacetDistribution* m_D;
        
        SampledSpectrum sampleInternal(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const override;
        SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const override;
        float evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, float* revPDF) const override;
        float weightInternal(const BSDFQuery &query) const override;
        SampledSpectrum getBaseColorInternal(DirectionType flags) const override;
    public:
        MicrofacetBTDF(const SampledSpectrum &etaExt, const SampledSpectrum &etaInt, bool dispersive, const MicrofacetDistribution* D) :
        BSDF(DirectionType::Reflection | DirectionType::HighFreq | (dispersive ? DirectionType::Dispersive : DirectionType())),
        m_F(etaExt, etaInt), m_D(D) { }
    };
}

#endif /* MicrofacetBSDF_hpp */

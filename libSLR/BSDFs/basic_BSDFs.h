//
//  basic_BSDFs.h
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__basic_BSDFs__
#define __SLR__basic_BSDFs__

#include "../defines.h"
#include "../references.h"
#include "../Core/directional_distribution_functions.h"

namespace SLR {
    class SLR_API LambertianBRDF : public BSDF {
        SampledSpectrum m_R;
        
        SampledSpectrum sampleInternal(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const override;
        SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const override;
        float evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, float* revPDF) const override;
        float weightInternal(const BSDFQuery &query, const BSDFSample &smp) const override;
        float weightInternal(const BSDFQuery &query, const Vector3D &dir, float* revWeight) const override;
        SampledSpectrum getBaseColorInternal(DirectionType flags) const override;
    public:
        LambertianBRDF(const SampledSpectrum &R) : BSDF(DirectionType::Reflection | DirectionType::LowFreq), m_R(R) { }
    };
    
    class SLR_API SpecularBRDF : public BSDF {
        SampledSpectrum m_coeffR;
        const Fresnel* m_fresnel;
        
        SampledSpectrum sampleInternal(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const override;
        SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const override;
        float evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, float* revPDF) const override;
        float weightInternal(const BSDFQuery &query, const BSDFSample &smp) const override;
        float weightInternal(const BSDFQuery &query, const Vector3D &dir, float* revWeight) const override;
        SampledSpectrum getBaseColorInternal(DirectionType flags) const override;
    public:
        SpecularBRDF(const SampledSpectrum &coeffR, const Fresnel* fresnel) : BSDF(DirectionType::Reflection | DirectionType::Delta0D), m_coeffR(coeffR), m_fresnel(fresnel) { }
    };
    
    class SLR_API SpecularBTDF : public BSDF {
        SampledSpectrum m_coeffT;
        FresnelDielectric m_fresnel;
        
        SampledSpectrum sampleInternal(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const override;
        SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const override;
        float evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, float* revPDF) const override;
        float weightInternal(const BSDFQuery &query, const BSDFSample &smp) const override;
        float weightInternal(const BSDFQuery &query, const Vector3D &dir, float* revWeight) const override;
        SampledSpectrum getBaseColorInternal(DirectionType flags) const override;
    public:
        SpecularBTDF(const SampledSpectrum &coeffT, const SampledSpectrum &etaExt, const SampledSpectrum &etaInt, bool dispersive) :
        BSDF(DirectionType::Transmission | DirectionType::Delta0D | (dispersive ? DirectionType::Dispersive : DirectionType())),
        m_coeffT(coeffT), m_fresnel(etaExt, etaInt) { }
    };
    
    class SLR_API InverseBSDF : public BSDF {
        const BSDF* m_baseBSDF;
        
        SampledSpectrum sampleInternal(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const override;
        SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const override;
        float evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, float* revPDF) const override;
        float weightInternal(const BSDFQuery &query, const BSDFSample &smp) const override;
        float weightInternal(const BSDFQuery &query, const Vector3D &dir, float* revWeight) const override;
        SampledSpectrum getBaseColorInternal(DirectionType flags) const override;
    public:
        InverseBSDF(const BSDF* baseBSDF) : BSDF(baseBSDF->m_type.flip()), m_baseBSDF(baseBSDF) { }
        
        bool matches(DirectionType flags) const override { return m_baseBSDF->matches(flags.flip()); }
    };
}

#endif /* defined(__SLR__basic_BSDFs__) */

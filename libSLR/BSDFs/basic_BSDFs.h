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
        
        SampledSpectrum sampleInternal(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult* result) const override;
        SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const override;
        SampledSpectrum evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* revPDF) const override;
        SampledSpectrum weightInternal(const BSDFQuery &query) const override;
        SampledSpectrum getBaseColorInternal(DirectionType flags) const override;
    public:
        LambertianBRDF(const SampledSpectrum &R) : BSDF(DirectionType::Reflection | DirectionType::LowFreq), m_R(R) { }
    };
    
    
    
    class SLR_API SpecularBRDF : public BSDF {
        SampledSpectrum m_coeffR;
        FresnelConductor m_fresnel;
        
        SampledSpectrum sampleInternal(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult* result) const override;
        SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const override;
        SampledSpectrum evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* revPDF) const override;
        SampledSpectrum weightInternal(const BSDFQuery &query) const override;
        SampledSpectrum getBaseColorInternal(DirectionType flags) const override;
    public:
        SpecularBRDF(const SampledSpectrum &coeffR, const SampledSpectrum &eta, const SampledSpectrum &k) : BSDF(DirectionType::Reflection | DirectionType::Delta0D), m_coeffR(coeffR), m_fresnel(eta, k) { }
    };
    
    
    
    class SLR_API SpecularBSDF : public BSDF {
        SampledSpectrum m_coeff;
        FresnelDielectric m_fresnel;
        
        SampledSpectrum sampleInternal(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult* result) const override;
        SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const override;
        SampledSpectrum evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* revPDF) const override;
        SampledSpectrum weightInternal(const BSDFQuery &query) const override;
        SampledSpectrum getBaseColorInternal(DirectionType flags) const override;
    public:
        SpecularBSDF(const SampledSpectrum &coeff, const SampledSpectrum &etaExt, const SampledSpectrum &etaInt) :
        BSDF(DirectionType::Transmission | DirectionType::Delta0D | DirectionType::Dispersive),
        m_coeff(coeff), m_fresnel(etaExt, etaInt) { }
    };
    
    
    
    class SLR_API InverseBSDF : public BSDF {
        const BSDF* m_baseBSDF;
        
        SampledSpectrum sampleInternal(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult* result) const override;
        SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const override;
        SampledSpectrum evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* revPDF) const override;
        SampledSpectrum weightInternal(const BSDFQuery &query) const override;
        SampledSpectrum getBaseColorInternal(DirectionType flags) const override;
    public:
        InverseBSDF(const BSDF* baseBSDF) : BSDF(baseBSDF->m_type.flip()), m_baseBSDF(baseBSDF) { }
    };
    
    
    
    class SLR_API NullBSDF : public BSDF {
        SampledSpectrum sampleInternal(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult* result) const override {
            SLRAssert(false, "NullBSDF's method should not be called.");
            return SampledSpectrum::Zero;
        };
        SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const override {
            SLRAssert(false, "NullBSDF's method should not be called.");
            return SampledSpectrum::Zero;
        };
        SampledSpectrum evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* revPDF) const override {
            SLRAssert(false, "NullBSDF's method should not be called.");
            return SampledSpectrum::Zero;
        };
        SampledSpectrum weightInternal(const BSDFQuery &query) const override {
            SLRAssert(false, "NullBSDF's method should not be called.");
            return SampledSpectrum::Zero;
        };
        SampledSpectrum getBaseColorInternal(DirectionType flags) const override {
            SLRAssert(false, "NullBSDF's method should not be called.");
            return SampledSpectrum::Zero;
        };
    public:
        NullBSDF() : BSDF(DirectionType()) { }
    };
}

#endif /* defined(__SLR__basic_BSDFs__) */

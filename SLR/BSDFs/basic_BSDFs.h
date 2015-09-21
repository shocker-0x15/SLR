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

class LambertianBRDF : public BSDF {
    SampledSpectrum m_R;
    
    SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir) const override;
public:
    LambertianBRDF(const SampledSpectrum &R) : BSDF(DirectionType::Reflection | DirectionType::LowFreq), m_R(R) { };
    
    SampledSpectrum sample(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const override;
    float evaluatePDF(const BSDFQuery &query, const Vector3D &dir) const override;
    float weight(const BSDFQuery &query, const BSDFSample &smp) const override;
    float weight(const BSDFQuery &query, const Vector3D &dir) const override;
    
    SampledSpectrum getBaseColor(DirectionType flags) const override;
};

class SpecularBRDF : public BSDF {
    SampledSpectrum m_coeffR;
    const Fresnel* m_fresnel;
    
    SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir) const override;
public:
    SpecularBRDF(const SampledSpectrum &coeffR, const Fresnel* fresnel) :
    BSDF(DirectionType::Reflection | DirectionType::Delta0D), m_coeffR(coeffR), m_fresnel(fresnel) { };
    
    SampledSpectrum sample(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const override;
    float evaluatePDF(const BSDFQuery &query, const Vector3D &dir) const override;
    float weight(const BSDFQuery &query, const BSDFSample &smp) const override;
    float weight(const BSDFQuery &query, const Vector3D &dir) const override;
    
    SampledSpectrum getBaseColor(DirectionType flags) const override;
};

class SpecularBTDF : public BSDF {
    SampledSpectrum m_coeffT;
    FresnelDielectric m_fresnel;
    bool m_dispersive;
    
    SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir) const override;
public:
    SpecularBTDF(const SampledSpectrum &coeffT, const SampledSpectrum &etaExt, const SampledSpectrum &etaInt, bool dispersive) :
    BSDF(DirectionType::Transmission | DirectionType::Delta0D), m_coeffT(coeffT), m_fresnel(etaExt, etaInt), m_dispersive(dispersive) { };
    
    SampledSpectrum sample(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const override;
    float evaluatePDF(const BSDFQuery &query, const Vector3D &dir) const override;
    float weight(const BSDFQuery &query, const BSDFSample &smp) const override;
    float weight(const BSDFQuery &query, const Vector3D &dir) const override;
    
    SampledSpectrum getBaseColor(DirectionType flags) const override;
};

#endif /* defined(__SLR__basic_BSDFs__) */

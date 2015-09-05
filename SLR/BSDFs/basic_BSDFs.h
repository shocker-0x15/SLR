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
    Spectrum m_R;
    
    Spectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir) const override;
public:
    LambertianBRDF(const Spectrum &R) : BSDF(DirectionType::Reflection | DirectionType::LowFreq), m_R(R) { };
    
    Spectrum sample(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const override;
    float evaluatePDF(const BSDFQuery &query, const Vector3D &dir) const override;
    float weight(const BSDFQuery &query) const override;
    
    Spectrum getBaseColor(DirectionType flags) const override;
};

class SpecularBRDF : public BSDF {
    Spectrum m_coeffR;
    const Fresnel* m_fresnel;
    
    Spectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir) const override;
public:
    SpecularBRDF(const Spectrum &coeffR, const Fresnel* fresnel) :
    BSDF(DirectionType::Reflection | DirectionType::Delta0D), m_coeffR(coeffR), m_fresnel(fresnel) { };
    
    Spectrum sample(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const override;
    float evaluatePDF(const BSDFQuery &query, const Vector3D &dir) const override;
    float weight(const BSDFQuery &query) const override;
    
    Spectrum getBaseColor(DirectionType flags) const override;
};

class SpecularBTDF : public BSDF {
    Spectrum m_coeffT;
    FresnelDielectric m_fresnel;
    
    Spectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir) const override;
public:
    SpecularBTDF(const Spectrum &coeffT, float etaExt, float etaInt) :
    BSDF(DirectionType::Transmission | DirectionType::Delta0D), m_coeffT(coeffT), m_fresnel(etaExt, etaInt) { };
    
    Spectrum sample(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const override;
    float evaluatePDF(const BSDFQuery &query, const Vector3D &dir) const override;
    float weight(const BSDFQuery &query) const override;
    
    Spectrum getBaseColor(DirectionType flags) const override;
};

#endif /* defined(__SLR__basic_BSDFs__) */

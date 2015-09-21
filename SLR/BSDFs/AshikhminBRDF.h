//
//  AshikhminBRDF.h
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__AshikhminBRDF__
#define __SLR__AshikhminBRDF__

#include "../defines.h"
#include "../references.h"
#include "../Core/directional_distribution_functions.h"

class AshikhminSpecularBRDF : public BSDF {
    SampledSpectrum m_Rs;
    float m_nu, m_nv;
    
    SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir) const override;
public:
    AshikhminSpecularBRDF(const SampledSpectrum &Rs, float nu, float nv) :
    BSDF(DirectionType::Reflection | DirectionType::HighFreq), m_Rs(Rs), m_nu(nu), m_nv(nv) { };
    
    SampledSpectrum sample(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const override;
    float evaluatePDF(const BSDFQuery &query, const Vector3D &dir) const override;
    float weight(const BSDFQuery &query, const BSDFSample &smp) const override;
    float weight(const BSDFQuery &query, const Vector3D &dir) const override;
    
    SampledSpectrum getBaseColor(DirectionType flags) const override;
};

class AshikhminDiffuseBRDF : public BSDF {
    SampledSpectrum m_Rs;
    SampledSpectrum m_Rd;
    
    SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir) const override;
public:
    AshikhminDiffuseBRDF(const SampledSpectrum &Rs, const SampledSpectrum &Rd) :
    BSDF(DirectionType::Reflection | DirectionType::LowFreq), m_Rs(Rs), m_Rd(Rd) { };
    
    SampledSpectrum sample(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const override;
    float evaluatePDF(const BSDFQuery &query, const Vector3D &dir) const override;
    float weight(const BSDFQuery &query, const BSDFSample &smp) const override;
    float weight(const BSDFQuery &query, const Vector3D &dir) const override;
    
    SampledSpectrum getBaseColor(DirectionType flags) const override;
};

#endif /* defined(__SLR__AshikhminBRDF__) */

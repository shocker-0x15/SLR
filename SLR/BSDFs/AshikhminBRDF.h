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
    Spectrum m_Rs;
    float m_nu, m_nv;
    
    Spectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir) const override;
public:
    AshikhminSpecularBRDF(const Spectrum &Rs, float nu, float nv) :
    BSDF(DirectionType::Reflection | DirectionType::HighFreq), m_Rs(Rs), m_nu(nu), m_nv(nv) { };
    
    Spectrum sample(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const override;
    float evaluatePDF(const BSDFQuery &query, const Vector3D &dir) const override;
    float weight(const BSDFQuery &query) const override;
    
    Spectrum getBaseColor(DirectionType flags) const override;
};

class AshikhminDiffuseBRDF : public BSDF {
    Spectrum m_Rs;
    Spectrum m_Rd;
    
    Spectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir) const override;
public:
    AshikhminDiffuseBRDF(const Spectrum &Rs, const Spectrum &Rd) :
    BSDF(DirectionType::Reflection | DirectionType::LowFreq), m_Rs(Rs), m_Rd(Rd) { };
    
    Spectrum sample(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const override;
    float evaluatePDF(const BSDFQuery &query, const Vector3D &dir) const override;
    float weight(const BSDFQuery &query) const override;
    
    Spectrum getBaseColor(DirectionType flags) const override;
};

#endif /* defined(__SLR__AshikhminBRDF__) */

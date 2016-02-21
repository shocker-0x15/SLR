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

namespace SLR {
    class SLR_API AshikhminSpecularBRDF : public BSDF {
        SampledSpectrum m_Rs;
        float m_nu, m_nv;
        
        SampledSpectrum sampleInternal(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const override;
        SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const override;
        float evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, float* revPDF) const override;
        float weightInternal(const BSDFQuery &query) const override;
        SampledSpectrum getBaseColorInternal(DirectionType flags) const override;
    public:
        AshikhminSpecularBRDF(const SampledSpectrum &Rs, float nu, float nv) : BSDF(DirectionType::Reflection | DirectionType::HighFreq),
        m_Rs(Rs), m_nu(nu), m_nv(nv) { }
    };
    
    class SLR_API AshikhminDiffuseBRDF : public BSDF {
        SampledSpectrum m_Rs;
        SampledSpectrum m_Rd;
        
        SampledSpectrum sampleInternal(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const override;
        SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const override;
        float evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, float* revPDF) const override;
        float weightInternal(const BSDFQuery &query) const override;
        SampledSpectrum getBaseColorInternal(DirectionType flags) const override;
    public:
        AshikhminDiffuseBRDF(const SampledSpectrum &Rs, const SampledSpectrum &Rd) : BSDF(DirectionType::Reflection | DirectionType::LowFreq), 
        m_Rs(Rs), m_Rd(Rd) { }
    };    
}

#endif /* defined(__SLR__AshikhminBRDF__) */

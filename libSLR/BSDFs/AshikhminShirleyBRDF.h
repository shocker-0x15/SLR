//
//  AshikhminShirleyBRDF.h
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

// References
// An Anisotropic Phong BRDF Model

#ifndef __SLR__AshikhminShirleyBRDF__
#define __SLR__AshikhminShirleyBRDF__

#include "../defines.h"
#include "../references.h"
#include "../Core/directional_distribution_functions.h"

namespace SLR {
    class SLR_API AshikhminShirleyBRDF : public BSDF {
        SampledSpectrum m_Rs;
        SampledSpectrum m_Rd;
        float m_nu, m_nv;
        
        SampledSpectrum sampleInternal(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult* result) const override;
        SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const override;
        SampledSpectrum evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* revPDF) const override;
        SampledSpectrum weightInternal(const BSDFQuery &query) const override;
        SampledSpectrum getBaseColorInternal(DirectionType flags) const override;
    public:
        AshikhminShirleyBRDF(const SampledSpectrum &Rs, const SampledSpectrum &Rd, float nu, float nv) : BSDF(DirectionType::Reflection | DirectionType::HighFreq | DirectionType::LowFreq),
        m_Rs(Rs), m_Rd(Rd), m_nu(nu), m_nv(nv) { }
    };    
}

#endif /* defined(__SLR__AshikhminBRDF__) */

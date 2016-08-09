//
//  OrenNayerBRDF.h
//
//  Created by 渡部 心 on 2016/06/17.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef OrenNayerBRDF_h
#define OrenNayerBRDF_h

#include "../defines.h"
#include "../references.h"
#include "../Core/directional_distribution_functions.h"

namespace SLR {
    class SLR_API OrenNayerBRDF : public BSDF {
        SampledSpectrum m_R;
        float m_A, m_B;
        
        SampledSpectrum sampleInternal(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult* result) const override;
        SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const override;
        float evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, float* revPDF) const override;
        float weightInternal(const BSDFQuery &query) const override;
        SampledSpectrum getBaseColorInternal(DirectionType flags) const override;
    public:
        OrenNayerBRDF(const SampledSpectrum &R, float sigma) :
        BSDF(DirectionType::Reflection | DirectionType::LowFreq), m_R(R),
        m_A(1.0f - 0.5f * sigma * sigma / (sigma * sigma + 0.33)), m_B(0.45 * sigma * sigma / (sigma * sigma + 0.09)) { }
    };
}

#endif /* OrenNayerBRDF_h */

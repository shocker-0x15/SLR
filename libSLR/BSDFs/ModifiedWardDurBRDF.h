//
//  ModifiedWardDurBRDF.h
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__ModifiedWardDurBRDF__
#define __SLR__ModifiedWardDurBRDF__

#include "../defines.h"
#include "../references.h"
#include "../Core/directional_distribution_functions.h"

namespace SLR {
    class SLR_API ModifiedWardDurBRDF : public BSDF {
        SampledSpectrum m_R;
        float m_anisoX, m_anisoY;
        
        SampledSpectrum sampleInternal(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const override;
        SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const override;
        SampledSpectrum evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* revPDF) const override;
        SampledSpectrum weightInternal(const BSDFQuery &query) const override;
        SampledSpectrum getBaseColorInternal(DirectionType flags) const override;
    public:
        ModifiedWardDurBRDF(const SampledSpectrum &R, float ax, float ay) : BSDF(DirectionType::Reflection | DirectionType::HighFreq),
        m_R(R), m_anisoX(ax), m_anisoY(ay) { }
    };    
}

#endif /* defined(__SLR__ModifiedWardDurBRDF__) */

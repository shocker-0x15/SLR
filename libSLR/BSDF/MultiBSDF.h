//
//  MultiBSDF.h
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_MultiBSDF__
#define __SLR_MultiBSDF__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/directional_distribution_functions.h"

namespace SLR {
    class SLR_API MultiBSDF : public BSDF {
        static const uint32_t MaxNumElems = 4;
        uint32_t m_numComponents;
        const BSDF* m_BSDFs[MaxNumElems];
        
        SampledSpectrum sampleInternalNoRev(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult *result) const;
        SampledSpectrum sampleInternalWithRev(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult *result) const;
        float evaluatePDFInternalNoRev(const BSDFQuery &query, const Vector3D &dirOut, float* revPDF) const;
        float evaluatePDFInternalWithRev(const BSDFQuery &query, const Vector3D &dirOut, float* revPDF) const;
        
        SampledSpectrum sampleInternal(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult* result) const override {
            return query.requestReverse ? sampleInternalWithRev(query, uComponent, uDir, result) : sampleInternalNoRev(query, uComponent, uDir, result);
        }
        SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dirOut, SampledSpectrum* rev_fs) const override;
        float evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dirOut, float* revPDF) const override {
            return revPDF ? evaluatePDFInternalWithRev(query, dirOut, revPDF) : evaluatePDFInternalNoRev(query, dirOut, revPDF);
        }
        float weightInternal(const BSDFQuery &query) const override;
        SampledSpectrum getBaseColorInternal(DirectionType flags) const override;
    public:
        MultiBSDF() : BSDF(DirectionType()), m_numComponents(0) { }
        
        void add(const BSDF* bsdf);
    };
}

#endif /* __SLR_MultiBSDF__ */

//
//  MultiBSDF.h
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__MultiBSDF__
#define __SLR__MultiBSDF__

#include "../defines.h"
#include "../references.h"
#include "../Core/directional_distribution_functions.h"

namespace SLR {
    class SLR_API MultiBSDF : public BSDF {
#define maxNumElems 4
        uint32_t m_numComponents;
        const BSDF* m_BSDFs[maxNumElems];
        
        SampledSpectrum sampleInternalNoRev(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult *result) const;
        SampledSpectrum sampleInternalWithRev(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult *result) const;
        SampledSpectrum evaluatePDFInternalNoRev(const BSDFQuery &query, const Vector3D &dirOut, SampledSpectrum* revPDF) const;
        SampledSpectrum evaluatePDFInternalWithRev(const BSDFQuery &query, const Vector3D &dirOut, SampledSpectrum* revPDF) const;
        
        SampledSpectrum sampleInternal(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const override {
            return result->reverse ? sampleInternalWithRev(query, smp, result) : sampleInternalNoRev(query, smp, result);
        }
        SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dirOut, SampledSpectrum* rev_fs) const override;
        SampledSpectrum evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dirOut, SampledSpectrum* revPDF) const override {
            return revPDF ? evaluatePDFInternalWithRev(query, dirOut, revPDF) : evaluatePDFInternalNoRev(query, dirOut, revPDF);
        }
        SampledSpectrum weightInternal(const BSDFQuery &query) const override;
        SampledSpectrum getBaseColorInternal(DirectionType flags) const override;
    public:
        MultiBSDF() : BSDF(DirectionType()), m_numComponents(0) { }
        
        void add(const BSDF* bsdf);
        
        bool matches(DirectionType flags) const override {
            for (int i = 0; i < m_numComponents; ++i)
                if (m_BSDFs[i]->matches(flags))
                    return true;
            return false;
        }
    };    
}

#endif /* defined(__SLR__MultiBSDF__) */

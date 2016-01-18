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
    class MultiBSDF : public BSDF {
#define maxNumElems 4
        uint32_t m_numComponents;
        const BSDF* m_BSDFs[maxNumElems];
        
        SampledSpectrum sampleInternal(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const override;
        SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dirOut) const override;
        float evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dirOut) const override;
        float weightInternal(const BSDFQuery &query, const BSDFSample &smp) const override;
        float weightInternal(const BSDFQuery &query, const Vector3D &dir) const override;
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

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

class MultiBSDF : public BSDF {
#define maxNumElems 4
    uint32_t m_numComponents;
    BSDF* m_BSDFs[maxNumElems];
    
    Spectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dirOut) const override;
public:
    MultiBSDF() : BSDF(DirectionType()), m_numComponents(0) { }
    
    void add(BSDF* bsdf);
    Spectrum sample(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const override;
    float evaluatePDF(const BSDFQuery &query, const Vector3D &dirOut) const override;
    float weight(const BSDFQuery &query) const override;
    
    Spectrum getBaseColor(DirectionType flags) const override;
    
    bool matches(DirectionType flags) const {
        for (int i = 0; i < m_numComponents; ++i)
            if (m_BSDFs[i]->matches(flags))
                return true;
        return false;
    };
};

#endif /* defined(__SLR__MultiBSDF__) */

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

class ModifiedWardDurBRDF : public BSDF {
    Spectrum m_R;
    float m_anisoX, m_anisoY;
    
    Spectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir) const override;
public:
    ModifiedWardDurBRDF(const Spectrum &R, float ax, float ay) :
    BSDF(DirectionType::Reflection | DirectionType::HighFreq), m_R(R), m_anisoX(ax), m_anisoY(ay) { };
    
    Spectrum sample(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const override;
    float evaluatePDF(const BSDFQuery &query, const Vector3D &dir) const override;
    float weight(const BSDFQuery &query) const override;
    
    Spectrum getBaseColor(DirectionType flags) const override;
};

#endif /* defined(__SLR__ModifiedWardDurBRDF__) */

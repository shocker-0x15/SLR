//
//  MultiBSDF.cpp
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "MultiBSDF.h"
#include "../Core/distributions.h"

void MultiBSDF::add(BSDF *bsdf) {
    if (!bsdf)
        return;
    SLRAssert(m_numComponents < maxNumElems, "Number of MultiBSDF elements exceeds the limit: %u", maxNumElems);
    m_BSDFs[m_numComponents++] = bsdf;
    m_type |= bsdf->m_type;
}

#define ACTUAL_WEIGHTS
Spectrum MultiBSDF::sample(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult *result) const {
#ifdef ACTUAL_WEIGHTS
    float weights[maxNumElems];
    for (int i = 0; i < m_numComponents; ++i)
        weights[i] = m_BSDFs[i]->sample(query, smp, result).maxValue() * std::fabs(result->dir_sn.z) / result->dirPDF;
#else
    float weights[maxNumElems];
    for (int i = 0; i < m_numComponents; ++i)
        weights[i] = m_BSDFs[i]->weight(query);
#endif
    
    float sumWeights, base;
    uint32_t idx = sampleDiscrete(weights, &sumWeights, &base, m_numComponents, smp.uComponent);
    BSDF* selectedBSDF = m_BSDFs[idx];
    
    Spectrum value = selectedBSDF->sample(query, smp, result);
    result->dirPDF *= weights[idx];
    
    if (!result->dirType.hasDelta()) {
        BSDFQuery mQuery = query;
        mQuery.flags &= sideTest(query.gNormal_sn, query.dir_sn, result->dir_sn);
        value = Spectrum::Zero;
        for (int i = 0; i < m_numComponents; ++i) {
            value += m_BSDFs[i]->evaluateInternal(mQuery, result->dir_sn);
            if (i != idx)
                result->dirPDF += m_BSDFs[i]->evaluatePDF(query, result->dir_sn) * weights[i];
        }
    }
    result->dirPDF /= sumWeights;
    
    return value;
}

Spectrum MultiBSDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dirOut) const {
    Spectrum retValue;
    for (int i = 0; i < m_numComponents; ++i)
        retValue += m_BSDFs[i]->evaluateInternal(query, dirOut);
    return retValue;
}

float MultiBSDF::evaluatePDF(const BSDFQuery &query, const Vector3D &dirOut) const {
    float sumWeights = 0.0f;
#ifdef ACTUAL_WEIGHTS
    float weights[maxNumElems];
    for (int i = 0; i < m_numComponents; ++i) {
        float PDF = m_BSDFs[i]->evaluatePDF(query, dirOut);
        if (PDF > 0)
            weights[i] = m_BSDFs[i]->evaluate(query, dirOut).maxValue() * std::fabs(dirOut.z) / PDF;
        else
            weights[i] = 0;
        sumWeights += weights[i];
    }
    if (sumWeights == 0.0f)
        return 0.0f;
#else
    float weights[maxNumElems];
    for (int i = 0; i < m_numComponents; ++i) {
        weights[i] = m_BSDFs[i]->weight(query);
        sumWeights += weights[i];
    }
#endif
    
    float retPDF = 0.0f;
    for (int i = 0; i < m_numComponents; ++i) {
        retPDF += m_BSDFs[i]->evaluatePDF(query, dirOut) * weights[i];
    }
    retPDF /= sumWeights;
    
    return retPDF;
}

float MultiBSDF::weight(const BSDFQuery &query) const {
    float sumWeights = 0.0f;
    for (int i = 0; i < m_numComponents; ++i)
        sumWeights += m_BSDFs[i]->weight(query);
    return sumWeights;
}

Spectrum MultiBSDF::getBaseColor(DirectionType flags) const {
    Spectrum baseColor = Spectrum::Zero;
    for (int i = 0; i < m_numComponents; ++i) {
        baseColor = m_BSDFs[i]->getBaseColor(flags);
        if (baseColor.hasNonZero())
            break;
    }
    return baseColor;
}

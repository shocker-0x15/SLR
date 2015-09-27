//
//  MultiBSDF.cpp
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "MultiBSDF.h"
#include "../Core/distributions.h"

namespace SLR {
    void MultiBSDF::add(BSDF *bsdf) {
        if (!bsdf)
            return;
        SLRAssert(m_numComponents < maxNumElems, "Number of MultiBSDF elements exceeds the limit: %u", maxNumElems);
        m_BSDFs[m_numComponents++] = bsdf;
        m_type |= bsdf->m_type;
    }
    
    SampledSpectrum MultiBSDF::sample(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult *result) const {
        float weights[maxNumElems];
        for (int i = 0; i < m_numComponents; ++i)
            weights[i] = m_BSDFs[i]->weight(query, smp);
        
        float sumWeights, base;
        uint32_t idx = sampleDiscrete(weights, &sumWeights, &base, m_numComponents, smp.uComponent);
        BSDF* selectedBSDF = m_BSDFs[idx];
        if (sumWeights == 0.0f) {
            result->dirPDF = 0.0f;
            return SampledSpectrum::Zero;
        }
        
        SampledSpectrum value = selectedBSDF->sample(query, smp, result);
        result->dirPDF *= weights[idx];
        
        if (!result->dirType.hasDelta()) {
            BSDFQuery mQuery = query;
            mQuery.flags &= sideTest(query.gNormal_sn, query.dir_sn, result->dir_sn);
            value = SampledSpectrum::Zero;
            for (int i = 0; i < m_numComponents; ++i) {
                value += m_BSDFs[i]->evaluateInternal(mQuery, result->dir_sn);
                if (i != idx)
                    result->dirPDF += m_BSDFs[i]->evaluatePDF(query, result->dir_sn) * weights[i];
            }
        }
        result->dirPDF /= sumWeights;
        
        return value;
    }
    
    SampledSpectrum MultiBSDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dirOut) const {
        SampledSpectrum retValue;
        for (int i = 0; i < m_numComponents; ++i)
            retValue += m_BSDFs[i]->evaluateInternal(query, dirOut);
        return retValue;
    }
    
    float MultiBSDF::evaluatePDF(const BSDFQuery &query, const Vector3D &dirOut) const {
        float sumWeights = 0.0f;
        float weights[maxNumElems];
        for (int i = 0; i < m_numComponents; ++i) {
            weights[i] = m_BSDFs[i]->weight(query, dirOut);
            sumWeights += weights[i];
        }
        if (sumWeights == 0.0f)
            return 0.0f;
        
        float retPDF = 0.0f;
        for (int i = 0; i < m_numComponents; ++i) {
            retPDF += m_BSDFs[i]->evaluatePDF(query, dirOut) * weights[i];
        }
        retPDF /= sumWeights;
        
        return retPDF;
    }
    
    float MultiBSDF::weight(const BSDFQuery &query, const BSDFSample &smp) const {
        float sumWeights = 0.0f;
        for (int i = 0; i < m_numComponents; ++i)
            sumWeights += m_BSDFs[i]->weight(query, smp);
        return sumWeights;
    }
    
    float MultiBSDF::weight(const BSDFQuery &query, const Vector3D &dir) const {
        float sumWeights = 0.0f;
        for (int i = 0; i < m_numComponents; ++i)
            sumWeights += m_BSDFs[i]->weight(query, dir);
        return sumWeights;
    }
    
    SampledSpectrum MultiBSDF::getBaseColor(DirectionType flags) const {
        SampledSpectrum baseColor = SampledSpectrum::Zero;
        for (int i = 0; i < m_numComponents; ++i) {
            baseColor = m_BSDFs[i]->getBaseColor(flags);
            if (baseColor.hasNonZero())
                break;
        }
        return baseColor;
    }    
}

//
//  MultiBSDF.cpp
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "MultiBSDF.h"
#include "../Core/distributions.h"

namespace SLR {
    void MultiBSDF::add(const BSDF *bsdf) {
        if (!bsdf)
            return;
        SLRAssert(m_numComponents < maxNumElems, "Number of MultiBSDF elements exceeds the limit: %u", maxNumElems);
        m_type |= bsdf->m_type;
        m_BSDFs[m_numComponents++] = bsdf;
    }
    
    SampledSpectrum MultiBSDF::sampleInternalNoRev(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult *result) const {
        SampledSpectrum weights[maxNumElems];
        for (int i = 0; i < m_numComponents; ++i)
            weights[i] = m_BSDFs[i]->weight(query);
        
        SampledSpectrum sumWeights;
        float base;
        uint32_t idx = sampleDiscrete(weights, query.heroIndex, &sumWeights, &base, m_numComponents, uComponent);
        const BSDF* selectedBSDF = m_BSDFs[idx];
        if (sumWeights[query.heroIndex] == 0.0f) {
            result->dirPDF = 0.0f;
            return SampledSpectrum::Zero;
        }
        
        uComponent = (uComponent * sumWeights[query.heroIndex] - base) / weights[idx][query.heroIndex];
        SampledSpectrum value = selectedBSDF->sampleInternal(query, uComponent, uDir, result);
        result->dirPDF *= weights[idx];
        if (result->dirPDF == 0.0f) {
            result->dirPDF = 0.0f;
            return SampledSpectrum::Zero;
        }
        
        if (!result->dirType.isDelta()) {
            for (int i = 0; i < m_numComponents; ++i) {
                if (i != idx && m_BSDFs[i]->matches(query.flags))
                    result->dirPDF += m_BSDFs[i]->evaluatePDFInternal(query, result->dir_sn, nullptr) * weights[i];
            }
            
            BSDFQuery mQuery = query;
            mQuery.flags &= sideTest(query.gNormal_sn, query.dir_sn, result->dir_sn);
            value = SampledSpectrum::Zero;
            for (int i = 0; i < m_numComponents; ++i) {
                if (!m_BSDFs[i]->matches(mQuery.flags))
                    continue;
                value += m_BSDFs[i]->evaluateInternal(mQuery, result->dir_sn, nullptr);
            }
        }
        result->dirPDF /= sumWeights;
        
        return value;
    }
    
    SampledSpectrum MultiBSDF::sampleInternalWithRev(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult *result) const {
        SampledSpectrum weights[maxNumElems];
        for (int i = 0; i < m_numComponents; ++i)
            weights[i] = m_BSDFs[i]->weight(query);
        
        SampledSpectrum sumWeights;
        float base;
        uint32_t idx = sampleDiscrete(weights, query.heroIndex, &sumWeights, &base, m_numComponents, uComponent);
        const BSDF* selectedBSDF = m_BSDFs[idx];
        if (sumWeights[query.heroIndex] == 0.0f) {
            result->dirPDF = 0.0f;
            return SampledSpectrum::Zero;
        }
        
        uComponent = (uComponent * sumWeights[query.heroIndex] - base) / weights[idx][query.heroIndex];
        SampledSpectrum value = selectedBSDF->sampleInternal(query, uComponent, uDir, result);
        if (result->dirPDF == 0.0f) {
            result->dirPDF = 0.0f;
            return SampledSpectrum::Zero;
        }
        
        BSDFQuery revQuery = query;// mQuery?
        Vector3D revDirIn = result->dir_sn;
        std::swap(revQuery.dir_sn, revDirIn);
        revQuery.adjoint ^= true;
        SampledSpectrum revWeights[maxNumElems];
        SampledSpectrumSum sumRevWeights = SampledSpectrum::Zero;
        for (int i = 0; i < m_numComponents; ++i) {
            revWeights[i] = m_BSDFs[i]->weight(revQuery);
            sumRevWeights += revWeights[i];
        }
        
        result->dirPDF *= weights[idx];
        result->reverse->dirPDF *= revWeights[idx];
        
        if (!result->dirType.isDelta()) {
            for (int i = 0; i < m_numComponents; ++i) {
                SampledSpectrum revPDF;
                if (i != idx && m_BSDFs[i]->matches(query.flags)) {
                    result->dirPDF += m_BSDFs[i]->evaluatePDFInternal(query, result->dir_sn, &revPDF) * weights[i];
                    result->reverse->dirPDF += revPDF * revWeights[i];
                }
            }
            
            BSDFQuery mQuery = query;
            mQuery.flags &= sideTest(query.gNormal_sn, query.dir_sn, result->dir_sn);
            value = SampledSpectrum::Zero;
            result->reverse->fs = SampledSpectrum::Zero;
            
            for (int i = 0; i < m_numComponents; ++i) {
                if (!m_BSDFs[i]->matches(mQuery.flags))
                    continue;
                SampledSpectrum eRev_fs;
                value += m_BSDFs[i]->evaluateInternal(mQuery, result->dir_sn, &eRev_fs);
                result->reverse->fs += eRev_fs;
            }
        }
        result->dirPDF /= sumWeights;
        result->reverse->dirPDF /= sumRevWeights;
        
        return value;
    }
    
    SampledSpectrum MultiBSDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dirOut, SampledSpectrum* rev_fs) const {
        if (rev_fs) {
            SampledSpectrum retValue = SampledSpectrum::Zero;
            *rev_fs = SampledSpectrum::Zero;
            for (int i = 0; i < m_numComponents; ++i) {
                if (!m_BSDFs[i]->matches(query.flags))
                    continue;
                SampledSpectrum eRev_fs;
                retValue += m_BSDFs[i]->evaluateInternal(query, dirOut, &eRev_fs);
                *rev_fs += eRev_fs;
            }
            return retValue;
        }
        else {
            SampledSpectrum retValue = SampledSpectrum::Zero;
            for (int i = 0; i < m_numComponents; ++i) {
                if (!m_BSDFs[i]->matches(query.flags))
                    continue;
                retValue += m_BSDFs[i]->evaluateInternal(query, dirOut, nullptr);
            }
            return retValue;
        }
    }
    
    SampledSpectrum MultiBSDF::evaluatePDFInternalNoRev(const BSDFQuery &query, const Vector3D &dirOut, SampledSpectrum* revPDF) const {
        SampledSpectrumSum sumWeights = SampledSpectrum::Zero;
        SampledSpectrum weights[maxNumElems];
        for (int i = 0; i < m_numComponents; ++i) {
            weights[i] = m_BSDFs[i]->weight(query);
            sumWeights += weights[i];
        }
        if (sumWeights.result[query.heroIndex] == 0.0f)
            return 0.0f;
        
        SampledSpectrum retPDF = 0.0f;
        for (int i = 0; i < m_numComponents; ++i) {
            retPDF += m_BSDFs[i]->evaluatePDFInternal(query, dirOut, nullptr) * weights[i];
        }
        retPDF /= sumWeights;
        
        return retPDF;
    }
    
    SampledSpectrum MultiBSDF::evaluatePDFInternalWithRev(const BSDFQuery &query, const Vector3D &dirOut, SampledSpectrum* revPDF) const {
        BSDFQuery revQuery = query;// mQuery?
        Vector3D revDirIn = dirOut;
        std::swap(revQuery.dir_sn, revDirIn);
        revQuery.adjoint ^= true;
        
        SampledSpectrumSum sumWeights = SampledSpectrum::Zero;
        SampledSpectrumSum sumRevWeights = SampledSpectrum::Zero;
        SampledSpectrum weights[maxNumElems];
        SampledSpectrum revWeights[maxNumElems];
        for (int i = 0; i < m_numComponents; ++i) {
            weights[i] = m_BSDFs[i]->weight(query);
            sumWeights += weights[i];
            revWeights[i] = m_BSDFs[i]->weight(revQuery);
            sumRevWeights += revWeights[i];
        }
        if (sumWeights.result[query.heroIndex] == 0.0f) {
            *revPDF = 0.0f;
            return 0.0f;
        }
        
        SampledSpectrum retPDF = 0.0f;
        *revPDF = 0.0f;
        for (int i = 0; i < m_numComponents; ++i) {
            SampledSpectrum eRevPDF;
            retPDF += m_BSDFs[i]->evaluatePDFInternal(query, dirOut, &eRevPDF) * weights[i];
            *revPDF += eRevPDF * revWeights[i];
        }
        retPDF /= sumWeights;
        *revPDF /= sumRevWeights;
        
        return retPDF;
    }
    
    SampledSpectrum MultiBSDF::weightInternal(const SLR::BSDFQuery &query) const {
        SampledSpectrumSum sumWeights = SampledSpectrum::Zero;
        for (int i = 0; i < m_numComponents; ++i)
            sumWeights += m_BSDFs[i]->weight(query);
        return sumWeights;
    }
    
    SampledSpectrum MultiBSDF::getBaseColorInternal(DirectionType flags) const {
        SampledSpectrum baseColor = SampledSpectrum::Zero;
        for (int i = 0; i < m_numComponents; ++i) {
            baseColor = m_BSDFs[i]->getBaseColor(flags);
            if (baseColor.hasNonZero())
                break;
        }
        return baseColor;
    }
}

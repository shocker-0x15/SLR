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
        m_BSDFs[m_numComponents++] = bsdf;
    }
    
    SampledSpectrum MultiBSDF::sampleInternalNoRev(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult *result) const {
        float weights[maxNumElems];
        for (int i = 0; i < m_numComponents; ++i)
            weights[i] = m_BSDFs[i]->weight(query, smp);
        
        float sumWeights, base;
        uint32_t idx = sampleDiscrete(weights, &sumWeights, &base, m_numComponents, smp.uComponent);
        const BSDF* selectedBSDF = m_BSDFs[idx];
        if (sumWeights == 0.0f) {
            result->dirPDF = 0.0f;
            return SampledSpectrum::Zero;
        }
        
        SampledSpectrum value = selectedBSDF->sampleInternal(query, smp, result);
        result->dirPDF *= weights[idx];
        
        if (!result->dirType.isDelta()) {
            BSDFQuery mQuery = query;
            mQuery.flags &= sideTest(query.gNormal_sn, query.dir_sn, result->dir_sn);
            value = SampledSpectrum::Zero;
            for (int i = 0; i < m_numComponents; ++i) {
                const BSDF* BSDF = m_BSDFs[i];
                if (!BSDF->matches(mQuery.flags))
                    continue;
                value += BSDF->evaluateInternal(mQuery, result->dir_sn, nullptr);
                if (i != idx)
                    result->dirPDF += BSDF->evaluatePDFInternal(mQuery, result->dir_sn, nullptr) * weights[i];
            }
        }
        result->dirPDF /= sumWeights;
        
        return value;
    }
    
    SampledSpectrum MultiBSDF::sampleInternalWithRev(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult *result) const {
        float weights[maxNumElems];
        for (int i = 0; i < m_numComponents; ++i)
            weights[i] = m_BSDFs[i]->weight(query, smp);
        
        float sumWeights, base;
        uint32_t idx = sampleDiscrete(weights, &sumWeights, &base, m_numComponents, smp.uComponent);
        const BSDF* selectedBSDF = m_BSDFs[idx];
        if (sumWeights == 0.0f) {
            result->dirPDF = 0.0f;
            return SampledSpectrum::Zero;
        }
        
        SampledSpectrum value = selectedBSDF->sampleInternal(query, smp, result);
        
        BSDFQuery revQuery = query;// mQuery?
        Vector3D revDirIn = result->dir_sn;
        std::swap(revQuery.dir_sn, revDirIn);
        revQuery.adjoint ^= true;
        float revWeights[maxNumElems];
        FloatSum sumRevWeights = 0;
        for (int i = 0; i < m_numComponents; ++i) {
            revWeights[i] = m_BSDFs[i]->weight(revQuery, revDirIn);
            sumRevWeights += revWeights[i];
        }
        
        result->dirPDF *= weights[idx];
        result->reverse->dirPDF *= revWeights[idx];
        
        if (!result->dirType.isDelta()) {
            BSDFQuery mQuery = query;
            mQuery.flags &= sideTest(query.gNormal_sn, query.dir_sn, result->dir_sn);
            value = SampledSpectrum::Zero;

            for (int i = 0; i < m_numComponents; ++i) {
                const BSDF* BSDF = m_BSDFs[i];
                if (!BSDF->matches(mQuery.flags))
                    continue;
                SampledSpectrum rev_fs;
                float revPDF;
                value += BSDF->evaluateInternal(mQuery, result->dir_sn, &rev_fs);
                result->reverse->fs += rev_fs;
                if (i != idx) {
                    result->dirPDF += BSDF->evaluatePDFInternal(mQuery, result->dir_sn, &revPDF) * weights[i];
                    result->reverse->dirPDF += revPDF * revWeights[i];
                }
            }
        }
        result->dirPDF /= sumWeights;
        result->reverse->dirPDF /= sumRevWeights;
        
        return value;
    }
    
    SampledSpectrum MultiBSDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dirOut, SampledSpectrum* rev_fs) const {
        if (rev_fs) {
            SampledSpectrum retValue;
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
            SampledSpectrum retValue;
            for (int i = 0; i < m_numComponents; ++i) {
                if (!m_BSDFs[i]->matches(query.flags))
                    continue;
                retValue += m_BSDFs[i]->evaluateInternal(query, dirOut, nullptr);
            }
            return retValue;
        }
    }
    
    float MultiBSDF::evaluatePDFInternalNoRev(const BSDFQuery &query, const Vector3D &dirOut, float* revPDF) const {
        FloatSum sumWeights = 0.0f;
        float weights[maxNumElems];
        for (int i = 0; i < m_numComponents; ++i) {
            weights[i] = m_BSDFs[i]->weight(query, dirOut);
            sumWeights += weights[i];
        }
        if (sumWeights == 0.0f)
            return 0.0f;
        
        float retPDF = 0.0f;
        for (int i = 0; i < m_numComponents; ++i) {
            if (weights[i] > 0)
                retPDF += m_BSDFs[i]->evaluatePDFInternal(query, dirOut, nullptr) * weights[i];
        }
        retPDF /= sumWeights;
        
        return retPDF;
    }
    
    float MultiBSDF::evaluatePDFInternalWithRev(const BSDFQuery &query, const Vector3D &dirOut, float* revPDF) const {
        FloatSum sumWeights = 0.0f;
        FloatSum sumRevWeights = 0.0f;
        float weights[maxNumElems];
        float revWeights[maxNumElems];
        for (int i = 0; i < m_numComponents; ++i) {
            weights[i] = m_BSDFs[i]->weight(query, dirOut, &revWeights[i]);
            sumWeights += weights[i];
            sumRevWeights += revWeights[i];
        }
        if (sumWeights == 0.0f) {
            *revPDF = 0.0f;
            return 0.0f;
        }
        
        float retPDF = 0.0f;
        *revPDF = 0.0f;
        for (int i = 0; i < m_numComponents; ++i) {
            if (weights[i] > 0) {
                float eRevPDF;
                retPDF += m_BSDFs[i]->evaluatePDFInternal(query, dirOut, &eRevPDF) * weights[i];
                *revPDF += eRevPDF * revWeights[i];
            }
        }
        retPDF /= sumWeights;
        *revPDF /= sumRevWeights;
        
        return retPDF;
    }
    
    float MultiBSDF::weightInternal(const BSDFQuery &query, const BSDFSample &smp) const {
        FloatSum sumWeights = 0.0f;
        for (int i = 0; i < m_numComponents; ++i)
            sumWeights += m_BSDFs[i]->weight(query, smp);
        return sumWeights;
    }
    
    float MultiBSDF::weightInternal(const BSDFQuery &query, const Vector3D &dir, float* revWeight) const {
        if (revWeight) {
            FloatSum sumWeights = 0.0f;
            FloatSum sumRevWeights = 0.0f;
            for (int i = 0; i < m_numComponents; ++i) {
                float eRevWeights = 0.0f;
                sumWeights += m_BSDFs[i]->weight(query, dir, &eRevWeights);
                sumRevWeights += eRevWeights;
            }
            *revWeight = sumRevWeights;
            return sumWeights;
        }
        else {
            FloatSum sumWeights = 0.0f;
            for (int i = 0; i < m_numComponents; ++i)
                sumWeights += m_BSDFs[i]->weight(query, dir);
            return sumWeights;
        }
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

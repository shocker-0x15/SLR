//
//  MultiEDF.cpp
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "MultiEDF.h"

#include "../Core/distributions.h"

namespace SLR {
    void MultiEDF::add(EDF *edf) {
        SLRAssert(m_numComponents < maxNumEDFElems, "Number of MultiBSDF elements exceeds the limit: %u", maxNumEDFElems);
        m_EDFs[m_numComponents++] = edf;
    }
    
    SampledSpectrum MultiEDF::sample(const EDFQuery &query, const EDFSample &smp, EDFQueryResult *result) const {
        float weights[maxNumEDFElems];
        for (int i = 0; i < m_numComponents; ++i)
            weights[i] = m_EDFs[i]->weight(query);
        
        float tempProb, tempRemapped;
        float sumWeights;
        uint32_t idx = sampleDiscrete(weights, m_numComponents, smp.uComponent, &tempProb, &sumWeights, &tempRemapped);
        EDF* selectedEDF = m_EDFs[idx];
        
        // TODO: add the sampleInternal interface as in BSDF.
        SampledSpectrum value = selectedEDF->sample(query, smp, result);
        result->dirPDF *= weights[idx];
        
        if (!result->dirType.isDelta()) {
            EDFQuery mQuery = query;
            for (int i = 0; i < m_numComponents; ++i) {
                value += m_EDFs[i]->evaluate(mQuery, result->dir_sn);
                if (i != idx)
                    result->dirPDF += m_EDFs[i]->evaluatePDF(query, result->dir_sn) * weights[i];
            }
        }
        result->dirPDF /= sumWeights;
        
        return value;
    }
    
    SampledSpectrum MultiEDF::evaluate(const EDFQuery &query, const Vector3D &dir) const {
        SampledSpectrum retValue;
        for (int i = 0; i < m_numComponents; ++i)
            retValue += m_EDFs[i]->evaluate(query, dir);
        return retValue;
    }
    
    float MultiEDF::evaluatePDF(const EDFQuery &query, const Vector3D &dir) const {
        float sumWeights = 0.0f;
        float weights[maxNumEDFElems];
        for (int i = 0; i < m_numComponents; ++i) {
            weights[i] = m_EDFs[i]->weight(query);
            sumWeights += weights[i];
        }
        
        float retPDF = 0.0f;
        for (int i = 0; i < m_numComponents; ++i) {
            retPDF += m_EDFs[i]->evaluatePDF(query, dir) * weights[i];
        }
        retPDF /= sumWeights;
        
        return retPDF;
    }
    
    float MultiEDF::weight(const EDFQuery &query) const {
        float sumWeights = 0.0f;
        for (int i = 0; i < m_numComponents; ++i)
            sumWeights += m_EDFs[i]->weight(query);
        return sumWeights;
    }    
}

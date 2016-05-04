//
//  MicrofacetBSDF.cpp
//  SLR
//
//  Created by 渡部 心 on 2016/05/03.
//  Copyright © 2016年 渡部 心. All rights reserved.
//

#include "MicrofacetBSDF.h"

namespace SLR {
    SampledSpectrum MicrofacetBRDF::sampleInternal(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const {
        Normal3D m;
        float mPDF;
        float D = m_D->sample(smp.uDir[0], smp.uDir[1], &m, &mPDF);
        if (query.dir_sn.z < 0)
            m = -m;
        float dotHV = dot(query.dir_sn, m);
        result->dir_sn = 2 * dotHV * m - query.dir_sn;
        if (result->dir_sn.z * query.dir_sn.z <= 0) {
            result->dirPDF = 0.0f;
            return SampledSpectrum::Zero;
        }
        SLRAssert(dotHV >= 0, "dot between half and given vectors must be greater or equal to 0 at this point.");
        result->dirPDF = mPDF / (4 * dotHV);
        result->dirType = m_type;
        
        SampledSpectrum F = m_F->evaluate(dotHV);
        float G = m_D->evaluateSmithG1(query.dir_sn, m) * m_D->evaluateSmithG1(result->dir_sn, m);
        SampledSpectrum fs = F * D * G / (4 * query.dir_sn.z * result->dir_sn.z);
        if (result->reverse) {
            result->reverse->dirPDF = result->dirPDF;
            result->reverse->fs = fs;
        }
        
        return fs;
    }
    
    SampledSpectrum MicrofacetBRDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const {
        if (dir.z * query.dir_sn.z <= 0) {
            if (rev_fs)
                *rev_fs = SampledSpectrum::Zero;
            return SampledSpectrum::Zero;
        }
        
        Normal3D m = halfVector(query.dir_sn, dir);
        float dotHV = dot(query.dir_sn, m);
        if (query.dir_sn.z < 0)
            m = -m;
        float D = m_D->evaluate(m);
        
        SampledSpectrum F = m_F->evaluate(dotHV);
        float G = m_D->evaluateSmithG1(query.dir_sn, m) * m_D->evaluateSmithG1(dir, m);
        SampledSpectrum fs = F * D * G / (4 * query.dir_sn.z * dir.z);
        
        if (rev_fs)
            *rev_fs = fs;
        return fs;
    }
    
    float MicrofacetBRDF::evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, float* revPDF) const {
        if (dir.z * query.dir_sn.z <= 0) {
            if (revPDF)
                *revPDF = 0.0f;
            return 0.0f;
        }
        
        Normal3D m = halfVector(query.dir_sn, dir);
        float dotHV = dot(query.dir_sn, m);
        if (query.dir_sn.z < 0)
            m = -m;
        float mPDF = m_D->evaluatePDF(m);
        float ret = mPDF / (4 * dotHV);
        
        if (revPDF)
            *revPDF = ret;
        return ret;
    }
    
    float MicrofacetBRDF::weightInternal(const BSDFQuery &query) const {
        float F = m_F->evaluate(query.dir_sn.z).importance(query.wlHint);
        return F;
    }
    
    SampledSpectrum MicrofacetBRDF::getBaseColorInternal(DirectionType flags) const {
        return m_F->evaluate(1.0f);
    }
    
    
    SampledSpectrum MicrofacetBTDF::sampleInternal(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const {
        return SampledSpectrum::Zero;
    }
    
    SampledSpectrum MicrofacetBTDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const {
        return SampledSpectrum::Zero;
    }
    
    float MicrofacetBTDF::evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, float* revPDF) const {
        return 0.0f;
    }
    
    float MicrofacetBTDF::weightInternal(const BSDFQuery &query) const {
        float F = m_F.evaluate(query.dir_sn.z).importance(query.wlHint);
        return (1 - F);
    }
    
    SampledSpectrum MicrofacetBTDF::getBaseColorInternal(DirectionType flags) const {
        return SampledSpectrum::One - m_F.evaluate(1.0f);
    }
}

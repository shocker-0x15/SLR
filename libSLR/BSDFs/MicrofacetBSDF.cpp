//
//  MicrofacetBSDF.cpp
//  SLR
//
//  Created by 渡部 心 on 2016/05/03.
//  Copyright © 2016年 渡部 心. All rights reserved.
//

#include "MicrofacetBSDF.h"

namespace SLR {
    SampledSpectrum MicrofacetBRDF::sampleInternal(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult* result) const {
        bool entering = query.dir_sn.z >= 0.0f;
        int32_t sign = entering ? 1 : -1;
        
        Normal3D m;
        float mPDF;
        float D = m_D->sample(uDir[0], uDir[1], &m, &mPDF);
        float dotHV = dot(query.dir_sn, m);
        if (dotHV * sign <= 0) {
            result->dirPDF = 0.0f;
            return SampledSpectrum::Zero;
        }
        result->dir_sn = 2 * dotHV * m - query.dir_sn;
        result->dirPDF = mPDF / (4 * dotHV * sign);
        result->dirType = m_type;
        
        SampledSpectrum F = m_F.evaluate(dotHV);
        float G = m_D->evaluateSmithG1(query.dir_sn * sign, m) * m_D->evaluateSmithG1(result->dir_sn * sign, m);
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
        
        bool entering = query.dir_sn.z >= 0.0f;
        int32_t sign = entering ? 1 : -1;
        
        Normal3D m = sign * halfVector(query.dir_sn, dir);
        float dotHV = dot(query.dir_sn, m);
        float D = m_D->evaluate(m);
        
        SampledSpectrum F = m_F.evaluate(dotHV);
        float G = m_D->evaluateSmithG1(query.dir_sn * sign, m) * m_D->evaluateSmithG1(dir * sign, m);
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
        
        bool entering = query.dir_sn.z >= 0.0f;
        int32_t sign = entering ? 1 : -1;
        
        Normal3D m = sign * halfVector(query.dir_sn, dir);
        float dotHV = dot(query.dir_sn, m);
        if (dotHV * sign <= 0) {
            if (revPDF)
                *revPDF = 0.0f;
            return 0.0f;
        }
        float mPDF = m_D->evaluatePDF(m);
        float ret = mPDF / (4 * dotHV * sign);
        
        if (revPDF)
            *revPDF = ret;
        return ret;
    }
    
    float MicrofacetBRDF::weightInternal(const BSDFQuery &query) const {
        bool entering = query.dir_sn.z >= 0.0f;
        int32_t sign = entering ? 1 : -1;
        return m_D->evaluateSmithG1(query.dir_sn * sign, Normal3D(0, 0, 1));
    }
    
    SampledSpectrum MicrofacetBRDF::getBaseColorInternal(DirectionType flags) const {
        return m_F.evaluate(1.0f);
    }
    
    
    SampledSpectrum MicrofacetBSDF::sampleInternal(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult* result) const {
        bool entering = query.dir_sn.z >= 0.0f;
        int32_t sign = entering ? 1 : -1;
        const SampledSpectrum &eEnter = entering ? m_F.etaExt() : m_F.etaInt();
        const SampledSpectrum &eExit = entering ? m_F.etaInt() : m_F.etaExt();
        
        Normal3D m;
        float mPDF;
        float D = m_D->sample(uDir[0], uDir[1], &m, &mPDF);
        float dotHV = dot(query.dir_sn, m);
        if (dotHV * sign <= 0) {
            result->dirPDF = 0.0f;
            return SampledSpectrum::Zero;
        }
        
        SampledSpectrum F = m_F.evaluate(dotHV);
        float reflectProb = F.importance(query.wlHint);
        if (uComponent < reflectProb) {
            result->dir_sn = 2 * dotHV * m - query.dir_sn;
            if (result->dir_sn.z * query.dir_sn.z <= 0) {
                result->dirPDF = 0.0f;
                return SampledSpectrum::Zero;
            }
            result->dirPDF = reflectProb * mPDF / (4 * dotHV * sign);
            result->dirType = DirectionType::Reflection | DirectionType::HighFreq;
            
            SampledSpectrum F = m_F.evaluate(dotHV);
            float G = m_D->evaluateSmithG1(query.dir_sn * sign, m) * m_D->evaluateSmithG1(result->dir_sn * sign, m);
            SampledSpectrum fs = F * D * G / (4 * query.dir_sn.z * result->dir_sn.z);
            if (result->reverse) {
                result->reverse->dirPDF = result->dirPDF;
                result->reverse->fs = fs;
            }
            
            return fs;
        }
        else {
            float recRelIOR = eEnter[query.wlHint] / eExit[query.wlHint];
            float innerRoot = 1 + recRelIOR * recRelIOR * (dotHV * dotHV - 1);
            if (innerRoot < 0) {
                result->dirPDF = 0.0f;
                return SampledSpectrum::Zero;
            }
            result->dir_sn = (recRelIOR * dotHV - sign * std::sqrt(innerRoot)) * m - recRelIOR * query.dir_sn;
            if (result->dir_sn.z * query.dir_sn.z >= 0) {
                result->dirPDF = 0.0f;
                return SampledSpectrum::Zero;
            }
            float dotHL = dot(result->dir_sn, m);
            float commonPDFTerm = (1 - reflectProb) * mPDF / std::pow(eEnter[query.wlHint] * dotHV + eExit[query.wlHint] * dotHL, 2);
            result->dirPDF = commonPDFTerm * eExit[query.wlHint] * eExit[query.wlHint] * std::fabs(dotHL);
            result->dirType = DirectionType::Transmission | DirectionType::HighFreq;
            
            SampledSpectrum ret = SampledSpectrum::Zero;
            for (int wlIdx = 0; wlIdx < SampledSpectrum::NumComponents; ++wlIdx) {
                Normal3D m_wl = normalize(-(eEnter[wlIdx] * query.dir_sn + eExit[wlIdx] * result->dir_sn));
                float dotHV_wl = dot(query.dir_sn, m_wl);
                float dotHL_wl = dot(result->dir_sn, m_wl);
                float F_wl = m_F.evaluate(dotHV_wl, wlIdx);
                float G_wl = m_D->evaluateSmithG1(query.dir_sn * sign, m_wl) * m_D->evaluateSmithG1(result->dir_sn * -sign, m_wl);
                float D_wl = m_D->evaluate(m_wl);
                ret[wlIdx] = std::fabs(dotHV_wl * dotHL_wl) * (1 - F_wl) * G_wl * D_wl / std::pow(eEnter[wlIdx] * dotHV_wl + eExit[wlIdx] * dotHL_wl, 2);
            }
            ret /= std::fabs(query.dir_sn.z * result->dir_sn.z);
            
            if (result->reverse) {
                result->reverse->dirPDF = commonPDFTerm * eEnter[query.wlHint] * eEnter[query.wlHint] * std::fabs(dotHV);
                result->reverse->fs = eEnter * eEnter * ret;
            }
            return eExit * eExit * ret;
        }
    }
    
    SampledSpectrum MicrofacetBSDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const {
        bool entering = query.dir_sn.z >= 0.0f;
        int32_t sign = entering ? 1 : -1;
        
        if (dir.z * query.dir_sn.z > 0) {
            Normal3D m = sign * halfVector(query.dir_sn, dir);
            float dotHV = dot(query.dir_sn, m);
            float D = m_D->evaluate(m);
            
            SampledSpectrum F = m_F.evaluate(dotHV);
            float G = m_D->evaluateSmithG1(query.dir_sn * sign, m) * m_D->evaluateSmithG1(dir * sign, m);
            SampledSpectrum fs = F * D * G / (4 * query.dir_sn.z * dir.z);
            
            if (rev_fs)
                *rev_fs = fs;
            return fs;
        }
        else {
            const SampledSpectrum &eEnter = entering ? m_F.etaExt() : m_F.etaInt();
            const SampledSpectrum &eExit = entering ? m_F.etaInt() : m_F.etaExt();
            
            SampledSpectrum ret = SampledSpectrum::Zero;
            for (int wlIdx = 0; wlIdx < SampledSpectrum::NumComponents; ++wlIdx) {
                Normal3D m_wl = normalize(-(eEnter[wlIdx] * query.dir_sn + eExit[wlIdx] * dir));
                float dotHV_wl = dot(query.dir_sn, m_wl);
                float dotHL_wl = dot(dir, m_wl);
                float F_wl = m_F.evaluate(dotHV_wl, wlIdx);
                float G_wl = m_D->evaluateSmithG1(query.dir_sn * sign, m_wl) * m_D->evaluateSmithG1(dir * -sign, m_wl);
                float D_wl = m_D->evaluate(m_wl);
                ret[wlIdx] = std::fabs(dotHV_wl * dotHL_wl) * (1 - F_wl) * G_wl * D_wl / std::pow(eEnter[wlIdx] * dotHV_wl + eExit[wlIdx] * dotHL_wl, 2);
            }
            ret /= std::fabs(query.dir_sn.z * dir.z);
            
            if (rev_fs)
                *rev_fs = eEnter * eEnter * ret;
            return eExit * eExit * ret;
        }
    }
    
    float MicrofacetBSDF::evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, float* revPDF) const {
        bool entering = query.dir_sn.z >= 0.0f;
        int32_t sign = entering ? 1 : -1;
        const SampledSpectrum &eEnter = entering ? m_F.etaExt() : m_F.etaInt();
        const SampledSpectrum &eExit = entering ? m_F.etaInt() : m_F.etaExt();
        
        Normal3D m;
        if (dir.z * query.dir_sn.z > 0)
            m = sign * halfVector(query.dir_sn, dir);
        else
            m = normalize(-(eEnter[query.wlHint] * query.dir_sn + eExit[query.wlHint] * dir));
        float dotHV = dot(query.dir_sn, m);
        if (dotHV * sign <= 0) {
            if (revPDF)
                *revPDF = 0.0f;
            return 0.0f;
        }
        float mPDF = m_D->evaluatePDF(m);
        
        SampledSpectrum F = m_F.evaluate(dotHV);
        float reflectProb = F.importance(query.wlHint);
        if (dir.z * query.dir_sn.z > 0) {
            float ret = reflectProb * mPDF / (4 * dotHV * sign);
            
            if (revPDF)
                *revPDF = ret;
            return ret;
        }
        else {
            float dotHL = dot(dir, m);
            float commonPDFTerm = (1 - reflectProb) * mPDF / std::pow(eEnter[query.wlHint] * dotHV + eExit[query.wlHint] * dotHL, 2);
            float ret = commonPDFTerm * eExit[query.wlHint] * eExit[query.wlHint] * std::fabs(dotHL);
            
            if (revPDF)
                *revPDF = commonPDFTerm * eEnter[query.wlHint] * eEnter[query.wlHint] * std::fabs(dotHV);
            return ret;
        }
    }
    
    float MicrofacetBSDF::weightInternal(const BSDFQuery &query) const {
        bool entering = query.dir_sn.z >= 0.0f;
        int32_t sign = entering ? 1 : -1;
        return m_D->evaluateSmithG1(query.dir_sn * sign, Normal3D(0, 0, 1));
    }
    
    SampledSpectrum MicrofacetBSDF::getBaseColorInternal(DirectionType flags) const {
        return SampledSpectrum::One;
    }
}

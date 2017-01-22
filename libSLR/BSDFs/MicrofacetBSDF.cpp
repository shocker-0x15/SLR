//
//  MicrofacetBSDF.cpp
//
//  Created by 渡部 心 on 2016/05/03.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "MicrofacetBSDF.h"

namespace SLR {
    SampledSpectrum MicrofacetBRDF::sampleInternal(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult* result) const {
        bool entering = query.dirLocal.z >= 0.0f;
        int32_t sign = entering ? 1 : -1;
        
        Normal3D m;
        float mPDF;
        float D = m_D->sample(sign * query.dirLocal, uDir[0], uDir[1], &m, &mPDF);
        float dotHV = dot(query.dirLocal, m);
        if (dotHV * sign <= 0) {
            result->dirPDF = 0.0f;
            return SampledSpectrum::Zero;
        }
        result->dirLocal = 2 * dotHV * m - query.dirLocal;
        if (result->dirLocal.z * query.dirLocal.z <= 0) {
            result->dirPDF = 0.0f;
            return SampledSpectrum::Zero;
        }
        float commonPDFTerm = 1.0f / (4 * dotHV * sign);
        result->dirPDF = commonPDFTerm * mPDF;
        result->sampledType = m_type;
        
        SampledSpectrum F = m_F.evaluate(dotHV);
        float G = m_D->evaluateSmithG1(query.dirLocal, m) * m_D->evaluateSmithG1(result->dirLocal, m);
        SampledSpectrum fs = F * D * G / (4 * query.dirLocal.z * result->dirLocal.z);
        if (result->reverse) {
            result->reverse->dirPDF = commonPDFTerm * m_D->evaluatePDF(sign * result->dirLocal, m);
            result->reverse->fs = fs;
        }
        
        SLRAssert(fs.allFinite(), "fs: %s, F: %s, G, %g, D: %g, wlIdx: %u, qDir: %s, rDir: %s",
                  fs.toString().c_str(), F.toString().c_str(), G, D, query.wlHint, query.dirLocal.toString().c_str(), result->dirLocal.toString().c_str());
        
        return fs;
    }
    
    SampledSpectrum MicrofacetBRDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const {
        if (dir.z * query.dirLocal.z <= 0) {
            if (rev_fs)
                *rev_fs = SampledSpectrum::Zero;
            return SampledSpectrum::Zero;
        }
        
        bool entering = query.dirLocal.z >= 0.0f;
        int32_t sign = entering ? 1 : -1;
        
        Normal3D m = sign * halfVector(query.dirLocal, dir);
        float dotHV = dot(query.dirLocal, m);
        float D = m_D->evaluate(m);
        
        SampledSpectrum F = m_F.evaluate(dotHV);
        float G = m_D->evaluateSmithG1(query.dirLocal, m) * m_D->evaluateSmithG1(dir, m);
        SampledSpectrum fs = F * D * G / (4 * query.dirLocal.z * dir.z);
        
        SLRAssert(fs.allFinite(), "fs: %s, F: %s, G, %g, D: %g, wlIdx: %u, qDir: %s, dir: %s",
                  fs.toString().c_str(), F.toString().c_str(), G, D, query.wlHint, query.dirLocal.toString().c_str(), dir.toString().c_str());
        
        if (rev_fs)
            *rev_fs = fs;
        return fs;
    }
    
    float MicrofacetBRDF::evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, float* revPDF) const {
        if (dir.z * query.dirLocal.z <= 0) {
            if (revPDF)
                *revPDF = 0.0f;
            return 0.0f;
        }
        
        bool entering = query.dirLocal.z >= 0.0f;
        int32_t sign = entering ? 1 : -1;
        
        Normal3D m = sign * halfVector(query.dirLocal, dir);
        float dotHV = dot(query.dirLocal, m);
        if (dotHV * sign <= 0) {
            if (revPDF)
                *revPDF = 0.0f;
            return 0.0f;
        }
        float mPDF = m_D->evaluatePDF(sign * query.dirLocal, m);
        float commonPDFTerm = 1.0f / (4 * dotHV * sign);
        float ret = commonPDFTerm * mPDF;
        
        SLRAssert(std::isfinite(commonPDFTerm) && std::isfinite(mPDF),
                  "commonPDFTerm: %g, mPDF: %g, wlIdx: %u, qDir: %s, dir: %s",
                  commonPDFTerm, mPDF, query.wlHint, query.dirLocal.toString().c_str(), dir.toString().c_str());
        
        if (revPDF)
            *revPDF = commonPDFTerm * m_D->evaluatePDF(sign * dir, m);
        return ret;
    }
    
    float MicrofacetBRDF::weightInternal(const BSDFQuery &query) const {
        bool entering = query.dirLocal.z >= 0.0f;
        int32_t sign = entering ? 1 : -1;
        return m_D->evaluateSmithG1(query.dirLocal * sign, Normal3D(0, 0, 1));
    }
    
    SampledSpectrum MicrofacetBRDF::getBaseColorInternal(DirectionType flags) const {
        return m_F.evaluate(1.0f);
    }
    
    
    SampledSpectrum MicrofacetBSDF::sampleInternal(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult* result) const {
        bool entering = query.dirLocal.z >= 0.0f;
        int32_t sign = entering ? 1 : -1;
        const SampledSpectrum &eEnter = entering ? m_F.etaExt() : m_F.etaInt();
        const SampledSpectrum &eExit = entering ? m_F.etaInt() : m_F.etaExt();
        
        Normal3D m;
        float mPDF;
        float D = m_D->sample(sign * query.dirLocal, uDir[0], uDir[1], &m, &mPDF);
        float dotHV = dot(query.dirLocal, m);
        if (dotHV * sign <= 0 || std::isnan(D)) {
            result->dirPDF = 0.0f;
            return SampledSpectrum::Zero;
        }
        
        SampledSpectrum F = m_F.evaluate(dotHV);
        float reflectProb = F.importance(query.wlHint);
        if (query.dirTypeFilter.isReflection())
            reflectProb = 1.0f;
        if (query.dirTypeFilter.isTransmission())
            reflectProb = 0.0f;
        if (uComponent < reflectProb) {
            result->dirLocal = 2 * dotHV * m - query.dirLocal;
            if (result->dirLocal.z * query.dirLocal.z <= 0) {
                result->dirPDF = 0.0f;
                return SampledSpectrum::Zero;
            }
            float commonPDFTerm = reflectProb / (4 * dotHV * sign);
            result->dirPDF = commonPDFTerm * mPDF;
            result->sampledType = DirectionType::Reflection | DirectionType::HighFreq;
            
            float G = m_D->evaluateSmithG1(query.dirLocal, m) * m_D->evaluateSmithG1(result->dirLocal, m);
            SampledSpectrum fs = F * D * G / (4 * query.dirLocal.z * result->dirLocal.z);
            
            SLRAssert(fs.allFinite(), "fs: %s, F: %s, G, %g, D: %g, wlIdx: %u, qDir: %s, rDir: %s",
                      fs.toString().c_str(), F.toString().c_str(), G, D, query.wlHint, query.dirLocal.toString().c_str(), result->dirLocal.toString().c_str());
            
            if (result->reverse) {
                result->reverse->dirPDF = commonPDFTerm * m_D->evaluatePDF(sign * result->dirLocal, m);
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
            result->dirLocal = (recRelIOR * dotHV - sign * std::sqrt(innerRoot)) * m - recRelIOR * query.dirLocal;
            if (result->dirLocal.z * query.dirLocal.z >= 0) {
                result->dirPDF = 0.0f;
                return SampledSpectrum::Zero;
            }
            float dotHL = dot(result->dirLocal, m);
            float commonPDFTerm = (1 - reflectProb) / std::pow(eEnter[query.wlHint] * dotHV + eExit[query.wlHint] * dotHL, 2);
            result->dirPDF = commonPDFTerm * mPDF * eExit[query.wlHint] * eExit[query.wlHint] * std::fabs(dotHL);
            result->sampledType = DirectionType::Transmission | DirectionType::HighFreq;
            
            SampledSpectrum ret = SampledSpectrum::Zero;
            for (int wlIdx = 0; wlIdx < SampledSpectrum::NumComponents; ++wlIdx) {
                Normal3D m_wl = normalize(-(eEnter[wlIdx] * query.dirLocal + eExit[wlIdx] * result->dirLocal));
                float dotHV_wl = dot(query.dirLocal, m_wl);
                float dotHL_wl = dot(result->dirLocal, m_wl);
                float F_wl = m_F.evaluate(dotHV_wl, wlIdx);
                float G_wl = m_D->evaluateSmithG1(query.dirLocal, m_wl) * m_D->evaluateSmithG1(result->dirLocal, m_wl);
                float D_wl = m_D->evaluate(m_wl);
                ret[wlIdx] = std::fabs(dotHV_wl * dotHL_wl) * (1 - F_wl) * G_wl * D_wl / std::pow(eEnter[wlIdx] * dotHV_wl + eExit[wlIdx] * dotHL_wl, 2);
                
                SLRAssert(std::isfinite(ret[wlIdx]), "fs: %g, F: %g, G, %g, D: %g, wlIdx: %u, qDir: %s",
                          ret[wlIdx], F_wl, G_wl, D_wl, query.wlHint, query.dirLocal.toString().c_str());
            }
            ret /= std::fabs(query.dirLocal.z * result->dirLocal.z);
            ret *= query.adjoint ? (eExit * eExit) : (eEnter * eEnter);// adjoint: need to cancel eEnter^2 / eExit^2 => eEnter^2 * (eExit^2 / eEnter^2)
            
            SLRAssert(ret.allFinite(), "fs: %s, wlIdx: %u, qDir: %s, rDir: %s",
                      ret.toString().c_str(), query.wlHint, query.dirLocal.toString().c_str(), result->dirLocal.toString().c_str());
            
            if (result->reverse) {
                result->reverse->dirPDF = commonPDFTerm * m_D->evaluatePDF(-sign * result->dirLocal, m) * eEnter[query.wlHint] * eEnter[query.wlHint] * std::fabs(dotHV);
                result->reverse->fs = ret;
            }
            return ret;
        }
    }
    
    SampledSpectrum MicrofacetBSDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const {
        bool entering = query.dirLocal.z >= 0.0f;
        int32_t sign = entering ? 1 : -1;
        float dotNVdotNL = dir.z * query.dirLocal.z;
        
        if (dotNVdotNL > 0 && query.dirTypeFilter.matches(DirectionType::Reflection | DirectionType::AllFreq)) {
            Normal3D m = sign * halfVector(query.dirLocal, dir);
            float dotHV = dot(query.dirLocal, m);
            float D = m_D->evaluate(m);
            
            SampledSpectrum F = m_F.evaluate(dotHV);
            float G = m_D->evaluateSmithG1(query.dirLocal, m) * m_D->evaluateSmithG1(dir, m);
            SampledSpectrum fs = F * D * G / (4 * dotNVdotNL);
            
            SLRAssert(fs.allFinite(), "fs: %s, F: %s, G, %g, D: %g, wlIdx: %u, qDir: %s, dir: %s",
                      fs.toString().c_str(), F.toString().c_str(), G, D, query.wlHint, query.dirLocal.toString().c_str(), dir.toString().c_str());
            
            if (rev_fs)
                *rev_fs = fs;
            return fs;
        }
        else if (dotNVdotNL < 0 && query.dirTypeFilter.matches(DirectionType::Transmission | DirectionType::AllFreq)) {
            const SampledSpectrum &eEnter = entering ? m_F.etaExt() : m_F.etaInt();
            const SampledSpectrum &eExit = entering ? m_F.etaInt() : m_F.etaExt();
            
            SampledSpectrum ret = SampledSpectrum::Zero;
            for (int wlIdx = 0; wlIdx < SampledSpectrum::NumComponents; ++wlIdx) {
                Normal3D m_wl = normalize(-(eEnter[wlIdx] * query.dirLocal + eExit[wlIdx] * dir));
                float dotHV_wl = dot(query.dirLocal, m_wl);
                float dotHL_wl = dot(dir, m_wl);
                float F_wl = m_F.evaluate(dotHV_wl, wlIdx);
                float G_wl = m_D->evaluateSmithG1(query.dirLocal, m_wl) * m_D->evaluateSmithG1(dir, m_wl);
                float D_wl = m_D->evaluate(m_wl);
                ret[wlIdx] = std::fabs(dotHV_wl * dotHL_wl) * (1 - F_wl) * G_wl * D_wl / std::pow(eEnter[wlIdx] * dotHV_wl + eExit[wlIdx] * dotHL_wl, 2);
                
                SLRAssert(std::isfinite(ret[wlIdx]), "fs: %g, F: %g, G, %g, D: %g, wlIdx: %u, qDir: %s, dir: %s",
                          ret[wlIdx], F_wl, G_wl, D_wl, query.wlHint, query.dirLocal.toString().c_str(), dir.toString().c_str());
            }
            ret /= std::fabs(dotNVdotNL);
            ret *= query.adjoint ? (eExit * eExit) : (eEnter * eEnter);// !adjoint: eExit^2 * (eEnter / eExit)^2
            
            SLRAssert(ret.allFinite(), "fs: %s, wlIdx: %u, qDir: %s, dir: %s",
                      ret.toString().c_str(), query.wlHint, query.dirLocal.toString().c_str(), dir.toString().c_str());
            
            if (rev_fs)
                *rev_fs = ret;
            return ret;
        }
        
        if (rev_fs)
            *rev_fs = SampledSpectrum::Zero;
        return SampledSpectrum::Zero;
    }
    
    float MicrofacetBSDF::evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, float* revPDF) const {
        bool entering = query.dirLocal.z >= 0.0f;
        int32_t sign = entering ? 1 : -1;
        float dotNVdotNL = dir.z * query.dirLocal.z;
        if (dotNVdotNL == 0)
            return 0.0f;
        const SampledSpectrum &eEnter = entering ? m_F.etaExt() : m_F.etaInt();
        const SampledSpectrum &eExit = entering ? m_F.etaInt() : m_F.etaExt();
        
        Normal3D m;
        if (dotNVdotNL > 0)
            m = sign * halfVector(query.dirLocal, dir);
        else
            m = normalize(-(eEnter[query.wlHint] * query.dirLocal + eExit[query.wlHint] * dir));
        float dotHV = dot(query.dirLocal, m);
        if (dotHV * sign <= 0) {
            if (revPDF)
                *revPDF = 0.0f;
            return 0.0f;
        }
        float mPDF = m_D->evaluatePDF(sign * query.dirLocal, m);
        
        SampledSpectrum F = m_F.evaluate(dotHV);
        float reflectProb = F.importance(query.wlHint);
        if (query.dirTypeFilter.isReflection())
            reflectProb = 1.0f;
        if (query.dirTypeFilter.isTransmission())
            reflectProb = 0.0f;
        if (dotNVdotNL > 0) {
            float commonPDFTerm = reflectProb / (4 * dotHV * sign);
            
            SLRAssert(std::isfinite(commonPDFTerm) && std::isfinite(mPDF),
                      "commonPDFTerm: %g, mPDF: %g, F: %s, wlIdx: %u, qDir: %s, dir: %s",
                      commonPDFTerm, mPDF, F.toString().c_str(), query.wlHint, query.dirLocal.toString().c_str(), dir.toString().c_str());
            
            if (revPDF)
                *revPDF = commonPDFTerm * m_D->evaluatePDF(sign * dir, m);
            return commonPDFTerm * mPDF;
        }
        else {
            float dotHL = dot(dir, m);
            float commonPDFTerm = (1 - reflectProb) / std::pow(eEnter[query.wlHint] * dotHV + eExit[query.wlHint] * dotHL, 2);
            
            SLRAssert(std::isfinite(commonPDFTerm) && std::isfinite(mPDF),
                      "commonPDFTerm: %g, mPDF: %g, F: %s, wlIdx: %u, qDir: %s, dir: %s",
                      commonPDFTerm, mPDF, F.toString().c_str(), query.wlHint, query.dirLocal.toString().c_str(), dir.toString().c_str());
            
            if (revPDF)
                *revPDF = commonPDFTerm * m_D->evaluatePDF(-sign * dir, m) * eEnter[query.wlHint] * eEnter[query.wlHint] * std::fabs(dotHV);
            return commonPDFTerm * mPDF * eExit[query.wlHint] * eExit[query.wlHint] * std::fabs(dotHL);
        }
    }
    
    float MicrofacetBSDF::weightInternal(const BSDFQuery &query) const {
        bool entering = query.dirLocal.z >= 0.0f;
        int32_t sign = entering ? 1 : -1;
        return m_D->evaluateSmithG1(query.dirLocal * sign, Normal3D(0, 0, 1));
    }
    
    SampledSpectrum MicrofacetBSDF::getBaseColorInternal(DirectionType flags) const {
        return SampledSpectrum::One;
    }
}

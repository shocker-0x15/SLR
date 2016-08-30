//
//  MicrofacetBSDF.cpp
//
//  Created by 渡部 心 on 2016/05/03.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "MicrofacetBSDF.h"

namespace SLR {
    SampledSpectrum MicrofacetBRDF::sampleInternal(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult* result) const {
        bool entering = query.dir_sn.z >= 0.0f;
        int32_t sign = entering ? 1 : -1;
        
        Normal3D m;
        float mPDF;
        float D = m_D->sample(sign * query.dir_sn, uDir[0], uDir[1], &m, &mPDF);
        float dotHV = dot(query.dir_sn, m);
        if (dotHV * sign <= 0) {
            result->dirPDF = SampledSpectrum::Zero;
            return SampledSpectrum::Zero;
        }
        result->dir_sn = 2 * dotHV * m - query.dir_sn;
        if (result->dir_sn.z * query.dir_sn.z <= 0) {
            result->dirPDF = SampledSpectrum::Zero;
            return SampledSpectrum::Zero;
        }
        float commonPDFTerm = 1.0f / (4 * dotHV * sign);
        result->dirPDF = commonPDFTerm * mPDF;
        result->dirType = m_type;
        
        SampledSpectrum F = m_F.evaluate(dotHV);
        float G = m_D->evaluateSmithG1(query.dir_sn, m) * m_D->evaluateSmithG1(result->dir_sn, m);
        SampledSpectrum fs = F * D * G / (4 * query.dir_sn.z * result->dir_sn.z);
        if (result->reverse) {
            result->reverse->dirPDF = commonPDFTerm * m_D->evaluatePDF(sign * result->dir_sn, m);
            result->reverse->fs = fs;
        }
        
        SLRAssert(!fs.hasInf() && !fs.hasNaN(), "fs: %s, F: %s, G, %g, D: %g, wlIdx: %u, qDir: %s, rDir: %s",
                  fs.toString().c_str(), F.toString().c_str(), G, D, query.wlHint, query.dir_sn.toString().c_str(), result->dir_sn.toString().c_str());
        
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
        float G = m_D->evaluateSmithG1(query.dir_sn, m) * m_D->evaluateSmithG1(dir, m);
        SampledSpectrum fs = F * D * G / (4 * query.dir_sn.z * dir.z);
        
        SLRAssert(!fs.hasInf() && !fs.hasNaN(), "fs: %s, F: %s, G, %g, D: %g, wlIdx: %u, qDir: %s, dir: %s",
                  fs.toString().c_str(), F.toString().c_str(), G, D, query.wlHint, query.dir_sn.toString().c_str(), dir.toString().c_str());
        
        if (rev_fs)
            *rev_fs = fs;
        return fs;
    }
    
    SampledSpectrum MicrofacetBRDF::evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* revPDF) const {
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
                *revPDF = SampledSpectrum::Zero;
            return SampledSpectrum::Zero;
        }
        float mPDF = m_D->evaluatePDF(sign * query.dir_sn, m);
        float commonPDFTerm = 1.0f / (4 * dotHV * sign);
        float ret = commonPDFTerm * mPDF;
        
        SLRAssert(!std::isnan(commonPDFTerm) && !std::isinf(commonPDFTerm) && !std::isnan(mPDF) && !std::isinf(mPDF),
                  "commonPDFTerm: %g, mPDF: %g, wlIdx: %u, qDir: %s, dir: %s",
                  commonPDFTerm, mPDF, query.wlHint, query.dir_sn.toString().c_str(), dir.toString().c_str());
        
        if (revPDF)
            *revPDF = commonPDFTerm * m_D->evaluatePDF(sign * dir, m);
        return ret;
    }
    
    SampledSpectrum MicrofacetBRDF::weightInternal(const BSDFQuery &query) const {
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
        float D = m_D->sample(sign * query.dir_sn, uDir[0], uDir[1], &m, &mPDF);
        float dotHV = dot(query.dir_sn, m);
        if (dotHV * sign <= 0 || std::isnan(D)) {
            result->dirPDF = 0.0f;
            return SampledSpectrum::Zero;
        }
        
        SampledSpectrum F = m_F.evaluate(dotHV);
        float reflectProb = F[query.heroIndex];
        if (query.flags.isReflection())
            reflectProb = 1.0f;
        if (query.flags.isTransmission())
            reflectProb = 0.0f;
        if (uComponent < reflectProb) {
            result->dir_sn = 2 * dotHV * m - query.dir_sn;
            if (result->dir_sn.z * query.dir_sn.z <= 0) {
                result->dirPDF = SampledSpectrum::Zero;
                return SampledSpectrum::Zero;
            }
            SampledSpectrum commonPDFTerm = F / (4 * dotHV * sign);
            result->dirPDF = commonPDFTerm * mPDF;
            result->dirType = DirectionType::Reflection | DirectionType::HighFreq;
            
            float G = m_D->evaluateSmithG1(query.dir_sn, m) * m_D->evaluateSmithG1(result->dir_sn, m);
            SampledSpectrum fs = F * D * G / (4 * query.dir_sn.z * result->dir_sn.z);
            
            SLRAssert(!fs.hasInf() && !fs.hasNaN(), "fs: %s, F: %s, G, %g, D: %g, wlIdx: %u, qDir: %s, rDir: %s",
                      fs.toString().c_str(), F.toString().c_str(), G, D, query.wlHint, query.dir_sn.toString().c_str(), result->dir_sn.toString().c_str());
            
            if (result->reverse) {
                result->reverse->dirPDF = commonPDFTerm * m_D->evaluatePDF(sign * result->dir_sn, m);
                result->reverse->fs = fs;
            }
            return fs;
        }
        else {
            float recRelIOR = eEnter[query.heroIndex] / eExit[query.heroIndex];
            float innerRoot = 1 + recRelIOR * recRelIOR * (dotHV * dotHV - 1);
            if (innerRoot < 0) {
                result->dirPDF = SampledSpectrum::Zero;
                return SampledSpectrum::Zero;
            }
            result->dir_sn = (recRelIOR * dotHV - sign * std::sqrt(innerRoot)) * m - recRelIOR * query.dir_sn;
            if (result->dir_sn.z * query.dir_sn.z >= 0) {
                result->dirPDF = SampledSpectrum::Zero;
                return SampledSpectrum::Zero;
            }
            float dotHL = dot(result->dir_sn, m);
            SampledSpectrum commonPDFTerm = SampledSpectrum::One / pow(eEnter * dotHV + eExit * dotHL, 2);
            result->dirPDF = (SampledSpectrum::One - F) * (commonPDFTerm * mPDF * eExit * eExit * std::fabs(dotHL));
            result->dirType = DirectionType::Transmission | DirectionType::HighFreq;
            
            SampledSpectrum ret = SampledSpectrum::Zero;
            for (int wlIdx = 0; wlIdx < SampledSpectrum::NumComponents; ++wlIdx) {
                Normal3D m_wl = normalize(-(eEnter[wlIdx] * query.dir_sn + eExit[wlIdx] * result->dir_sn));
                float dotHV_wl = dot(query.dir_sn, m_wl);
                float dotHL_wl = dot(result->dir_sn, m_wl);
                float F_wl = m_F.evaluate(dotHV_wl, wlIdx);
                float G_wl = m_D->evaluateSmithG1(query.dir_sn, m_wl) * m_D->evaluateSmithG1(result->dir_sn, m_wl);
                float D_wl = m_D->evaluate(m_wl);
                ret[wlIdx] = std::fabs(dotHV_wl * dotHL_wl) * (1 - F_wl) * G_wl * D_wl / std::pow(eEnter[wlIdx] * dotHV_wl + eExit[wlIdx] * dotHL_wl, 2);
                
                SLRAssert(!std::isnan(ret[wlIdx]) && !std::isinf(ret[wlIdx]), "fs: %g, F: %g, G, %g, D: %g, wlIdx: %u, qDir: %s",
                          ret[wlIdx], F_wl, G_wl, D_wl, query.wlHint, query.dir_sn.toString().c_str());
            }
            ret /= std::fabs(query.dir_sn.z * result->dir_sn.z);
            ret *= query.adjoint ? (eExit * eExit) : (eEnter * eEnter);// adjoint: need to cancel eEnter^2 / eExit^2 => eEnter^2 * (eExit^2 / eEnter^2)
            
            SLRAssert(!ret.hasInf() && !ret.hasNaN(), "fs: %s, wlIdx: %u, qDir: %s, rDir: %s",
                      ret.toString().c_str(), query.wlHint, query.dir_sn.toString().c_str(), result->dir_sn.toString().c_str());
            
            if (result->reverse) {
                SampledSpectrum revF = m_F.evaluate(dotHL);
                result->reverse->dirPDF = (SampledSpectrum::One - revF) * (commonPDFTerm * m_D->evaluatePDF(-sign * result->dir_sn, m) * eEnter * eEnter * std::fabs(dotHV));
                result->reverse->fs = ret;
            }
            return ret;
        }
    }
    
    SampledSpectrum MicrofacetBSDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const {
        bool entering = query.dir_sn.z >= 0.0f;
        int32_t sign = entering ? 1 : -1;
        float dotNVdotNL = dir.z * query.dir_sn.z;
        
        if (dotNVdotNL > 0 && query.flags.matches(DirectionType::Reflection | DirectionType::AllFreq)) {
            Normal3D m = sign * halfVector(query.dir_sn, dir);
            float dotHV = dot(query.dir_sn, m);
            float D = m_D->evaluate(m);
            
            SampledSpectrum F = m_F.evaluate(dotHV);
            float G = m_D->evaluateSmithG1(query.dir_sn, m) * m_D->evaluateSmithG1(dir, m);
            SampledSpectrum fs = F * D * G / (4 * dotNVdotNL);
            
            SLRAssert(!fs.hasInf() && !fs.hasNaN(), "fs: %s, F: %s, G, %g, D: %g, wlIdx: %u, qDir: %s, dir: %s",
                      fs.toString().c_str(), F.toString().c_str(), G, D, query.wlHint, query.dir_sn.toString().c_str(), dir.toString().c_str());
            
            if (rev_fs)
                *rev_fs = fs;
            return fs;
        }
        else if (dotNVdotNL < 0 && query.flags.matches(DirectionType::Transmission | DirectionType::AllFreq)) {
            const SampledSpectrum &eEnter = entering ? m_F.etaExt() : m_F.etaInt();
            const SampledSpectrum &eExit = entering ? m_F.etaInt() : m_F.etaExt();
            
            SampledSpectrum ret = SampledSpectrum::Zero;
            for (int wlIdx = 0; wlIdx < SampledSpectrum::NumComponents; ++wlIdx) {
                Normal3D m_wl = normalize(-(eEnter[wlIdx] * query.dir_sn + eExit[wlIdx] * dir));
                float dotHV_wl = dot(query.dir_sn, m_wl);
                float dotHL_wl = dot(dir, m_wl);
                float F_wl = m_F.evaluate(dotHV_wl, wlIdx);
                float G_wl = m_D->evaluateSmithG1(query.dir_sn, m_wl) * m_D->evaluateSmithG1(dir, m_wl);
                float D_wl = m_D->evaluate(m_wl);
                ret[wlIdx] = std::fabs(dotHV_wl * dotHL_wl) * (1 - F_wl) * G_wl * D_wl / std::pow(eEnter[wlIdx] * dotHV_wl + eExit[wlIdx] * dotHL_wl, 2);
                
                SLRAssert(!std::isnan(ret[wlIdx]) && !std::isinf(ret[wlIdx]), "fs: %g, F: %g, G, %g, D: %g, wlIdx: %u, qDir: %s, dir: %s",
                          ret[wlIdx], F_wl, G_wl, D_wl, query.wlHint, query.dir_sn.toString().c_str(), dir.toString().c_str());
            }
            ret /= std::fabs(dotNVdotNL);
            ret *= query.adjoint ? (eExit * eExit) : (eEnter * eEnter);// !adjoint: eExit^2 * (eEnter / eExit)^2
            
            SLRAssert(!ret.hasInf() && !ret.hasNaN(), "fs: %s, wlIdx: %u, qDir: %s, dir: %s",
                      ret.toString().c_str(), query.wlHint, query.dir_sn.toString().c_str(), dir.toString().c_str());
            
            if (rev_fs)
                *rev_fs = ret;
            return ret;
        }
        
        if (rev_fs)
            *rev_fs = SampledSpectrum::Zero;
        return SampledSpectrum::Zero;
    }
    
    SampledSpectrum MicrofacetBSDF::evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* revPDF) const {
        bool entering = query.dir_sn.z >= 0.0f;
        int32_t sign = entering ? 1 : -1;
        float dotNVdotNL = dir.z * query.dir_sn.z;
        if (dotNVdotNL == 0)
            return SampledSpectrum::Zero;
        const SampledSpectrum &eEnter = entering ? m_F.etaExt() : m_F.etaInt();
        const SampledSpectrum &eExit = entering ? m_F.etaInt() : m_F.etaExt();
        
        Normal3D m;
        if (dotNVdotNL > 0)
            m = sign * halfVector(query.dir_sn, dir);
        else
            m = normalize(-(eEnter[query.heroIndex] * query.dir_sn + eExit[query.heroIndex] * dir));
        float dotHV = dot(query.dir_sn, m);
        if (dotHV * sign <= 0) {
            if (revPDF)
                *revPDF = SampledSpectrum::Zero;
            return SampledSpectrum::Zero;
        }
        float mPDF = m_D->evaluatePDF(sign * query.dir_sn, m);
        
        SampledSpectrum F = m_F.evaluate(dotHV);
        float reflectProb = F[query.heroIndex];
        if (query.flags.isReflection())
            reflectProb = 1.0f;
        if (query.flags.isTransmission())
            reflectProb = 0.0f;
        if (dotNVdotNL > 0) {
            SampledSpectrum commonPDFTerm = F / (4 * dotHV * sign);
            
            SLRAssert(!std::isnan(commonPDFTerm) && !std::isinf(commonPDFTerm) && !std::isnan(mPDF) && !std::isinf(mPDF),
                      "commonPDFTerm: %g, mPDF: %g, F: %s, wlIdx: %u, qDir: %s, dir: %s",
                      commonPDFTerm, mPDF, F.toString().c_str(), query.wlHint, query.dir_sn.toString().c_str(), dir.toString().c_str());
            
            if (revPDF)
                *revPDF = commonPDFTerm * m_D->evaluatePDF(sign * dir, m);
            return commonPDFTerm * mPDF;
        }
        else {
            float dotHL = dot(dir, m);
            SampledSpectrum commonPDFTerm = SampledSpectrum::One / pow(eEnter * dotHV + eExit * dotHL, 2);
            
            SLRAssert(!std::isnan(commonPDFTerm) && !std::isinf(commonPDFTerm) && !std::isnan(mPDF) && !std::isinf(mPDF),
                      "commonPDFTerm: %g, mPDF: %g, F: %s, wlIdx: %u, qDir: %s, dir: %s",
                      commonPDFTerm, mPDF, F.toString().c_str(), query.wlHint, query.dir_sn.toString().c_str(), dir.toString().c_str());
            
            if (revPDF) {
                SampledSpectrum revF = m_F.evaluate(dotHL);
                *revPDF = (SampledSpectrum::One - revF) * (commonPDFTerm * m_D->evaluatePDF(-sign * dir, m) * eEnter * eEnter * std::fabs(dotHV));
            }
            return (SampledSpectrum::One - F) * (commonPDFTerm * mPDF * eExit * eExit * std::fabs(dotHL));
        }
    }
    
    SampledSpectrum MicrofacetBSDF::weightInternal(const BSDFQuery &query) const {
        bool entering = query.dir_sn.z >= 0.0f;
        int32_t sign = entering ? 1 : -1;
        return m_D->evaluateSmithG1(query.dir_sn * sign, Normal3D(0, 0, 1));
    }
    
    SampledSpectrum MicrofacetBSDF::getBaseColorInternal(DirectionType flags) const {
        return SampledSpectrum::One;
    }
}

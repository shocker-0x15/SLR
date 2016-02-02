//
//  AshikhminBRDF.cpp
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "AshikhminBRDF.h"
#include "../Core/distributions.h"

namespace SLR {
    SampledSpectrum AshikhminSpecularBRDF::sampleInternal(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult *result) const {
        float quad = 2 * M_PI * smp.uDir[1];
        float phi_h = std::atan2(std::sqrt(m_nu + 1) * std::sin(quad), std::sqrt(m_nv + 1) * std::cos(quad));
        float cosphi = std::cos(phi_h);
        float sinphi = std::sin(phi_h);
        float theta_h = std::acos(std::pow(1 - smp.uDir[0], 1.0f / (m_nu * cosphi * cosphi + m_nv * sinphi * sinphi + 1)));
        if (query.dir_sn.z < 0)
            theta_h = M_PI - theta_h;
        Vector3D halfv = Vector3D(std::sin(theta_h) * std::cos(phi_h), std::sin(theta_h) * std::sin(phi_h), std::cos(theta_h));
        result->dir_sn = 2 * dot(query.dir_sn, halfv) * halfv - query.dir_sn;
        result->dirType = m_type;
        
        float dotHV = dot(halfv, query.dir_sn);
        float exp = (m_nu * halfv.x * halfv.x + m_nv * halfv.y * halfv.y) / (1 - halfv.z * halfv.z);
        SampledSpectrum F = m_Rs + (SampledSpectrum::One - m_Rs) * std::pow(1.0f - dotHV, 5);
        float commonTerm = std::sqrt((m_nu + 1) * (m_nv + 1)) / (8 * M_PI * dotHV) * std::pow(std::fabs(halfv.z), exp);
        SampledSpectrum fs = commonTerm / std::fmax(std::fabs(query.dir_sn.z), std::fabs(result->dir_sn.z)) * F;
        result->dirPDF = commonTerm;
        if (result->reverse) {
            result->reverse->fs = fs;
            result->reverse->dirPDF = commonTerm;
        }
        return fs;
    }
    
    SampledSpectrum AshikhminSpecularBRDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const {
        Vector3D halfv = halfVector(query.dir_sn, dir);
        float dotHV = dot(halfv, query.dir_sn);
        float exp = (m_nu * halfv.x * halfv.x + m_nv * halfv.y * halfv.y) / (1 - halfv.z * halfv.z);
        SampledSpectrum F = m_Rs + (SampledSpectrum::One - m_Rs) * std::pow(1.0f - dotHV, 5);
        float commonTerm = std::sqrt((m_nu + 1) * (m_nv + 1)) / (8 * M_PI * dotHV) * std::pow(std::fabs(halfv.z), exp);
        SampledSpectrum fs = commonTerm / std::fmax(std::fabs(query.dir_sn.z), std::fabs(dir.z)) * F;
        if (rev_fs)
            *rev_fs = fs;
        return fs;
    }
    
    float AshikhminSpecularBRDF::evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, float* revPDF) const {
        Vector3D halfv = halfVector(query.dir_sn, dir);
        float dotHV = dot(halfv, query.dir_sn);
        float exp = (m_nu * halfv.x * halfv.x + m_nv * halfv.y * halfv.y) / (1 - halfv.z * halfv.z);
        float ret = std::sqrt((m_nu + 1) * (m_nv + 1)) / (8 * M_PI * dotHV) * std::pow(std::fabs(halfv.z), exp);
        if (revPDF)
            *revPDF = ret;
        return ret;
    }
    
    float AshikhminSpecularBRDF::weightInternal(const BSDFQuery &query, const BSDFSample &smp) const {
#ifdef Use_BSDF_Actual_Weights
        BSDFQueryResult result;
        float fs = sample(query, smp, &result).maxValue();
        return result.dirPDF > 0.0f ? fs * std::fabs(result.dir_sn.z) / result.dirPDF : 0.0f;
#else
        SampledSpectrum F = m_Rs + (SampledSpectrum::One - m_Rs) * std::pow(1.0f - std::fabs(query.dir_sn.z), 5);
        return F.maxValue();
        //    return m_Rs.maxValue();
#endif
    }
    
    float AshikhminSpecularBRDF::weightInternal(const BSDFQuery &query, const Vector3D &dir, float* revWeight) const {
#ifdef Use_BSDF_Actual_Weights
        if (revWeight) {
            SampledSpectrum rev_fs;
            float revDirPDF;
            float fs = evaluate(query, dir, &rev_fs).maxValue();
            float dirPDF = evaluatePDF(query, dir, &revDirPDF);
            if (dirPDF > 0.0f) {
                *revWeight = rev_fs.maxValue() * std::fabs(query.dir_sn.z) / revDirPDF;
                return fs * std::fabs(dir.z) / dirPDF;
            }
            else {
                *revWeight = 0.0f;
                return 0.0f;
            }
        }
        else {
            BSDFQueryResult result;
            float fs = evaluate(query, dir).maxValue();
            float dirPDF = evaluatePDF(query, dir);
            return dirPDF > 0 ? fs * std::fabs(dir.z) / dirPDF : 0.0f;
        }
#else
        return weight(query, BSDFSample(0, 0, 0));
#endif
    }
    
    SampledSpectrum AshikhminSpecularBRDF::getBaseColorInternal(DirectionType flags) const {
        return m_Rs;
    }
    
    
    SampledSpectrum AshikhminDiffuseBRDF::sampleInternal(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult *result) const {
        result->dir_sn = cosineSampleHemisphere(smp.uDir[0], smp.uDir[1]);
        result->dirPDF = result->dir_sn.z / M_PI;
        result->dirType = m_type;
        SampledSpectrum fs = evaluateInternal(query, result->dir_sn, nullptr);
        if (result->reverse) {
            result->reverse->fs = fs;
            result->reverse->dirPDF = std::fabs(query.dir_sn.z) / M_PI;
        }
        return fs;
    }
    
    SampledSpectrum AshikhminDiffuseBRDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const {
        SampledSpectrum fs = (28 * m_Rd / (23 * M_PI) * (SampledSpectrum::One - m_Rs) *
                              (1.0f - std::pow(1.0f - std::fabs(query.dir_sn.z) / 2, 5)) *
                              (1.0f - std::pow(1.0f - std::fabs(dir.z) / 2, 5))
                              );
        if (rev_fs)
            *rev_fs = fs;
        return fs;
    }
    
    float AshikhminDiffuseBRDF::evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, float* revPDF) const {
        if (query.dir_sn.z * dir.z >= 0.0f) {
            if (revPDF)
                *revPDF = std::fabs(query.dir_sn.z) / M_PI;
            return std::fabs(dir.z) / M_PI;
        }
        else {
            if (revPDF)
                *revPDF = 0.0f;
            return 0.0f;
        }
    }
    
    float AshikhminDiffuseBRDF::weightInternal(const BSDFQuery &query, const BSDFSample &smp) const {
#ifdef Use_BSDF_Actual_Weights
        BSDFQueryResult result;
        float fs = sample(query, smp, &result).maxValue();
        return result.dirPDF > 0.0f ? fs * std::fabs(result.dir_sn.z) / result.dirPDF : 0.0f;
#else
        float F = std::pow((1.0f - std::pow(1.0f - std::fabs(query.dir_sn.z) / 2, 5)), 2);
        return (28 * m_Rd / (23 * M_PI) * (SampledSpectrum::One - m_Rs) * F).maxValue();
        //    return m_Rd.maxValue();
#endif
    }
    
    float AshikhminDiffuseBRDF::weightInternal(const BSDFQuery &query, const Vector3D &dir, float* revWeight) const {
#ifdef Use_BSDF_Actual_Weights
        BSDFQueryResult result;
        if (revWeight) {
            SampledSpectrum rev_fs;
            float revDirPDF;
            float fs = evaluate(query, dir, &rev_fs).maxValue();
            float dirPDF = evaluatePDF(query, dir, &revDirPDF);
            if (dirPDF > 0.0f) {
                *revWeight = rev_fs.maxValue() * std::fabs(query.dir_sn.z) / revDirPDF;
                return fs * std::fabs(dir.z) / dirPDF;
            }
            else {
                *revWeight = 0.0f;
                return 0.0f;
            }
        }
        else {
            float fs = evaluate(query, dir).maxValue();
            float dirPDF = evaluatePDF(query, dir);
            return dirPDF > 0.0f ? fs * std::fabs(dir.z) / dirPDF : 0.0f;
        }
#else
        return weight(query, BSDFSample(0, 0, 0));
#endif
    }
    
    SampledSpectrum AshikhminDiffuseBRDF::getBaseColorInternal(DirectionType flags) const {
        return m_Rd;
    }    
}

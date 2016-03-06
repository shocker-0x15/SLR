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
        if (result->dir_sn.z * query.dir_sn.z <= 0) {
            result->dirPDF = 0.0f;
            return SampledSpectrum::Zero;
        }
        
        float dotHV = dot(halfv, query.dir_sn);
        float exp = (m_nu * halfv.x * halfv.x + m_nv * halfv.y * halfv.y) / (1 - halfv.z * halfv.z);
        SampledSpectrum F = m_Rs + (SampledSpectrum::One - m_Rs) * std::pow(1.0f - dotHV, 5);
        float commonTerm = std::sqrt((m_nu + 1) * (m_nv + 1)) / (8 * M_PI * dotHV) * std::pow(std::fabs(halfv.z), exp);
        
        result->dirPDF = commonTerm;
        result->dirType = m_type;
        SampledSpectrum fs = commonTerm / std::fmax(std::fabs(query.dir_sn.z), std::fabs(result->dir_sn.z)) * F;
        if (result->reverse) {
            result->reverse->fs = fs;
            result->reverse->dirPDF = commonTerm;
        }
        return fs;
    }
    
    SampledSpectrum AshikhminSpecularBRDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const {
        if (dir.z * query.dir_sn.z <= 0) {
            if (rev_fs)
                *rev_fs = SampledSpectrum::Zero;
            return SampledSpectrum::Zero;
        }
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
    
    SampledSpectrum AshikhminSpecularBRDF::evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* revPDF) const {
        if (dir.z * query.dir_sn.z <= 0) {
            if (revPDF)
                *revPDF = 0.0f;
            return 0.0f;
        }
        Vector3D halfv = halfVector(query.dir_sn, dir);
        float dotHV = dot(halfv, query.dir_sn);
        float exp = (m_nu * halfv.x * halfv.x + m_nv * halfv.y * halfv.y) / (1 - halfv.z * halfv.z);
        float ret = std::sqrt((m_nu + 1) * (m_nv + 1)) / (8 * M_PI * dotHV) * std::pow(std::fabs(halfv.z), exp);
        if (revPDF)
            *revPDF = ret;
        return ret;
    }
    
    SampledSpectrum AshikhminSpecularBRDF::weightInternal(const SLR::BSDFQuery &query) const {
        float dotHV = std::fabs(query.dir_sn.z);
        return m_Rs + (SampledSpectrum::One - m_Rs) * std::pow(1.0f - dotHV, 5);
    }
    
    SampledSpectrum AshikhminSpecularBRDF::getBaseColorInternal(DirectionType flags) const {
        return m_Rs;
    }
    
    
    SampledSpectrum AshikhminDiffuseBRDF::sampleInternal(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult *result) const {
        result->dir_sn = cosineSampleHemisphere(smp.uDir[0], smp.uDir[1]);
        result->dirPDF = result->dir_sn.z / M_PI;
        result->dirType = m_type;
        SampledSpectrum fs = evaluateInternal(query, result->dir_sn, nullptr);
        result->dir_sn.z *= dot(query.dir_sn, query.gNormal_sn) > 0 ? 1 : -1;
        if (result->reverse) {
            result->reverse->fs = fs;
            result->reverse->dirPDF = std::fabs(query.dir_sn.z) / M_PI;
        }
        return fs;
    }
    
    SampledSpectrum AshikhminDiffuseBRDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const {
        if (dir.z * query.dir_sn.z <= 0) {
            if (rev_fs)
                *rev_fs = SampledSpectrum::Zero;
            return SampledSpectrum::Zero;
        }
        SampledSpectrum fs = (28 * m_Rd / (23 * M_PI) * (SampledSpectrum::One - m_Rs) *
                              (1.0f - std::pow(1.0f - std::fabs(query.dir_sn.z) / 2, 5)) *
                              (1.0f - std::pow(1.0f - std::fabs(dir.z) / 2, 5))
                              );
        if (rev_fs)
            *rev_fs = fs;
        return fs;
    }
    
    SampledSpectrum AshikhminDiffuseBRDF::evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* revPDF) const {
        if (query.dir_sn.z * dir.z <= 0.0f) {
            if (revPDF)
                *revPDF = 0.0f;
            return 0.0f;
        }
        if (revPDF)
            *revPDF = std::fabs(query.dir_sn.z) / M_PI;
        return std::fabs(dir.z) / M_PI;
    }
    
    SampledSpectrum AshikhminDiffuseBRDF::weightInternal(const SLR::BSDFQuery &query) const {
        return (2954.0f / 1863.0f * m_Rd * (SampledSpectrum::One - m_Rs) * (1.0f - std::pow(1.0f - std::fabs(query.dir_sn.z) / 2, 5)));
    }
    
    SampledSpectrum AshikhminDiffuseBRDF::getBaseColorInternal(DirectionType flags) const {
        return m_Rd;
    }    
}

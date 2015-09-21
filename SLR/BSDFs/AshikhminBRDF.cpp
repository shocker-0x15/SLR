//
//  AshikhminBRDF.cpp
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "AshikhminBRDF.h"
#include "../Core/distributions.h"

SampledSpectrum AshikhminSpecularBRDF::sample(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult *result) const {
    if (!matches(query.flags, result))
        return SampledSpectrum::Zero;
    
    float quad = 2 * M_PI * smp.uDir[1];
    float phi_h = std::atan2(std::sqrt(m_nu + 1) * std::sin(quad), std::sqrt(m_nv + 1) * std::cos(quad));
    float cosphi = std::cos(phi_h);
    float sinphi = std::sin(phi_h);
    float theta_h = std::acos(std::pow(1 - smp.uDir[0], 1.0f / (m_nu * cosphi * cosphi + m_nv * sinphi * sinphi + 1)));
    if (query.dir_sn.z < 0)
        theta_h = M_PI - theta_h;
    Vector3D halfv = Vector3D(std::sin(theta_h) * std::cos(phi_h), std::sin(theta_h) * std::sin(phi_h), std::cos(theta_h));
    result->dir_sn = 2 * dot(query.dir_sn, halfv) * halfv - query.dir_sn;
    result->dirPDF = evaluatePDF(query, result->dir_sn);
    result->dirType = m_type;
    return evaluateInternal(query, result->dir_sn);
}

SampledSpectrum AshikhminSpecularBRDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir) const {
    if (!query.flags.matches(m_type))
        return SampledSpectrum::Zero;
    
    Vector3D halfv = halfVector(query.dir_sn, dir);
    float dotHV = dot(halfv, query.dir_sn);
    float exp = (m_nu * halfv.x * halfv.x + m_nv * halfv.y * halfv.y) / (1 - halfv.z * halfv.z);
    SampledSpectrum fr = m_Rs + (SampledSpectrum::One - m_Rs) * std::pow(1.0f - dotHV, 5);
    SampledSpectrum ret = std::sqrt((m_nu + 1) * (m_nv + 1)) / (8 * M_PI) *
    std::pow(std::fabs(halfv.z), exp) / (dotHV * std::fmax(std::fabs(query.dir_sn.z), std::fabs(dir.z))) * fr;
    BSDF_EVALUATE_ASSERT;
    return ret;
}

float AshikhminSpecularBRDF::evaluatePDF(const BSDFQuery &query, const Vector3D &dir) const {
    if (!query.flags.matches(m_type))
        return 0;
    Vector3D halfv = halfVector(query.dir_sn, dir);
    float exp = (m_nu * halfv.x * halfv.x + m_nv * halfv.y * halfv.y) / (1 - halfv.z * halfv.z);
    float ret = std::sqrt((m_nu + 1) * (m_nv + 1)) / (2 * M_PI) * std::pow(std::fabs(halfv.z), exp) / (4 * dot(query.dir_sn, halfv));
    BSDF_EVALUATE_PDF_ASSERT;
    return ret;
}

float AshikhminSpecularBRDF::weight(const BSDFQuery &query, const BSDFSample &smp) const {
    if (!query.flags.matches(m_type))
        return 0;
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

float AshikhminSpecularBRDF::weight(const BSDFQuery &query, const Vector3D &dir) const {
    if (!query.flags.matches(m_type))
        return 0;
#ifdef Use_BSDF_Actual_Weights
    BSDFQueryResult result;
    float fs = evaluate(query, dir).maxValue();
    float dirPDF = evaluatePDF(query, dir);
    return dirPDF > 0 ? fs * std::fabs(dir.z) / dirPDF : 0.0f;
#else
    return weight(query, BSDFSample(0, 0, 0));
#endif
}

SampledSpectrum AshikhminSpecularBRDF::getBaseColor(DirectionType flags) const {
    if (!flags.matches(m_type))
        return SampledSpectrum::Zero;
    return m_Rs;
}


SampledSpectrum AshikhminDiffuseBRDF::sample(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult *result) const {
    if (!matches(query.flags, result))
        return SampledSpectrum::Zero;
    
    result->dir_sn = cosineSampleHemisphere(smp.uDir[0], smp.uDir[1]);
    result->dirPDF = result->dir_sn.z / M_PI;
    result->dirType = m_type;
    return evaluateInternal(query, result->dir_sn);
}

SampledSpectrum AshikhminDiffuseBRDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir) const {
    if (!query.flags.matches(m_type))
        return SampledSpectrum::Zero;
    
    SampledSpectrum ret = (28 * m_Rd / (23 * M_PI) * (SampledSpectrum::One - m_Rs) *
                    (1.0f - std::pow(1.0f - std::fabs(query.dir_sn.z) / 2, 5)) *
                    (1.0f - std::pow(1.0f - std::fabs(dir.z) / 2, 5))
                    );
    BSDF_EVALUATE_ASSERT;
    return ret;
}

float AshikhminDiffuseBRDF::evaluatePDF(const BSDFQuery &query, const Vector3D &dir) const {
    if (!query.flags.matches(m_type))
        return 0;
    return query.dir_sn.z * dir.z >= 0.0f ? dir.z / M_PI : 0.0f;
}

float AshikhminDiffuseBRDF::weight(const BSDFQuery &query, const BSDFSample &smp) const {
    if (!query.flags.matches(m_type))
        return 0;
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

float AshikhminDiffuseBRDF::weight(const BSDFQuery &query, const Vector3D &dir) const {
    if (!query.flags.matches(m_type))
        return 0;
#ifdef Use_BSDF_Actual_Weights
    BSDFQueryResult result;
    float fs = evaluate(query, dir).maxValue();
    float dirPDF = evaluatePDF(query, dir);
    return dirPDF > 0.0f ? fs * std::fabs(dir.z) / dirPDF : 0.0f;
#else
    return weight(query, BSDFSample(0, 0, 0));
#endif
}

SampledSpectrum AshikhminDiffuseBRDF::getBaseColor(DirectionType flags) const {
    if (!flags.matches(m_type))
        return SampledSpectrum::Zero;
    return m_Rd;
}

//
//  basic_BSDFs.cpp
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "basic_BSDFs.h"
#include "../Core/distributions.h"

SampledSpectrum LambertianBRDF::sample(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const {
    if (!matches(query.flags, result))
        return SampledSpectrum::Zero;
    result->dir_sn = cosineSampleHemisphere(smp.uDir[0], smp.uDir[1]);
    result->dirPDF = result->dir_sn.z / M_PI;
    result->dirType = m_type;
    return m_R / M_PI;
}

SampledSpectrum LambertianBRDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir) const {
    if (!query.flags.matches(m_type))
        return SampledSpectrum::Zero;
    return m_R / M_PI;
}

float LambertianBRDF::evaluatePDF(const BSDFQuery &query, const Vector3D &dir) const {
    if (!query.flags.matches(m_type))
        return 0;
    return query.dir_sn.z * dir.z >= 0.0f ? std::abs(dir.z) / M_PI : 0.0f;
}

float LambertianBRDF::weight(const BSDFQuery &query, const BSDFSample &smp) const {
    if (!query.flags.matches(m_type))
        return 0;
#ifdef Use_BSDF_Actual_Weights
    return m_R.maxValue();
#else
    return m_R.maxValue();
#endif
}

float LambertianBRDF::weight(const BSDFQuery &query, const Vector3D &dir) const {
    if (!query.flags.matches(m_type))
        return 0;
#ifdef Use_BSDF_Actual_Weights
    return m_R.maxValue();
#else
    return m_R.maxValue();
#endif
}

SampledSpectrum LambertianBRDF::getBaseColor(DirectionType flags) const {
    if (!flags.matches(m_type))
        return SampledSpectrum::Zero;
    return m_R;
}

SampledSpectrum SpecularBRDF::sample(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const {
    if (!matches(query.flags, result))
        return SampledSpectrum::Zero;
    result->dir_sn = Vector3D(-query.dir_sn.x, -query.dir_sn.y, query.dir_sn.z);
    result->dirPDF = 1.0f;
    result->dirType = m_type;
    SampledSpectrum ret = m_coeffR * m_fresnel->evaluate(query.dir_sn.z) / std::fabs(query.dir_sn.z);
    BSDF_SAMPLE_ASSERT;
    return ret;
}

SampledSpectrum SpecularBRDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir) const {
    if (!query.flags.matches(m_type))
        return SampledSpectrum::Zero;
    return SampledSpectrum::Zero;
}

float SpecularBRDF::evaluatePDF(const BSDFQuery &query, const Vector3D &dir) const {
    if (!query.flags.matches(m_type))
        return 0.0f;
    return 0.0f;
}

float SpecularBRDF::weight(const BSDFQuery &query, const BSDFSample &smp) const {
    if (!query.flags.matches(m_type))
        return 0.0f;
#ifdef Use_BSDF_Actual_Weights
    BSDFQueryResult result;
    float fs = sample(query, smp, &result).maxValue();
    return fs * std::fabs(result.dir_sn.z) / result.dirPDF;
#else
    return m_coeffR.maxValue();
#endif
}

float SpecularBRDF::weight(const BSDFQuery &query, const Vector3D &dir) const {
    return 0.0f;
}

SampledSpectrum SpecularBRDF::getBaseColor(DirectionType flags) const {
    if (!flags.matches(m_type))
        return SampledSpectrum::Zero;
    return m_coeffR;
}

SampledSpectrum SpecularBTDF::sample(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const {
    if (!matches(query.flags, result))
        return SampledSpectrum::Zero;
    
    bool entering = query.dir_sn.z > 0.0f;
    float eEnter = entering ? m_fresnel.etaExt()[query.wlHint] : m_fresnel.etaInt()[query.wlHint];
    float eExit = entering ? m_fresnel.etaInt()[query.wlHint] : m_fresnel.etaExt()[query.wlHint];
    
    float sinEnter2 = 1.0f - query.dir_sn.z * query.dir_sn.z;
    float rrEta = eEnter / eExit;// reciprocal of relative IOR.
    float sinExit2 = rrEta * rrEta * sinEnter2;
    
    if (sinExit2 >= 1.0f) {
        result->dirPDF = 0.0f;
        return SampledSpectrum::Zero;
    }
    float cosExit = std::sqrt(std::fmax(0.0f, 1.0f - sinExit2));
    if (entering)
        cosExit = -cosExit;
    result->dir_sn = Vector3D(rrEta * -query.dir_sn.x, rrEta * -query.dir_sn.y, cosExit);
    result->dirPDF = 1.0f;
    result->dirType = m_type | (m_dispersive ? DirectionType::Dispersive : DirectionType());
    cosExit = std::fabs(cosExit);
    float F = FresnelDielectric::evalF(eEnter, eExit, std::fabs(query.dir_sn.z), cosExit);
    SampledSpectrum ret = SampledSpectrum::Zero;
    ret[query.wlHint] = m_coeffT[query.wlHint] * (1.0f - F) / cosExit;
    BSDF_SAMPLE_ASSERT;
    return ret;
}

SampledSpectrum SpecularBTDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir) const {
    if (!query.flags.matches(m_type))
        return SampledSpectrum::Zero;
    return SampledSpectrum::Zero;
}

float SpecularBTDF::evaluatePDF(const BSDFQuery &query, const Vector3D &dir) const {
    if (!query.flags.matches(m_type))
        return 0.0f;
    return 0.0f;
}

float SpecularBTDF::weight(const BSDFQuery &query, const BSDFSample &smp) const {
    if (!query.flags.matches(m_type))
        return 0.0f;
#ifdef Use_BSDF_Actual_Weights
    BSDFQueryResult result;
    float fs = sample(query, smp, &result)[query.wlHint];
    return result.dirPDF > 0.0f ? fs * std::fabs(result.dir_sn.z) / result.dirPDF : 0.0f;
#else
    return m_coeffT.maxValue();
#endif
}

float SpecularBTDF::weight(const BSDFQuery &query, const Vector3D &dir) const {
    return 0.0f;
}

SampledSpectrum SpecularBTDF::getBaseColor(DirectionType flags) const {
    if (!flags.matches(m_type))
        return SampledSpectrum::Zero;
    return m_coeffT;
}

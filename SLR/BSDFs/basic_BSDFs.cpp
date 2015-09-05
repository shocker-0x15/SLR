//
//  basic_BSDFs.cpp
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "basic_BSDFs.h"
#include "../Core/distributions.h"

Spectrum LambertianBRDF::sample(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const {
    if (!matches(query.flags, result))
        return Spectrum::Zero;
    result->dir_sn = cosineSampleHemisphere(smp.uDir[0], smp.uDir[1]);
    result->dirPDF = result->dir_sn.z / M_PI;
    result->dirType = m_type;
    return m_R / M_PI;
}

Spectrum LambertianBRDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir) const {
    if (!query.flags.matches(m_type))
        return Spectrum::Zero;
    return m_R / M_PI;
}

float LambertianBRDF::evaluatePDF(const BSDFQuery &query, const Vector3D &dir) const {
    if (!query.flags.matches(m_type))
        return 0;
    return query.dir_sn.z * dir.z >= 0.0f ? std::abs(dir.z) / M_PI : 0.0f;
}

float LambertianBRDF::weight(const BSDFQuery &query) const {
    if (!query.flags.matches(m_type))
        return 0;
    return m_R.maxValue();
}

Spectrum LambertianBRDF::getBaseColor(DirectionType flags) const {
    if (!flags.matches(m_type))
        return Spectrum::Zero;
    return m_R;
}

Spectrum SpecularBRDF::sample(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const {
    if (!matches(query.flags, result))
        return Spectrum::Zero;
    result->dir_sn = Vector3D(-query.dir_sn.x, -query.dir_sn.y, query.dir_sn.z);
    result->dirPDF = 1.0f;
    result->dirType = m_type;
    Spectrum ret = m_coeffR * m_fresnel->evaluate(query.dir_sn.z) / std::fabs(query.dir_sn.z);
    BSDF_SAMPLE_ASSERT;
    return ret;
}

Spectrum SpecularBRDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir) const {
    if (!query.flags.matches(m_type))
        return Spectrum::Zero;
    return Spectrum::Zero;
}

float SpecularBRDF::evaluatePDF(const BSDFQuery &query, const Vector3D &dir) const {
    if (!query.flags.matches(m_type))
        return 0.0f;
    return 0.0f;
}

float SpecularBRDF::weight(const BSDFQuery &query) const {
    if (!query.flags.matches(m_type))
        return 0.0f;
    return m_coeffR.maxValue();
}

Spectrum SpecularBRDF::getBaseColor(DirectionType flags) const {
    if (!flags.matches(m_type))
        return Spectrum::Zero;
    return m_coeffR;
}

Spectrum SpecularBTDF::sample(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const {
    if (!matches(query.flags, result))
        return Spectrum::Zero;
    
    bool entering = query.dir_sn.z > 0.0f;
    float eEnter = m_fresnel.etaExt();
    float eExit = m_fresnel.etaInt();
    if (!entering)
        std::swap(eEnter, eExit);
    
    float sinEnter2 = 1.0f - query.dir_sn.z * query.dir_sn.z;
    float rrEta = eEnter / eExit;// reciprocal of relative IOR.
    float sinExit2 = rrEta * rrEta * sinEnter2;
    
    if (sinExit2 >= 1.0f) {
        result->dirPDF = 0.0f;
        return Spectrum::Zero;
    }
    float cosExit = std::sqrt(std::fmax(0.0f, 1.0f - sinExit2));
    if (entering)
        cosExit = -cosExit;
    result->dir_sn = Vector3D(rrEta * -query.dir_sn.x, rrEta * -query.dir_sn.y, cosExit);
    result->dirPDF = 1.0f;
    result->dirType = m_type;
    Spectrum F = Spectrum::One - m_fresnel.evaluate(query.dir_sn.z);
    Spectrum ret = m_coeffT * F / std::fabs(query.dir_sn.z);
    BSDF_SAMPLE_ASSERT;
    return ret;
}

Spectrum SpecularBTDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir) const {
    if (!query.flags.matches(m_type))
        return Spectrum::Zero;
    return Spectrum::Zero;
}

float SpecularBTDF::evaluatePDF(const BSDFQuery &query, const Vector3D &dir) const {
    if (!query.flags.matches(m_type))
        return 0.0f;
    return 0.0f;
}

float SpecularBTDF::weight(const BSDFQuery &query) const {
    if (!query.flags.matches(m_type))
        return 0.0f;
    return m_coeffT.maxValue();
}

Spectrum SpecularBTDF::getBaseColor(DirectionType flags) const {
    if (!flags.matches(m_type))
        return Spectrum::Zero;
    return m_coeffT;
}

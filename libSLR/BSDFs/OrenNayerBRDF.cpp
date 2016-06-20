//
//  OrenNayerBRDF.cpp
//  SLR
//
//  Created by 渡部 心 on 2016/06/17.
//  Copyright © 2016年 渡部 心. All rights reserved.
//

#include "OrenNayerBRDF.h"
#include "../Core/distributions.h"

namespace SLR {
    SampledSpectrum OrenNayerBRDF::sampleInternal(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult* result) const {
        bool frontSide = dot(query.dir_sn, query.gNormal_sn) > 0;
        result->dir_sn = cosineSampleHemisphere(uDir[0], uDir[1]);
        result->dirPDF = result->dir_sn.z / M_PI;
        result->dirType = m_type;
        result->dir_sn.z *= frontSide ? 1 : -1;
        
        float sinThetaI = 1.0f - result->dir_sn.z * result->dir_sn.z;
        float sinThetaO = 1.0f - query.dir_sn.z * query.dir_sn.z;
        float absTanThetaI = sinThetaI / std::abs(result->dir_sn.z);
        float absTanThetaO = sinThetaO / std::abs(query.dir_sn.z);
        float sinAlpha = std::max(sinThetaI, sinThetaO);
        float tanBeta = std::min(absTanThetaI, absTanThetaO);
        float cos_dAzimuth = (result->dir_sn.x * query.dir_sn.x + result->dir_sn.y * query.dir_sn.y) / (sinThetaI * sinThetaO);
        if (!std::isfinite(cos_dAzimuth))
            cos_dAzimuth = 0.0f;
        SampledSpectrum fs = m_R * ((m_A + m_B * std::max(0.0f, cos_dAzimuth) * sinAlpha * tanBeta) / M_PI);
        if (result->reverse) {
            result->reverse->fs = fs;
            result->reverse->dirPDF = std::fabs(query.dir_sn.z) / M_PI;
        }
        return fs;
    }
    
    SampledSpectrum OrenNayerBRDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const {
        if (query.dir_sn.z * dir.z <= 0.0f) {
            SampledSpectrum fs = SampledSpectrum::Zero;
            if (rev_fs)
                *rev_fs = fs;
            return fs;
        }
        float sinThetaI = 1.0f - dir.z * dir.z;
        float sinThetaO = 1.0f - query.dir_sn.z * query.dir_sn.z;
        float absTanThetaI = sinThetaI / std::abs(dir.z);
        float absTanThetaO = sinThetaO / std::abs(query.dir_sn.z);
        float sinAlpha = std::max(sinThetaI, sinThetaO);
        float tanBeta = std::min(absTanThetaI, absTanThetaO);
        float cos_dAzimuth = (dir.x * query.dir_sn.x + dir.y * query.dir_sn.y) / (sinThetaI * sinThetaO);
        SampledSpectrum fs = m_R * ((m_A + m_B * std::max(0.0f, cos_dAzimuth) * sinAlpha * tanBeta) / M_PI);
        if (rev_fs)
            *rev_fs = fs;
        return fs;
    }
    
    SampledSpectrum OrenNayerBRDF::evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* revPDF) const {
        if (query.dir_sn.z * dir.z <= 0.0f) {
            if (revPDF)
                *revPDF = 0.0f;
            return 0.0f;
        }
        if (revPDF)
            *revPDF = std::fabs(query.dir_sn.z) / M_PI;
        return std::abs(dir.z) / M_PI;
    }
    
    SampledSpectrum OrenNayerBRDF::weightInternal(const SLR::BSDFQuery &query) const {
        return m_R;
    }
    
    SampledSpectrum OrenNayerBRDF::getBaseColorInternal(DirectionType flags) const {
        return m_R;
    }
}
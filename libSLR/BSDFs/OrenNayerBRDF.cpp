//
//  OrenNayerBRDF.cpp
//
//  Created by 渡部 心 on 2016/06/17.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "OrenNayerBRDF.h"
#include "../Core/distributions.h"

namespace SLR {
    SampledSpectrum OrenNayerBRDF::sampleInternal(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult* result) const {
        bool frontSide = dot(query.dirLocal, query.gNormalLocal) > 0;
        result->dirLocal = cosineSampleHemisphere(uDir[0], uDir[1]);
        result->dirPDF = result->dirLocal.z / M_PI;
        result->sampledType = m_type;
        result->dirLocal.z *= frontSide ? 1 : -1;
        
        float sinThetaI = 1.0f - result->dirLocal.z * result->dirLocal.z;
        float sinThetaO = 1.0f - query.dirLocal.z * query.dirLocal.z;
        float absTanThetaI = sinThetaI / std::abs(result->dirLocal.z);
        float absTanThetaO = sinThetaO / std::abs(query.dirLocal.z);
        float sinAlpha = std::max(sinThetaI, sinThetaO);
        float tanBeta = std::min(absTanThetaI, absTanThetaO);
        float cos_dAzimuth = (result->dirLocal.x * query.dirLocal.x + result->dirLocal.y * query.dirLocal.y) / (sinThetaI * sinThetaO);
        if (!std::isfinite(cos_dAzimuth))
            cos_dAzimuth = 0.0f;
        SampledSpectrum fs = m_R * ((m_A + m_B * std::max(0.0f, cos_dAzimuth) * sinAlpha * tanBeta) / M_PI);
        if (result->reverse) {
            result->reverse->fs = fs;
            result->reverse->dirPDF = std::fabs(query.dirLocal.z) / M_PI;
        }
        return fs;
    }
    
    SampledSpectrum OrenNayerBRDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const {
        if (query.dirLocal.z * dir.z <= 0.0f) {
            SampledSpectrum fs = SampledSpectrum::Zero;
            if (rev_fs)
                *rev_fs = fs;
            return fs;
        }
        float sinThetaI = 1.0f - dir.z * dir.z;
        float sinThetaO = 1.0f - query.dirLocal.z * query.dirLocal.z;
        float absTanThetaI = sinThetaI / std::abs(dir.z);
        float absTanThetaO = sinThetaO / std::abs(query.dirLocal.z);
        float sinAlpha = std::max(sinThetaI, sinThetaO);
        float tanBeta = std::min(absTanThetaI, absTanThetaO);
        float cos_dAzimuth = (dir.x * query.dirLocal.x + dir.y * query.dirLocal.y) / (sinThetaI * sinThetaO);
        SampledSpectrum fs = m_R * ((m_A + m_B * std::max(0.0f, cos_dAzimuth) * sinAlpha * tanBeta) / M_PI);
        if (rev_fs)
            *rev_fs = fs;
        return fs;
    }
    
    float OrenNayerBRDF::evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, float* revPDF) const {
        if (query.dirLocal.z * dir.z <= 0.0f) {
            if (revPDF)
                *revPDF = 0.0f;
            return 0.0f;
        }
        if (revPDF)
            *revPDF = std::fabs(query.dirLocal.z) / M_PI;
        return std::abs(dir.z) / M_PI;
    }
    
    float OrenNayerBRDF::weightInternal(const BSDFQuery &query) const {
        return m_R.luminance();
    }
    
    SampledSpectrum OrenNayerBRDF::getBaseColorInternal(DirectionType flags) const {
        return m_R;
    }
}

//
//  directional_distribution_functions.cpp
//
//  Created by 渡部 心 on 2015/06/28.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "directional_distribution_functions.h"
#include "distributions.h"

namespace SLR {
    SLR_API const DirectionType DirectionType::LowFreq = DirectionType::IE_LowFreq;
    SLR_API const DirectionType DirectionType::HighFreq = DirectionType::IE_HighFreq;
    SLR_API const DirectionType DirectionType::Delta0D = DirectionType::IE_Delta0D;
    SLR_API const DirectionType DirectionType::Delta1D = DirectionType::IE_Delta1D;
    SLR_API const DirectionType DirectionType::NonDelta = DirectionType::IE_NonDelta;
    SLR_API const DirectionType DirectionType::Delta = DirectionType::IE_Delta;
    SLR_API const DirectionType DirectionType::AllFreq = DirectionType::IE_AllFreq;
    SLR_API const DirectionType DirectionType::Reflection = DirectionType::IE_Reflection;
    SLR_API const DirectionType DirectionType::Transmission = DirectionType::IE_Transmission;
    SLR_API const DirectionType DirectionType::Emission = DirectionType::IE_Emission;
    SLR_API const DirectionType DirectionType::Acquisition = DirectionType::IE_Acquisition;
    SLR_API const DirectionType DirectionType::WholeSphere = DirectionType::IE_WholeSphere;
    SLR_API const DirectionType DirectionType::All = DirectionType::IE_All;
    SLR_API const DirectionType DirectionType::Dispersive = DirectionType::IE_Dispersive;
    SLR_API const DirectionType DirectionType::LowFreqReflection = DirectionType::IE_LowFreqReflection;
    SLR_API const DirectionType DirectionType::LowFreqTransmission = DirectionType::IE_LowFreqTransmission;
    SLR_API const DirectionType DirectionType::LowFreqScattering = DirectionType::IE_LowFreqScattering;
    SLR_API const DirectionType DirectionType::HighFreqReflection = DirectionType::IE_HighFreqReflection;
    SLR_API const DirectionType DirectionType::HighFreqTransmission = DirectionType::IE_HighFreqTransmission;
    SLR_API const DirectionType DirectionType::HighFreqScattering = DirectionType::IE_HighFreqScattering;
    SLR_API const DirectionType DirectionType::Delta0DReflection = DirectionType::IE_Delta0DReflection;
    SLR_API const DirectionType DirectionType::Delta0DTransmission = DirectionType::IE_Delta0DTransmission;
    SLR_API const DirectionType DirectionType::Delta0DScattering = DirectionType::IE_Delta0DScattering;
    
    
    SampledSpectrum BSDF::rho(uint32_t numSamples, BSDFSample* samples, float* uDir0, float* uDir1, float* uWl, DirectionType flags, bool fromUpper) const {
        SampledSpectrumSum ret(SampledSpectrum::Zero);
        for (int i = 0; i < numSamples; ++i) {
            Vector3D dir_sn = cosineSampleHemisphere(uDir0[i], uDir1[i]);
            dir_sn.z *= fromUpper ? 1 : -1;
            float dirPDF = dir_sn.z / M_PI;
            
            int16_t wlIdx = std::min(int16_t(SampledSpectrum::NumComponents * uWl[i]), int16_t(SampledSpectrum::NumComponents - 1));
            BSDFQuery query{dir_sn, Normal3D(0, 0, 1), wlIdx, flags, false};
            
            SampledSpectrumSum dirSum(SampledSpectrum::Zero);
            for (int j = 0; j < numSamples; ++j) {
                BSDFQueryResult fsResult;
                SampledSpectrum fs = sample(query, samples[j], &fsResult);
                dirSum += fs * std::fabs(fsResult.dir_sn.z) / fsResult.dirPDF;
            }
            ret += dirSum.result * std::fabs(dir_sn.z) / dirPDF;
        }
        return ret.result / M_PI;
    }
    
    
    SampledSpectrum BSSRDF::sample(const SLR::BSSRDFQuery &query, const SLR::BSSRDFSample &smp, SLR::BSSRDFQueryResult *result) const {
        Vector3D s, t, u;
        float uRadius;
        if (smp.uPos[0] < 0.5f) {
            uRadius = smp.uPos[0] / 0.5f;
            s = query.surfPt.shadingFrame.x;
            t = query.surfPt.shadingFrame.y;
            u = query.surfPt.shadingFrame.z;
        }
        else if (smp.uPos[1] < 0.75f) {
            uRadius = (smp.uPos[1] - 0.5f) / 0.25f;
            s = query.surfPt.shadingFrame.y;
            t = query.surfPt.shadingFrame.z;
            u = query.surfPt.shadingFrame.x;
        }
        else {
            uRadius = (smp.uPos[1] - 0.75f) / 0.25f;
            s = query.surfPt.shadingFrame.z;
            t = query.surfPt.shadingFrame.x;
            u = query.surfPt.shadingFrame.y;
        }
        
        return SampledSpectrum::Zero;
    }
    
    SampledSpectrum FresnelNoOp::evaluate(float cosEnter) const {
        return SampledSpectrum::One;
    }
    
    float FresnelNoOp::evaluate(float cosEnter, uint32_t wlIdx) const {
        return 1.0f;
    }
    
    SampledSpectrum FresnelConductor::evaluate(float cosEnter) const {
        cosEnter = std::fabs(cosEnter);
        float cosEnter2 = cosEnter * cosEnter;
        SampledSpectrum _2EtaCosEnter = 2.0f * m_eta * cosEnter;
        SampledSpectrum tmp_f = m_eta * m_eta + m_k * m_k;
        SampledSpectrum tmp = tmp_f * cosEnter2;
        SampledSpectrum Rparl2 = (tmp - _2EtaCosEnter + 1) / (tmp + _2EtaCosEnter + 1);
        SampledSpectrum Rperp2 = (tmp_f - _2EtaCosEnter + cosEnter2) / (tmp_f + _2EtaCosEnter + cosEnter2);
        return (Rparl2 + Rperp2) / 2.0f;
    }
    
    float FresnelConductor::evaluate(float cosEnter, uint32_t wlIdx) const {
        cosEnter = std::fabs(cosEnter);
        float cosEnter2 = cosEnter * cosEnter;
        float _2EtaCosEnter = 2.0f * m_eta[wlIdx] * cosEnter;
        float tmp_f = m_eta[wlIdx] * m_eta[wlIdx] + m_k[wlIdx] * m_k[wlIdx];
        float tmp = tmp_f * cosEnter2;
        float Rparl2 = (tmp - _2EtaCosEnter + 1) / (tmp + _2EtaCosEnter + 1);
        float Rperp2 = (tmp_f - _2EtaCosEnter + cosEnter2) / (tmp_f + _2EtaCosEnter + cosEnter2);
        return (Rparl2 + Rperp2) / 2.0f;
    }
    
    SampledSpectrum FresnelDielectric::evaluate(float cosEnter) const {
        cosEnter = std::clamp(cosEnter, -1.0f, 1.0f);
        
        bool entering = cosEnter > 0.0f;
        const SampledSpectrum &eEnter = entering ? m_etaExt : m_etaInt;
        const SampledSpectrum &eExit = entering ? m_etaInt : m_etaExt;
        
        SampledSpectrum sinExit = eEnter / eExit * std::sqrt(std::fmax(0.0f, 1.0f - cosEnter * cosEnter));
        SampledSpectrum ret = SampledSpectrum::Zero;
        cosEnter = std::fabs(cosEnter);
        for (int i = 0; i < SampledSpectrum::NumComponents; ++i) {
            if (sinExit[i] >= 1.0f) {
                ret[i] = 1.0f;
            }
            else {
                float cosExit = std::sqrt(std::fmax(0.0f, 1.0f - sinExit[i] * sinExit[i]));
                ret[i] = evalF(eEnter[i], eExit[i], cosEnter, cosExit);
            }
        }
        return ret;
        
//        cosEnter = std::fabs(std::clamp(cosEnter, -1.0f, 1.0f));
//        
//        bool entering = cosEnter > 0.0f;
//        const SampledSpectrum &eEnter = entering ? m_etaExt : m_etaInt;
//        const SampledSpectrum &eExit = entering ? m_etaInt : m_etaExt;
//        
//        SampledSpectrum relIOR = eExit / eEnter;
//        SampledSpectrum g2 = relIOR * relIOR - SampledSpectrum::One + cosEnter * cosEnter;
//        
//        SampledSpectrum ret = SampledSpectrum::Zero;
//        for (int i = 0; i < SampledSpectrum::NumComponents; ++i) {
//            if (g2[i] < 0) {
//                ret[i] = 1.0f;
//            }
//            else {
//                float g_minus_c = std::sqrt(g2[i]) - cosEnter;
//                float g_plus_c = std::sqrt(g2[i]) + cosEnter;
//                float termN = cosEnter * g_plus_c - 1;
//                float termD = cosEnter * g_minus_c + 1;
//                
//                ret[i] = 0.5f * (g_minus_c * g_minus_c) / (g_plus_c * g_plus_c) * (1 + (termN * termN) / (termD * termD));
//            }
//        }
//        return ret;
    }
    
    float FresnelDielectric::evaluate(float cosEnter, uint32_t wlIdx) const {
        cosEnter = std::clamp(cosEnter, -1.0f, 1.0f);
        
        bool entering = cosEnter > 0.0f;
        const float &eEnter = entering ? m_etaExt[wlIdx] : m_etaInt[wlIdx];
        const float &eExit = entering ? m_etaInt[wlIdx] : m_etaExt[wlIdx];
        
        float sinExit = eEnter / eExit * std::sqrt(std::fmax(0.0f, 1.0f - cosEnter * cosEnter));
        cosEnter = std::fabs(cosEnter);
        if (sinExit >= 1.0f) {
            return 1.0f;
        }
        else {
            float cosExit = std::sqrt(std::fmax(0.0f, 1.0f - sinExit * sinExit));
            return evalF(eEnter, eExit, cosEnter, cosExit);
        }
    }
    
    float FresnelDielectric::evalF(float etaEnter, float etaExit, float cosEnter, float cosExit) {
        float Rparl = ((etaExit * cosEnter) - (etaEnter * cosExit)) / ((etaExit * cosEnter) + (etaEnter * cosExit));
        float Rperp = ((etaEnter * cosEnter) - (etaExit * cosExit)) / ((etaEnter * cosEnter) + (etaExit * cosExit));
        return (Rparl * Rparl + Rperp * Rperp) / 2.0f;
    }
    
    
    float GGX::sample(float u0, float u1, Normal3D* m, float* normalPDF) const {
        float theta_m = std::atan(m_alpha_g * std::sqrt(u0 / (1 - u0)));
        float phi_m = 2 * M_PI * u1;
        *m = Normal3D(std::sin(theta_m) * std::cos(phi_m), std::sin(theta_m) * std::sin(phi_m), std::cos(theta_m));
        
        //        float chi = m->z > 0 ? 1 : 0;
        float cosTheta_m = m->z;
        float tanTheta_m = std::tan(theta_m);
        float ret = m_alpha_g * m_alpha_g / (M_PI * std::pow(cosTheta_m, 4) * std::pow(m_alpha_g * m_alpha_g + tanTheta_m * tanTheta_m, 2));
        *normalPDF = ret * m->z;
        
        return ret;
    }
    
    float GGX::evaluate(const Normal3D &m) const {
        if (m.z <= 0)
            return 0.0f;
        float theta_m = std::acos(m.z);
        float cosTheta_m = m.z;
        float tanTheta_m = std::tan(theta_m);
        return m_alpha_g * m_alpha_g / (M_PI * std::pow(cosTheta_m, 4) * std::pow(m_alpha_g * m_alpha_g + tanTheta_m * tanTheta_m, 2));
    }
    
    float GGX::evaluatePDF(const Normal3D &m) const {
        if (m.z <= 0)
            return 0.0f;
        float theta_m = std::acos(m.z);
        float cosTheta_m = m.z;
        float tanTheta_m = std::tan(theta_m);
        float D = m_alpha_g * m_alpha_g / (M_PI * std::pow(cosTheta_m, 4) * std::pow(m_alpha_g * m_alpha_g + tanTheta_m * tanTheta_m, 2));
        return D * m.z;
    }
    
    float GGX::evaluateSmithG1(const Vector3D &v, const Normal3D &m) const {
        float chi = (dot(v, m) / std::copysign(v.z, m.z)) > 0 ? 1 : 0;
        float theta_v = std::acos(std::clamp(v.z, -1.0f, 1.0f));
        return chi * 2 / (1 + std::sqrt(1 + std::pow(m_alpha_g * std::tan(theta_v), 2)));
    }
}

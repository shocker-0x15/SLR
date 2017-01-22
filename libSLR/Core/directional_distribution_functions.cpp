//
//  directional_distribution_functions.cpp
//
//  Created by 渡部 心 on 2015/06/28.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "directional_distribution_functions.h"
#include "distributions.h"
#include "../Memory/ArenaAllocator.h"
#include "../Core/light_path_samplers.h"

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
        uint32_t numFails = 0;
        for (int i = 0; i < numSamples; ++i) {
            Vector3D dir_sn = uniformSampleHemisphere(uDir0[i], uDir1[i]);
            dir_sn.z *= fromUpper ? 1 : -1;
            float dirPDF = 1.0f / (2 * M_PI);
            
            int16_t wlIdx = std::min(int16_t(SampledSpectrum::NumComponents * uWl[i]), int16_t(SampledSpectrum::NumComponents - 1));
            BSDFQuery query{dir_sn, Normal3D(0, 0, 1), wlIdx, flags, false};
            
            BSDFQueryResult fsResult;
            SampledSpectrum fs = sample(query, samples[i], &fsResult);
            if (fsResult.dirPDF == 0.0f) {
                ++numFails;
                continue;
            }
            ret += fs * std::fabs(fsResult.dirLocal.z) * std::fabs(dir_sn.z) / (dirPDF * fsResult.dirPDF);
        }
        return ret.result / (M_PI * (numSamples - numFails));
    }
    
    SampledSpectrum BSDF::sample(const ABDFQuery* query, LightPathSampler &sampler, ArenaAllocator &mem, ABDFQueryResult** result) const {
        BSDFQueryResult* concreteResult = mem.create<BSDFQueryResult>();
        // MEMO: *(const BSDFQuery*)query doesn't seem a good way.
        SampledSpectrum ret = sample(*(const BSDFQuery*)query, sampler.getBSDFSample(), concreteResult);
        *result = concreteResult;
        return ret;
    }
    
    SampledSpectrum BSDF::evaluate(const ABDFQuery* query, const Vector3D &dir, SampledSpectrum* rev_fs) const {
        return evaluate(*(const BSDFQuery*)query, dir, rev_fs);
    }
    
    float BSDF::evaluatePDF(const ABDFQuery* query, const Vector3D &dir, float* revPDF) const {
        return evaluatePDF(*(const BSDFQuery*)query, dir, revPDF);
    }
    
    
    
    SampledSpectrum VolumetricBSDF::sample(const ABDFQuery* query, LightPathSampler &sampler, ArenaAllocator &mem, ABDFQueryResult** result) const {
        PFQueryResult concreteResult;
        PFQuery pfQuery(query->dirLocal, query->wlHint, query->dirTypeFilter);
        SampledSpectrum ret = m_pf->sample(pfQuery, sampler.getPFSample(), &concreteResult);
        *result = mem.create<VolumetricBSDFQueryResult>(concreteResult);
        return m_albedo * ret;
    }
    
    SampledSpectrum VolumetricBSDF::evaluate(const ABDFQuery* query, const Vector3D &dir, SampledSpectrum* rev_fs) const {
        PFQuery pfQuery(query->dirLocal, query->wlHint, query->dirTypeFilter);
        SampledSpectrum ret = m_pf->evaluate(pfQuery, dir);
        if (rev_fs)
            *rev_fs *= m_albedo;
        return m_albedo * ret;
    }
    
    float VolumetricBSDF::evaluatePDF(const ABDFQuery* query, const Vector3D &dir, float* revPDF) const {
        PFQuery pfQuery(query->dirLocal, query->wlHint, query->dirTypeFilter);
        float ret = m_pf->evaluatePDF(pfQuery, dir);
        Vector3D dirIn = pfQuery.dirLocal;
        pfQuery.dirLocal = dir;
        if (revPDF)
            *revPDF = m_pf->evaluatePDF(pfQuery, dirIn);
        return ret;
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
        return evaluate(m) * m.z;
    }
    
    // References
    // Importance Sampling Microfacet-Based BSDFs using the Distribution of Visible Normals
    float GGX::sample(const Vector3D &v, float u0, float u1, Normal3D* m, float* normalPDF) const {
        float alpha_gx = m_alpha_g;
        float alpha_gy = m_alpha_g;
        
        // stretch the input vector
        Vector3D sv = Vector3D(alpha_gx * v.x, alpha_gy * v.y, v.z).normalize();
        float theta_sv = std::acos(sv.z);
        float phi_sv = std::atan2(sv.y, sv.x);
        if (sv.z > 0.99999f) {
            theta_sv = 0.0f;
            phi_sv = 0.0f;
        }
        
        // sample slopes
        float slope_x, slope_y;
        if (theta_sv < 0.0001) {
            // special case (normal incidence)
            const float r = std::sqrt(u0 / (1 - u0));
            const float phi = 2 * M_PI * u1;
            slope_x = r * cos(phi);
            slope_y = r * sin(phi);
        }
        else {
            // precomputations
            const float tan_theta_i = tan(theta_sv);
            const float a = 1 / tan_theta_i;
            const float G1 = 2 / (1 + std::sqrt(1.0 + 1.0 / (a * a)));
            
            // sample slope_x
            const float A = 2.0 * u0 / G1 - 1.0;
            const float tmp = 1.0 / (A * A - 1.0);
            const float B = tan_theta_i;
            const float D = std::sqrt(B * B * tmp * tmp - (A * A - B * B) * tmp);
            const float slope_x_1 = B * tmp - D;
            const float slope_x_2 = B * tmp + D;
            slope_x = (A < 0 || slope_x_2 > 1.0 / tan_theta_i) ? slope_x_1 : slope_x_2;
            if (u0 == 0)
                slope_x = 0;
            
            // sample slope_y
            float S;
            if (u1 > 0.5) {
                S = 1.0;
                u1 = 2.0 * (u1 - 0.5);
            }
            else {
                S = -1.0;
                u1 = 2.0 * (0.5 - u1);
            }
            const float z = (u1 * (u1 * (u1 * 0.27385 - 0.73369) + 0.46341)) / (u1 * (u1 * (u1 * 0.093073 + 0.309420) - 1.000000) + 0.597999);
            slope_y = S * z * std::sqrt(1.0 + slope_x * slope_x);
        }
        
        // rotate
        float tmp = std::cos(phi_sv) * slope_x - std::sin(phi_sv) * slope_y;
        slope_y = std::sin(phi_sv) * slope_x + std::cos(phi_sv) * slope_y;
        slope_x = tmp;
        
        // unstretch
        slope_x *= alpha_gx;
        slope_y *= alpha_gy;
        
        *m = normalize(Normal3D(-slope_x, -slope_y, 1));
        float D = evaluate(*m);
        *normalPDF = evaluateSmithG1(v, *m) * absDot(v, *m) * D / std::abs(v.z);
        
        return D;
    }
    
    float GGX::evaluatePDF(const Vector3D &v, const Normal3D &m) const {
        return evaluateSmithG1(v, m) * absDot(v, m) * evaluate(m) / std::abs(v.z);
    }
    
    float GGX::evaluateSmithG1(const Vector3D &v, const Normal3D &m) const {
        float chi = (dot(v, m) / v.z) > 0 ? 1 : 0;
        float theta_v = std::acos(std::clamp(v.z, -1.0f, 1.0f));
        return chi * 2 / (1 + std::sqrt(1 + std::pow(m_alpha_g * std::tan(theta_v), 2)));
    }
}

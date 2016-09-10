//
//  directional_distribution_functions.cpp
//
//  Created by 渡部 心 on 2015/06/28.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "directional_distribution_functions.h"
#include "distributions.h"
#include "SurfaceObject.h"

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
            ret += fs * std::fabs(fsResult.dir_sn.z) * std::fabs(dir_sn.z) / (dirPDF * fsResult.dirPDF);
        }
        return ret.result / (M_PI * (numSamples - numFails));
    }
    
    
    float BSSRDF::sample_distance(float rMax, int16_t wlIdx, float u, float* distPDF) const {
        const float sigma_tr = m_sigma_tr[wlIdx];
        const float z_r = m_z_r[wlIdx];
        const float z_v = m_z_v[wlIdx];
        const float d_r_max = std::sqrt(rMax * rMax + z_r * z_r);
        const float d_v_max = std::sqrt(rMax * rMax + z_v * z_v);
        const float k_d_r_max = z_r / d_r_max * std::exp(-sigma_tr * d_r_max);
        const float k_d_v_max = z_v / d_v_max * std::exp(-sigma_tr * d_v_max);
        const float k = 1.0f / (std::exp(-sigma_tr * z_r) + std::exp(-sigma_tr * z_v) - k_d_r_max - k_d_v_max);
        
        float searchLowBound = 0;
        float searchHighBound = rMax;
        
        float r = 10 * z_r;
        uint32_t count = 0;
        while (true) {
            ++count;
            if (r < searchLowBound || r > searchHighBound)
                r = (searchLowBound + searchHighBound) * 0.5f;
            
            float d_r = std::sqrt(r * r + z_r * z_r);
            float d_v = std::sqrt(r * r + z_v * z_v);
            float term_d_r = z_r / d_r * std::exp(-sigma_tr * d_r);
            float term_d_v = z_v / d_v * std::exp(-sigma_tr * d_v);
            float value = 1 - k * (term_d_r - k_d_r_max + term_d_v - k_d_v_max) - u;
            
            float deriv = k * r * (term_d_r * (1 + sigma_tr * d_r) / (d_r * d_r) + term_d_v * (1 + sigma_tr * d_v) / (d_v * d_v));
            if (std::fabs(value) < 1e-5) {
                *distPDF = deriv;
//                printf("%u, ", count);
                break;
            }
            if (count == 100) {
                *distPDF = 0.0f;
                break;
            }
            
            if (value > 0)
                searchHighBound = r;
            else
                searchLowBound = r;
            
            r -= value / deriv;
        }
        return r;
    }
    
    float BSSRDF::evaludateDistancePDF(float rMax, int16_t wlIdx, float r) const {
        if (r > rMax)
            return 0.0f;
        const float sigma_tr = m_sigma_tr[wlIdx];
        const float z_r = m_z_r[wlIdx];
        const float z_v = m_z_v[wlIdx];
        const float d_r_max = std::sqrt(rMax * rMax + z_r * z_r);
        const float d_v_max = std::sqrt(rMax * rMax + z_v * z_v);
        const float k_d_r_max = z_r / d_r_max * std::exp(-sigma_tr * d_r_max);
        const float k_d_v_max = z_v / d_v_max * std::exp(-sigma_tr * d_v_max);
        const float k = 1.0f / (std::exp(-sigma_tr * z_r) + std::exp(-sigma_tr * z_v) - k_d_r_max - k_d_v_max);
        
        float d_r = std::sqrt(r * r + z_r * z_r);
        float d_v = std::sqrt(r * r + z_v * z_v);
        float term_d_r = z_r / d_r * std::exp(-sigma_tr * d_r);
        float term_d_v = z_v / d_v * std::exp(-sigma_tr * d_v);
        return k * r * (term_d_r * (1 + sigma_tr * d_r) / (d_r * d_r) + term_d_v * (1 + sigma_tr * d_v) / (d_v * d_v));
    }
    
    SampledSpectrum BSSRDF::Rd(float r) const {
        SampledSpectrum d_r = sqrt(m_z_r * m_z_r + r * r);
        SampledSpectrum d_v = sqrt(m_z_v * m_z_v + r * r);
        SampledSpectrum sigma_tr_d_r = m_sigma_tr * d_r;
        SampledSpectrum sigma_tr_d_v = m_sigma_tr * d_v;
        
        SampledSpectrum realTerm = m_z_r * (SampledSpectrum::One + sigma_tr_d_r) * exp(-sigma_tr_d_r) / (d_r * d_r * d_r);
        SampledSpectrum virtualTerm = m_z_v * (SampledSpectrum::One + sigma_tr_d_v) * exp(-sigma_tr_d_v) / (d_v * d_v * d_v);
        
        return 1.0f / (4 * M_PI) * (realTerm + virtualTerm);
    }
    
    SampledSpectrum BSSRDF::sample(const SLR::BSSRDFQuery &query, const SLR::BSSRDFSample &smp, SLR::BSSRDFQueryResult *result) const {
        if (smp.uAuxiliary < 0.25f)  {
            float so = -std::log(smp.uPos[0]) / m_sigma_e[query.wlHint];
            float distPDF = m_sigma_e[query.wlHint] * std::exp(-m_sigma_e[query.wlHint] * so);
            Ray ray{query.surfPt.p, query.dir, query.time, Ray::Epsilon, so};
            Intersection isect;
            if (query.scene.intersect(ray, &isect)) {
                isect.getSurfacePoint(&result->surfPt);
                result->dir = -ray.dir;
                float distProb = std::exp(-m_sigma_e[query.wlHint] * isect.dist);
                result->areaPDF = 1.0f * 0.25f * distProb;
                result->dirPDF = 1.0f;
                return exp(-m_sigma_e * isect.dist);
            }
            auto sampleHenyeyGreenstein = [](float g, float u, float* theta) {
                float cosTheta;
                if (g == 0)
                    cosTheta = 1 - 2 * u;
                else
                    cosTheta = 1 / (2 * g) * (1 + g * g - std::pow((1 - g * g) / (1 - g + 2 * g * u), 2));
                cosTheta = std::clamp(cosTheta, -1.0f, 1.0f);
                *theta = std::acos(cosTheta);
                return (1 - g * g) / (4 * M_PI * std::pow(1 + g * g - 2 * g * cosTheta, 1.5));
            };
            float theta;
            float fp = sampleHenyeyGreenstein(m_g, smp.uDir[0], &theta);
            float phi = 2 * M_PI * smp.uDir[1];
            
            ReferenceFrame scatterFrame;
            scatterFrame.z = ray.dir;
            scatterFrame.z.makeCoordinateSystem(&scatterFrame.x, &scatterFrame.y);
            float sinTheta = std::sin(theta);
            Vector3D scatterDir = scatterFrame.fromLocal(Vector3D(std::cos(phi) * sinTheta, std::sin(phi) * sinTheta, std::cos(theta)));
            
            ray = Ray(ray.org + ray.dir * so, scatterDir, ray.time, 0.0f, INFINITY);
            isect = Intersection();
            if (!query.scene.intersect(ray, &isect)) {
                result->areaPDF = 0.0f;
                result->dirPDF = 0.0f;
                return SampledSpectrum::Zero;
            }
            isect.getSurfacePoint(&result->surfPt);
            result->dir = -ray.dir;
            result->areaPDF = 1.0f * 0.25f * fp;
            result->dirPDF = 1.0f * distPDF;
            
            return m_sigma_s * fp * exp(-m_sigma_e * (so + isect.dist));
        }
        
        float uAxis = (smp.uAuxiliary - 0.25f) / 0.75f;
        
        Vector3D s, t, u;
        float uSelect;
        Vector3D bases[3];
        bases[2] = query.surfPt.gNormal;
        bases[2].makeCoordinateSystem(&bases[0], &bases[1]);
        const float axisProbs[] = {0.5f, 0.25f, 0.25f};
        int8_t selectedAxis = 0;
        // Should it use shading frame for this?
        if (uAxis < 0.5f) {
            selectedAxis = 0;
            uSelect = uAxis / 0.5f;
            s = bases[0];//query.surfPt.shadingFrame.x;
            t = bases[1];//query.surfPt.shadingFrame.y;
            u = bases[2];//query.surfPt.shadingFrame.z;
        }
        else if (smp.uAuxiliary < 0.75f) {
            selectedAxis = 1;
            uSelect = (uAxis - 0.5f) / 0.25f;
            s = bases[1];//query.surfPt.shadingFrame.y;
            t = bases[2];//query.surfPt.shadingFrame.z;
            u = bases[0];//query.surfPt.shadingFrame.x;
        }
        else {
            selectedAxis = 2;
            uSelect = (uAxis - 0.75f) / 0.25f;
            s = bases[2];//query.surfPt.shadingFrame.z;
            t = bases[0];//query.surfPt.shadingFrame.x;
            u = bases[1];//query.surfPt.shadingFrame.y;
        }
        
        float distPDF;
        float rMax = sample_distance(1e+6 * m_z_r[query.wlHint], query.wlHint, 0.999f, &distPDF);
        float dist = sample_distance(rMax, query.wlHint, smp.uPos[0], &distPDF);
        SLRAssert(!std::isnan(dist) && !std::isinf(dist) && !std::isnan(distPDF) && !std::isinf(distPDF),
                  "dist: %g, distPDF: %g, rMax: %g, u: %g, wlIdx: %u", dist, distPDF, rMax, smp.uPos[0], query.wlHint);
        if (distPDF == 0.0f) {
            result->areaPDF = 0.0f;
            result->dirPDF = 0.0f;
            return SampledSpectrum::Zero;
        }
        float theta = 2 * M_PI * smp.uPos[1];
        float areaPDF = distPDF / (2 * M_PI * dist);
        
        float probeLength = 2 * std::sqrt(rMax * rMax - dist * dist);
        // is it safer to spawn a probe ray from beneath the surface?
        Point3D probeOrg = query.surfPt.p + dist * (s * std::cos(theta) + t * std::sin(theta)) - 0.5f * probeLength * u;
        Ray probeRay{probeOrg, u, query.time, 0.0f, probeLength};
        
        std::vector<Intersection> isects;
        while (true) {
            Intersection isect;
            if (!query.scene.intersect(probeRay, &isect))
                break;
            if (isect.getSurfaceMaterial() == query.surfPt.obj->getSurfaceMaterial())
                isects.push_back(isect);
            probeRay.distMin = isect.dist + Ray::Epsilon;
            probeRay.distMax = probeLength;
        }
        if (isects.empty()) {
            result->areaPDF = 0.0f;
            result->dirPDF = 0.0f;
            return SampledSpectrum::Zero;
        }
        
//        uint32_t idxIsect = std::min(uint32_t(uSelect * isects.size()), uint32_t(isects.size()) - 1);
        uint32_t idxIsect = 0;
        float minDiff = INFINITY;
        for (int i = 0; i < isects.size(); ++i) {
            float diff = std::fabs(0.5f * probeLength - isects[i].dist);
            if (diff < minDiff) {
                idxIsect = i;
                minDiff = diff;
            }
        }
        isects[idxIsect].getSurfacePoint(&result->surfPt);
        areaPDF *= axisProbs[selectedAxis] * absDot(u, result->surfPt.gNormal);
        
        Vector3D distVec = result->surfPt.p - query.surfPt.p;
        for (int axis = 0; axis < 3; ++axis) {
            if (axis == selectedAxis)
                continue;
            Vector3D probeAxis = bases[(axis + 2) % 3];
            Vector3D projectedDistVec = distVec - dot(distVec, probeAxis) * probeAxis;
            float otherDist = projectedDistVec.length();
            float otherAxisPDF = evaludateDistancePDF(rMax, query.wlHint, otherDist) / (2 * M_PI * otherDist);
            otherAxisPDF *= axisProbs[axis] * absDot(probeAxis, result->surfPt.gNormal) / 1;// TODO: consider multiple hits.
            SLRAssert(!std::isnan(otherAxisPDF) && !std::isinf(otherAxisPDF),
                      "otherDist: %g, otherAxisPDF: %g, rMax: %g, wlIdx: %u", otherDist, otherAxisPDF, rMax, query.wlHint);
            areaPDF += otherAxisPDF;
        }
        result->areaPDF = 0.75f * areaPDF;
        
        // it might be good to sample this direction by sampling BSDF with randomly sampled the opposite side direction.
        Vector3D dirLocal = cosineSampleHemisphere(smp.uDir[0], smp.uDir[1]);
        result->dirPDF = dirLocal.z / M_PI;
        dirLocal.z *= -1;
        
        ReferenceFrame geometricFrame;
        geometricFrame.z = result->surfPt.gNormal;
        geometricFrame.z.makeCoordinateSystem(&geometricFrame.x, &geometricFrame.y);
        result->dir = geometricFrame.fromLocal(dirLocal);
        
//        float z = dot(distVec, query.surfPt.gNormal);
        SampledSpectrum ret = Rd(distVec.length()) * (absDot(result->surfPt.gNormal, result->dir) / M_PI);
        SLRAssert(!ret.hasInf() && !ret.hasNaN(),
                  "R: %s, u(%g, %g, %g, %g, %g), p: %s, gNormal: %s",
                  ret.toString().c_str(), smp.uPos[0], smp.uPos[1], smp.uDir[0], smp.uDir[1], smp.uAuxiliary,
                  query.surfPt.p.toString().c_str(), query.surfPt.gNormal.toString().c_str());
        
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

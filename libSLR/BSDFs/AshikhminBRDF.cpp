//
//  AshikhminBRDF.cpp
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "AshikhminBRDF.h"
#include "../Core/distributions.h"

namespace SLR {
    SampledSpectrum AshikhminShirleyBRDF::sampleInternal(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult *result) const {
        float iRs = m_Rs.importance(query.wlHint);
        float iRd = m_Rd.importance(query.wlHint);
        float vDotHV = std::fabs(query.dir_sn.z);
        float specularWeight = iRs + (1 - iRs) * std::pow(1.0f - vDotHV, 5);
        float transmissionTerm = 1 - std::pow(1 - vDotHV * 0.5f, 5);
        float diffuseWeght = 28 * iRd / 23 * (1 - iRs) * transmissionTerm * transmissionTerm;
        float sumWeights = specularWeight + diffuseWeght;
        
        float specularDirPDF, diffuseDirPDF;
        float revSpecularDirPDF, revDiffuseDirPDF;
        SampledSpectrum fs;
        if (uComponent * sumWeights < specularWeight) {
            result->dirType = DirectionType::Reflection | DirectionType::HighFreq;
            
            float quad = 2 * M_PI * uDir[1];
            float phi_h = std::atan2(std::sqrt(m_nu + 1) * std::sin(quad), std::sqrt(m_nv + 1) * std::cos(quad));
            float cosphi = std::cos(phi_h);
            float sinphi = std::sin(phi_h);
            float theta_h = std::acos(std::pow(1 - uDir[0], 1.0f / (m_nu * cosphi * cosphi + m_nv * sinphi * sinphi + 1)));
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
            float commonTerm = std::sqrt((m_nu + 1) * (m_nv + 1)) / (8 * M_PI * dotHV) * std::pow(std::fabs(halfv.z), exp);
            SampledSpectrum F = m_Rs + (SampledSpectrum::One - m_Rs) * std::pow(1.0f - dotHV, 5);
            specularDirPDF = commonTerm;
            revSpecularDirPDF = specularDirPDF;
            SampledSpectrum specular_fs = commonTerm / std::fmax(std::fabs(query.dir_sn.z), std::fabs(result->dir_sn.z)) * F;
            
            // calculate PDF to generate the result direction by sampling from diffuse component.
            diffuseDirPDF = std::fabs(result->dir_sn.z) / M_PI;
            revDiffuseDirPDF = std::fabs(query.dir_sn.z) / M_PI;
            SampledSpectrum diffuse_fs = (28 * m_Rd / (23 * M_PI) * (SampledSpectrum::One - m_Rs) *
                                          (1.0f - std::pow(1.0f - std::fabs(query.dir_sn.z) / 2, 5)) *
                                          (1.0f - std::pow(1.0f - std::fabs(result->dir_sn.z) / 2, 5))
                                          );
            
            fs = specular_fs + diffuse_fs;
        }
        else {
            result->dirType = DirectionType::Reflection | DirectionType::LowFreq;
            
            result->dir_sn = cosineSampleHemisphere(uDir[0], uDir[1]);
            diffuseDirPDF = result->dir_sn.z / M_PI;
            revDiffuseDirPDF = std::fabs(query.dir_sn.z) / M_PI;
            result->dir_sn.z *= dot(query.dir_sn, query.gNormal_sn) > 0 ? 1 : -1;
            SampledSpectrum diffuse_fs = (28 * m_Rd / (23 * M_PI) * (SampledSpectrum::One - m_Rs) *
                                          (1.0f - std::pow(1.0f - std::fabs(query.dir_sn.z) / 2, 5)) *
                                          (1.0f - std::pow(1.0f - std::fabs(result->dir_sn.z) / 2, 5))
                                          );
            
            // calculate PDF to generate the result direction by sampling from specular component.
            Vector3D halfv = halfVector(query.dir_sn, result->dir_sn);
            float dotHV = dot(halfv, query.dir_sn);
            float exp = (m_nu * halfv.x * halfv.x + m_nv * halfv.y * halfv.y) / (1 - halfv.z * halfv.z);
            float commonTerm = std::sqrt((m_nu + 1) * (m_nv + 1)) / (8 * M_PI * dotHV) * std::pow(std::fabs(halfv.z), exp);
            SampledSpectrum F = m_Rs + (SampledSpectrum::One - m_Rs) * std::pow(1.0f - dotHV, 5);
            specularDirPDF = commonTerm;
            revSpecularDirPDF = specularDirPDF;
            SampledSpectrum specular_fs = commonTerm / std::fmax(std::fabs(query.dir_sn.z), std::fabs(result->dir_sn.z)) * F;
            
            fs = specular_fs + diffuse_fs;
        }
        
        result->dirPDF = (specularDirPDF * specularWeight + diffuseDirPDF * diffuseWeght) / sumWeights;
        if (result->reverse) {
            float revVDotHV = std::fabs(result->dir_sn.z);
            float revSpecularWeight = iRs + (1 - iRs) * std::pow(1.0f - revVDotHV, 5);
            float revTransmissionTerm = 1 - std::pow(1 - revVDotHV * 0.5f, 5);
            float revDiffuseWeght = 28 * iRd / 23 * (1 - iRs) * revTransmissionTerm * revTransmissionTerm;
            float revSumWeights = revSpecularWeight + revDiffuseWeght;
            
            result->reverse->fs = fs;
            result->reverse->dirPDF = (revSpecularDirPDF * revSpecularWeight + revDiffuseDirPDF * revDiffuseWeght) / revSumWeights;
        }
        return fs;
    }
    
    SampledSpectrum AshikhminShirleyBRDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const {
        if (dir.z * query.dir_sn.z <= 0) {
            if (rev_fs)
                *rev_fs = SampledSpectrum::Zero;
            return SampledSpectrum::Zero;
        }
        Vector3D halfv = halfVector(query.dir_sn, dir);
        float dotHV = dot(halfv, query.dir_sn);
        float exp = (m_nu * halfv.x * halfv.x + m_nv * halfv.y * halfv.y) / (1 - halfv.z * halfv.z);
        float commonTerm = std::sqrt((m_nu + 1) * (m_nv + 1)) / (8 * M_PI * dotHV) * std::pow(std::fabs(halfv.z), exp);
        SampledSpectrum F = m_Rs + (SampledSpectrum::One - m_Rs) * std::pow(1.0f - dotHV, 5);
        SampledSpectrum specular_fs = commonTerm / std::fmax(std::fabs(query.dir_sn.z), std::fabs(dir.z)) * F;
        SampledSpectrum diffuse_fs = (28 * m_Rd / (23 * M_PI) * (SampledSpectrum::One - m_Rs) *
                                      (1.0f - std::pow(1.0f - std::fabs(query.dir_sn.z) / 2, 5)) *
                                      (1.0f - std::pow(1.0f - std::fabs(dir.z) / 2, 5))
                                      );
        SampledSpectrum fs = specular_fs + diffuse_fs;
        if (rev_fs)
            *rev_fs = fs;
        return fs;
    }
    
    float AshikhminShirleyBRDF::evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, float* revPDF) const {
        if (dir.z * query.dir_sn.z <= 0) {
            if (revPDF)
                *revPDF = 0.0f;
            return 0.0f;
        }
        
        Vector3D halfv = halfVector(query.dir_sn, dir);
        float dotHV = dot(halfv, query.dir_sn);
        float exp = (m_nu * halfv.x * halfv.x + m_nv * halfv.y * halfv.y) / (1 - halfv.z * halfv.z);
        float specularDirPDF = std::sqrt((m_nu + 1) * (m_nv + 1)) / (8 * M_PI * dotHV) * std::pow(std::fabs(halfv.z), exp);
        float revSpecularDirPDF = specularDirPDF;
        
        float diffuseDirPDF = std::fabs(dir.z) / M_PI;
        float revDiffuseDirPDF = std::fabs(query.dir_sn.z) / M_PI;
        
        float iRs = m_Rs.importance(query.wlHint);
        float iRd = m_Rd.importance(query.wlHint);
        float vDotHV = std::fabs(query.dir_sn.z);
        float specularWeight = iRs + (1 - iRs) * std::pow(1.0f - vDotHV, 5);
        float transmissionTerm = 1 - std::pow(1 - vDotHV * 0.5f, 5);
        float diffuseWeght = 28 * iRd / 23 * (1 - iRs) * transmissionTerm * transmissionTerm;
        float sumWeights = specularWeight + diffuseWeght;
        
        if (revPDF) {
            float revVDotHV = std::fabs(dir.z);
            float revSpecularWeight = iRs + (1 - iRs) * std::pow(1.0f - revVDotHV, 5);
            float revTransmissionTerm = 1 - std::pow(1 - revVDotHV * 0.5f, 5);
            float revDiffuseWeght = 28 * iRd / 23 * (1 - iRs) * revTransmissionTerm * revTransmissionTerm;
            float revSumWeights = revSpecularWeight + revDiffuseWeght;
            
            *revPDF = (revSpecularDirPDF * revSpecularWeight + revDiffuseDirPDF * revDiffuseWeght) / revSumWeights;
        }
        return (specularDirPDF * specularWeight + diffuseDirPDF * diffuseWeght) / sumWeights;
    }
    
    float AshikhminShirleyBRDF::weightInternal(const SLR::BSDFQuery &query) const {
        // assume the half vector is the same to the normal.
        float iRs = m_Rs.importance(query.wlHint);
        float iRd = m_Rd.importance(query.wlHint);
        float dotHV = std::fabs(query.dir_sn.z);
        float specularWeight = iRs + (1 - iRs) * std::pow(1.0f - dotHV, 5);
        float transmissionTerm = 1 - std::pow(1 - dotHV * 0.5f, 5);
        float diffuseWeght = 28 * iRd / 23 * (1 - iRs) * transmissionTerm * transmissionTerm;
        return specularWeight + diffuseWeght;
    }
    
    SampledSpectrum AshikhminShirleyBRDF::getBaseColorInternal(DirectionType flags) const {
        return m_Rs + (SampledSpectrum::One - m_Rs) * m_Rd;
    }
}

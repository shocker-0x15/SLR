//
//  disney_bsdfs.cpp
//
//  Created by 渡部 心 on 2017/05/04.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "disney_bsdfs.h"

#include "../Core/distributions.h"

namespace SLR {
    DisneyBRDF::DisneyBRDF(const SampledSpectrum &baseColor, float baseColorLuminance,
                           float subsurface, float metallic, float specular, float specularTint, float roughness, float anisotropic, 
                           float sheen, float sheenTint, float clearCoat, float clearCoatGloss) :
    BSDF(DirectionType::Reflection | DirectionType::HighFreq), 
    m_baseColor(baseColor), 
    m_subsurface(subsurface), m_metallic(metallic), m_specular(specular), m_specularTint(specularTint), m_roughness(roughness), m_anisotropic(anisotropic), 
    m_sheen(sheen), m_sheenTint(sheenTint), m_clearCoat(clearCoat), m_clearCoatGloss(clearCoatGloss), 
    m_baseColorLuminance(baseColorLuminance)
    {
        float aspect = std::sqrt(1 - 0.9f * m_anisotropic);
        float alpha_x = std::max(0.001f, m_roughness * m_roughness / aspect);
        float alpha_y = std::max(0.001f, m_roughness * m_roughness * aspect);
        m_base_D = new GGXMicrofacetDistribution(alpha_x, alpha_y);
        
        m_clearcoat_D = new BerryMicrofacetDistribution(0.1f * (1 - clearCoatGloss) + 0.001f * clearCoatGloss);
        m_clearcoat_forG = new GGXMicrofacetDistribution(0.25f, 0.25f);
    }
    
    DisneyBRDF::~DisneyBRDF() {
        delete m_clearcoat_forG;
        delete m_clearcoat_D;
        delete m_base_D;
    }
    
    SampledSpectrum DisneyBRDF::sampleInternal(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult* result) const {
        bool entering = query.dirLocal.z >= 0.0f;
        Vector3D dirL;
        Vector3D dirV = entering ? query.dirLocal : -query.dirLocal;
        
        float expectedFresnelH = std::pow(1 - dirV.z, 5);
        SampledSpectrum tintColor = m_baseColorLuminance > 0 ? m_baseColor / m_baseColorLuminance : SampledSpectrum::One;
        float iBaseColor = m_baseColor.importance(query.wlHint);
        float iSpecularColor = lerp(SampledSpectrum::One, tintColor, m_specularTint).importance(query.wlHint);
        float iSpecularF0 = 0.08f * m_specular * iSpecularColor * (1 - m_metallic) + iBaseColor * m_metallic;
        
        float diffuseWeight = iBaseColor * (1 - m_metallic);
        float specularWeight = iSpecularF0 * (1 - expectedFresnelH) + 1.0f * expectedFresnelH;
        float clearCoatWeight = 0.25f * m_clearCoat * (0.04f * (1 - expectedFresnelH) + 1.0f * expectedFresnelH);
        
        float weights[] = {diffuseWeight, specularWeight, clearCoatWeight};
        float probSelection;
        float sumWeights = 0.0f;
        uint32_t component = sampleDiscrete(weights, 3, uComponent, &probSelection, &sumWeights, &uComponent);
        
        float diffuseDirPDF, specularDirPDF, clearCoatDirPDF;
        float revDiffuseDirPDF, revSpecularDirPDF, revClearCoatDirPDF;
        SampledSpectrum fs, rev_fs;
        if (component == 0) {
            result->sampledType = DirectionType::Reflection | DirectionType::LowFreq;
            
            // sample based on cosine distribution.
            dirL = cosineSampleHemisphere(uDir[0], uDir[1]);
            result->dirLocal = entering ? dirL : -dirL;
            diffuseDirPDF = dirL.z / M_PI;
            revDiffuseDirPDF = dirV.z / M_PI;
            
            Normal3D m = halfVector(dirL, dirV);
            float dotHV = dot(dirV, m);
            float commonPDFTerm = 1.0f / (4 * dotHV);
            
            specularDirPDF = commonPDFTerm * m_base_D->evaluatePDF(dirV, m);
            revSpecularDirPDF = commonPDFTerm * m_base_D->evaluatePDF(dirL, m);
            
            clearCoatDirPDF = commonPDFTerm * m_clearcoat_D->evaluatePDF(dirV, m);
            revClearCoatDirPDF = commonPDFTerm * m_clearcoat_D->evaluatePDF(dirL, m);
            
            fs = evaluateInternal(query, result->dirLocal, &rev_fs);
        }
        else if (component == 1) { 
            result->sampledType = DirectionType::Reflection | DirectionType::HighFreq;
            
            // sample based on the base specular microfacet distribution.
            Normal3D m;
            float mPDF;
            m_base_D->sample(dirV, uDir[0], uDir[1], &m, &mPDF);
            float dotHV = dot(dirV, m);
            dirL = 2 * dotHV * m - query.dirLocal;
            result->dirLocal = entering ? dirL : -dirL;
            if (dirL.z * dirV.z <= 0) {
                result->dirPDF = 0.0f;
                return SampledSpectrum::Zero;
            }
            float commonPDFTerm = 1.0f / (4 * dotHV);
            specularDirPDF = commonPDFTerm * mPDF;
            revSpecularDirPDF = commonPDFTerm * m_base_D->evaluatePDF(dirL, m);
            
            diffuseDirPDF = dirL.z / M_PI;
            revDiffuseDirPDF = dirV.z / M_PI;
            
            clearCoatDirPDF = commonPDFTerm * m_clearcoat_D->evaluatePDF(dirV, m);
            revClearCoatDirPDF = commonPDFTerm * m_clearcoat_D->evaluatePDF(dirL, m);
            
            fs = evaluateInternal(query, result->dirLocal, &rev_fs);
        }
        else {
            result->sampledType = DirectionType::Reflection | DirectionType::HighFreq;
            
            // sample based on the clear coat microfacet distribution.
            Normal3D m;
            float mPDF;
            
            m_clearcoat_D->sample(dirV, uDir[0], uDir[1], &m, &mPDF);
            float dotHV = dot(dirV, m);
            float commonPDFTerm = 1.0f / (4 * dotHV);
            
            dirL = 2 * dotHV * m - query.dirLocal;
            result->dirLocal = entering ? dirL : -dirL;
            if (dirL.z * dirV.z <= 0) {
                result->dirPDF = 0.0f;
                return SampledSpectrum::Zero;
            }
            clearCoatDirPDF = commonPDFTerm * mPDF;
            revClearCoatDirPDF = commonPDFTerm * m_clearcoat_D->evaluatePDF(dirL, m);
            
            diffuseDirPDF = dirL.z / M_PI;
            revDiffuseDirPDF = dirV.z / M_PI;
            
            specularDirPDF = commonPDFTerm * m_base_D->evaluatePDF(dirV, m);
            revSpecularDirPDF = commonPDFTerm * m_base_D->evaluatePDF(dirL, m);
            
            fs = evaluateInternal(query, result->dirLocal, &rev_fs);
        }
        
        result->dirPDF = (diffuseDirPDF * diffuseWeight + 
                          specularDirPDF * specularWeight + 
                          clearCoatDirPDF * clearCoatWeight) / sumWeights;
        if (query.requestReverse) {
            float revExpectedFresnelH = std::pow(1 - dirL.z, 5);
            
            float revDiffuseWeight = diffuseWeight;
            float revSpecularWeight = iSpecularF0 * (1 - revExpectedFresnelH) + 1.0f * revExpectedFresnelH;
            float revClearCoatWeight = 0.25f * m_clearCoat * (0.04f * (1 - revExpectedFresnelH) + 1.0f * revExpectedFresnelH);
            float revSumWeights = revDiffuseWeight + revSpecularWeight + revClearCoatWeight;
            
            result->reverse.value = fs;
            result->reverse.dirPDF = (revDiffuseDirPDF * revDiffuseWeight + 
                                      revSpecularDirPDF * revSpecularWeight + 
                                      revClearCoatDirPDF * revClearCoatWeight) / revSumWeights;
        }
        return fs;
    }
    
    SampledSpectrum DisneyBRDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const {
        if (dir.z * query.dirLocal.z <= 0) {
            if (query.requestReverse)
                *rev_fs = SampledSpectrum::Zero;
            return SampledSpectrum::Zero;
        }
        
        bool entering = query.dirLocal.z >= 0.0f;
        Vector3D dirL = entering ? dir : -dir;
        Vector3D dirV = entering ? query.dirLocal : -query.dirLocal;
        
        Normal3D m = halfVector(dirL, dirV);
        float dotLH = dot(dirL, m);
        
        // Fresnel factors
        float fresnelL = std::pow(1 - dirL.z, 5);
        float fresnelV = std::pow(1 - dirV.z, 5);
        float fresnelH = std::pow(1 - dotLH, 5);
        
        // tintColor is calculated by normalizing luminance component of the baseColor.
        SampledSpectrum tintColor = m_baseColorLuminance > 0 ? m_baseColor / m_baseColorLuminance : SampledSpectrum::One;
        SampledSpectrum specularColor = lerp(SampledSpectrum::One, tintColor, m_specularTint);
        SampledSpectrum sheenColor = lerp(SampledSpectrum::One, tintColor, m_sheenTint);
        SampledSpectrum specularF0Color = lerp(0.08f * m_specular * specularColor, m_baseColor, m_metallic);
        
        // Base - Diffuse Term
        // go from 1 at normal incidence to 0.5 at grazing and mix in diffuse retro-reflection based on roughness.
        float F_D90 = 0.5f + 2 * m_roughness * dotLH * dotLH;
        float diffuseFresnel = (1 + (F_D90 - 1) * fresnelL) * (1 + (F_D90 - 1) * fresnelV); // new name suggustion: diffuseFresnelThroughput?
        //                     ((1 - fresnelL) * 1 + fresnelL * F_D90) ... 1 at normal incidence and F_D90 at grazing.
        
        // Base - Subsurface Term
        // based on Hanrahan-Krueger BRDF approximation of isotropic BSSRDF.
        // F_SS90 used to "flatten" retroreflection based on roughness.
        // 1.25 scale is used to (roughly) preserve albedo.
        float F_SS90 = m_roughness * dotLH * dotLH;
        float ssFresnel = (1 + (F_SS90 - 1) * fresnelL) * (1 + (F_SS90 - 1) * fresnelV);
        float ssCoeff = 1.25f * (ssFresnel * (1.0f / (dir.z + query.dirLocal.z) - 0.5f) + 0.5f);
        
        // Base - Specular Term
        float base_D = m_base_D->evaluate(m);
        float base_G = m_base_D->evaluateSmithG1(dirL, m) * m_base_D->evaluateSmithG1(dirV, m);
        SampledSpectrum base_F = lerp(specularF0Color, SampledSpectrum::One, fresnelH); // grazing specular is always achromatic.
        
        // Base - Sheen Term
        // stronger at grazing angle, primarily intended for cloth.
        SampledSpectrum sheenValue = fresnelH * m_sheen * sheenColor;
        
        // Clear coat (IoR = 1.5 -> F0 = 0.04)
        float coat_D = m_clearcoat_D->evaluate(m);
        float coat_G = m_clearcoat_forG->evaluateSmithG1(dirL, m) * m_clearcoat_forG->evaluateSmithG1(dirV, m);
        float coat_F = 0.04f * (1 - fresnelH) + 1.0f * fresnelH;
        
        // Disney's implementation in BRDF explorer implicitly includes this denominator's value in shadowing-masking factors.
        float microfacetDenom = 4 * dir.z * query.dirLocal.z;
        float clearCoatValue = 0.25f * m_clearCoat * coat_F * coat_D * coat_G / microfacetDenom;
        SampledSpectrum baseSpecularValue = base_F * ((base_D * base_G) / microfacetDenom);
        SampledSpectrum baseDiffuseValue = (((1 - m_subsurface) * diffuseFresnel + m_subsurface * ssCoeff) / M_PI * m_baseColor + sheenValue) * (1 - m_metallic);
        
        SampledSpectrum ret = baseDiffuseValue + baseSpecularValue + clearCoatValue;
        if (query.requestReverse)
            *rev_fs = ret;
        
        return ret;
    }
    
    float DisneyBRDF::evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, float* revPDF) const {
        if (dir.z * query.dirLocal.z <= 0) {
            if (query.requestReverse)
                *revPDF = 0.0f;
            return 0.0f;
        }
        
        bool entering = query.dirLocal.z >= 0.0f;
        Vector3D dirL = entering ? dir : -dir;
        Vector3D dirV = entering ? query.dirLocal : -query.dirLocal;
        
        Normal3D m = halfVector(dirL, dirV);
        float dotHV = dot(dirV, m);
        float commonPDFTerm = 1.0f / (4 * dotHV);
                
        float expectedFresnelH = std::pow(1 - dirV.z, 5);
        SampledSpectrum tintColor = m_baseColorLuminance > 0 ? m_baseColor / m_baseColorLuminance : SampledSpectrum::One;
        float iBaseColor = m_baseColor.importance(query.wlHint);
        float iSpecularColor = lerp(SampledSpectrum::One, tintColor, m_specularTint).importance(query.wlHint);
        float iSpecularF0 = 0.08f * m_specular * iSpecularColor * (1 - m_metallic) + iBaseColor * m_metallic;
        
        float diffuseWeight = iBaseColor * (1 - m_metallic);
        float specularWeight = iSpecularF0 * (1 - expectedFresnelH) + 1.0f * expectedFresnelH;
        float clearCoatWeight = 0.25f * m_clearCoat * (0.04f * (1 - expectedFresnelH) + 1.0f * expectedFresnelH);
        
        float sumWeights = diffuseWeight + specularWeight + clearCoatWeight;
        
        float diffuseDirPDF = dirL.z / M_PI;
        float specularDirPDF = commonPDFTerm * m_base_D->evaluatePDF(dirV, m);
        float clearCoatDirPDF = commonPDFTerm * m_clearcoat_D->evaluatePDF(dirV, m);
        
        if (query.requestReverse) {
            float revExpectedFresnelH = std::pow(1 - dirL.z, 5);
            
            float revDiffuseWeight = diffuseWeight;
            float revSpecularWeight = iSpecularF0 * (1 - revExpectedFresnelH) + 1.0f * revExpectedFresnelH;
            float revClearCoatWeight = 0.25f * m_clearCoat * (0.04f * (1 - revExpectedFresnelH) + 1.0f * revExpectedFresnelH);
            float revSumWeights = revDiffuseWeight + revSpecularWeight + revClearCoatWeight;
            
            float revDiffuseDirPDF = dirV.z / M_PI;
            float revSpecularDirPDF = commonPDFTerm * m_base_D->evaluatePDF(dirL, m);
            float revClearCoatDirPDF = commonPDFTerm * m_clearcoat_D->evaluatePDF(dirL, m);
            
            *revPDF = (revDiffuseDirPDF * revDiffuseWeight + 
                       revSpecularDirPDF * revSpecularWeight + 
                       revClearCoatDirPDF * revClearCoatWeight) / revSumWeights;
        }
        
        return (diffuseDirPDF * diffuseWeight + 
                specularDirPDF * specularWeight + 
                clearCoatDirPDF * clearCoatWeight) / sumWeights;
    }
    
    float DisneyBRDF::weightInternal(const SLR::BSDFQuery &query) const {        
        float iBaseColor = m_baseColor.importance(query.wlHint);
        float diffuseWeight = iBaseColor * (1 - m_metallic);
        float specularWeight = 0.08f * m_specular * (1 - m_metallic) + iBaseColor * m_metallic;
        float clearCoatWeight = 0.25f * m_clearCoat;
        
        return diffuseWeight + specularWeight + clearCoatWeight;
    }
    
    SampledSpectrum DisneyBRDF::getBaseColorInternal(SLR::DirectionType flags) const {
        return m_baseColor;
    }
}

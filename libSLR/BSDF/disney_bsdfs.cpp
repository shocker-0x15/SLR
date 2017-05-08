//
//  disney_bsdfs.cpp
//
//  Created by 渡部 心 on 2017/05/04.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "disney_bsdfs.h"

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
        SLRAssert_NotImplemented();
        return SampledSpectrum::Zero;
    }
    
    SampledSpectrum DisneyBRDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const {
        Vector3D halfVector = normalize(query.dirLocal + dir);
        float dotLH = dot(dir, halfVector);
        
        // Fresnel factors
        float fresnelL = std::pow(1 - dir.z, 5);
        float fresnelV = std::pow(1 - query.dirLocal.z, 5);
        float fresnelH = std::pow(1 - dotLH, 5);
        
        // tintColor is calculated by normalizing luminance component of the baseColor.
        SampledSpectrum tintColor = m_baseColorLuminance > 0 ? m_baseColor / m_baseColorLuminance : SampledSpectrum::One;
        SampledSpectrum specularColor = lerp(0.08f * m_specular * lerp(SampledSpectrum::One, tintColor, m_specularTint), m_baseColor, m_metallic);
        SampledSpectrum sheenColor = lerp(SampledSpectrum::One, tintColor, m_sheenTint);
        
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
        float base_D = m_base_D->evaluate(halfVector);
        float base_G = m_base_D->evaluateSmithG1(dir, halfVector) * m_base_D->evaluateSmithG1(query.dirLocal, halfVector);
        SampledSpectrum base_F = lerp(specularColor, SampledSpectrum::One, fresnelH); // grazing specular is always achromatic.
        
        // Base - Sheen Term
        // stronger at grazing angle, primarily intended for cloth.
        SampledSpectrum sheenValue = fresnelH * m_sheen * sheenColor;
        
        // Clear coat (IoR = 1.5 -> F0 = 0.04)
        float coat_D = m_clearcoat_D->evaluate(halfVector);
        float coat_G = m_clearcoat_forG->evaluateSmithG1(dir, halfVector) * m_clearcoat_forG->evaluateSmithG1(query.dirLocal, halfVector);
        float coat_F = 0.04f * (1 - fresnelH) + 1.0f * fresnelH;
        
        // Disney's implementation in BRDF explorer implicitly includes this denominator's value in shadowing-masking factors.
        float microfacetDenom = 4 * dir.z * query.dirLocal.z;
        float clearCoatValue = 0.25f * m_clearCoat * coat_F * coat_D * coat_G / microfacetDenom;
        SampledSpectrum baseSpecularValue = base_F * ((base_D * base_G) / microfacetDenom);
        SampledSpectrum baseDiffuseValue = (((1 - m_subsurface) * diffuseFresnel + m_subsurface * ssCoeff) / M_PI * m_baseColor + sheenValue) * (1 - m_metallic);
        
        return baseDiffuseValue + baseSpecularValue + clearCoatValue;
    }
    
    float DisneyBRDF::evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, float* revPDF) const {
        SLRAssert_NotImplemented();
        return 0.0f;
    }
    
    float DisneyBRDF::weightInternal(const SLR::BSDFQuery &query) const {
        SLRAssert_NotImplemented();
        return 0.0f;
    }
    
    SampledSpectrum DisneyBRDF::getBaseColorInternal(SLR::DirectionType flags) const {
        return m_baseColor;
    }
}

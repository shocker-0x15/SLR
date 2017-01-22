//
//  surface_materials.cpp
//
//  Created by 渡部 心 on 2015/10/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "surface_materials.hpp"
#include <libSLR/SurfaceMaterials/surface_material_headers.h>
#include "textures.hpp"
#include "Scene.h"

namespace SLRSceneGraph {
    SurfaceMaterial::~SurfaceMaterial() {
        delete m_rawData;
    }
    
    
    
    EmitterSurfaceProperty::~EmitterSurfaceProperty() {
        delete m_rawData;
    }
    
    EmitterSurfaceMaterial::EmitterSurfaceMaterial(const SurfaceMaterialRef &mat, const EmitterSurfacePropertyRef &emit) :
    m_mat(mat), m_emit(emit) {
        m_rawData = new SLR::EmitterSurfaceMaterial(mat->getRaw(), emit->getRaw());
    }
    
    
    
    SVFresnel::~SVFresnel() {
        delete m_rawData;
    }
    
    SVFresnelNoOp::SVFresnelNoOp() {
        m_rawData = new SLR::SVFresnelNoOp();
    }
    
    SVFresnelConductor::SVFresnelConductor(const SpectrumTextureRef &eta, const SpectrumTextureRef &k) : m_eta(eta), m_k(k) {
        m_rawData = new SLR::SVFresnelConductor(eta->getRaw(), k->getRaw());
    }
    
    SVFresnelDielectric::SVFresnelDielectric(const SpectrumTextureRef &etaExt, const SpectrumTextureRef &etaInt) : m_etaExt(etaExt), m_etaInt(etaInt) {
        m_rawData = new SLR::SVFresnelDielectric(etaExt->getRaw(), etaInt->getRaw());
    }
    
    
    SVMicrofacetDistribution::~SVMicrofacetDistribution() {
        delete m_rawData;
    }
    
    SVGGX::SVGGX(const FloatTextureRef &alpha_g) : m_alpha_g(alpha_g) {
        m_rawData = new SLR::SVGGX(alpha_g->getRaw());
    }
    
    
    
    DiffuseReflection::DiffuseReflection(const SpectrumTextureRef &reflectance, const FloatTextureRef &sigma) :
    m_reflectance(reflectance), m_sigma(sigma) {
        m_rawData = new SLR::DiffuseReflection(reflectance->getRaw(), sigma ? sigma->getRaw() : nullptr);
    }
    
    SpecularReflection::SpecularReflection(const SpectrumTextureRef &coeffR, const SpectrumTextureRef &eta, const SpectrumTextureRef &k) :
    m_coeffR(coeffR), m_eta(eta), m_k(k) {
        m_rawData = new SLR::SpecularReflection(coeffR->getRaw(), eta->getRaw(), k->getRaw());
    }
    
    SpecularScattering::SpecularScattering(const SpectrumTextureRef &coeff, const SpectrumTextureRef &etaExt, const SpectrumTextureRef &etaInt) :
    m_coeff(coeff), m_etaExt(etaExt), m_etaInt(etaInt) {
        m_rawData = new SLR::SpecularScattering(coeff->getRaw(), etaExt->getRaw(), etaInt->getRaw());
    }
    
    InverseSurfaceMaterial::InverseSurfaceMaterial(const SurfaceMaterialRef &baseMat) :
    m_baseMat(baseMat) {
        m_rawData = new SLR::InverseSurfaceMaterial(baseMat->getRaw());
    }
    
    AshikhminShirleyReflection::AshikhminShirleyReflection(const SpectrumTextureRef &Rs, const SpectrumTextureRef &Rd, const FloatTextureRef &nu, const FloatTextureRef &nv) :
    m_Rs(Rs), m_Rd(Rd), m_nu(nu), m_nv(nv) {
        m_rawData = new SLR::AshikhminShirleyReflection(Rs->getRaw(), Rd->getRaw(), nu->getRaw(), nv->getRaw());
    }
    
    ModifiedWardDurReflection::ModifiedWardDurReflection(const SpectrumTextureRef &reflectance, const FloatTextureRef &anisoX, const FloatTextureRef &anisoY) :
    m_reflectance(reflectance), m_anisoX(anisoX), m_anisoY(anisoY) {
        m_rawData = new SLR::ModifiedWardDurReflection(reflectance->getRaw(), anisoX->getRaw(), anisoY->getRaw());
    }
    
    MicrofacetReflection::MicrofacetReflection(const SpectrumTextureRef &eta, const SpectrumTextureRef &k, const SVMicrofacetDistributionRef &D) :
    m_eta(eta), m_k(k), m_D(D) {
        m_rawData = new SLR::MicrofacetReflection(eta->getRaw(), k->getRaw(), D->getRaw());
    }
    
    MicrofacetScattering::MicrofacetScattering(const SpectrumTextureRef &etaExt, const SpectrumTextureRef &etaInt, const SVMicrofacetDistributionRef &D) :
    m_etaExt(etaExt), m_etaInt(etaInt), m_D(D) {
        m_rawData = new SLR::MicrofacetScattering(etaExt->getRaw(), etaInt->getRaw(), D->getRaw());
    }
    
    SummedSurfaceMaterial::SummedSurfaceMaterial(const SurfaceMaterialRef &m0, const SurfaceMaterialRef &m1) :
    m_mat0(m0), m_mat1(m1) {
        m_rawData = new SLR::SummedSurfaceMaterial(m0->getRaw(), m1->getRaw());
    }
    
    MixedSurfaceMaterial::MixedSurfaceMaterial(const SurfaceMaterialRef &m0, const SurfaceMaterialRef &m1, const FloatTextureRef &factor) :
    m_mat0(m0), m_mat1(m1), m_factor(factor) {
        m_rawData = new SLR::MixedSurfaceMaterial(m0->getRaw(), m1->getRaw(), factor->getRaw());
    }
    
    
    
    DiffuseEmission::DiffuseEmission(const SpectrumTextureRef &emittance) :
    m_emittance(emittance) {
        m_rawData = new SLR::DiffuseEmission(emittance->getRaw());
    };
    
    IBLEmission::IBLEmission(const SceneWRef &scene, const SpectrumTextureRef &coeffM, float scale) :
    m_scene(scene), m_coeffM(coeffM), m_scale(scale) {
        SceneRef sceneRef = m_scene.lock();
        m_rawData = new SLR::IBLEmission(sceneRef->getRaw(), coeffM->getRaw(), scale);
    };
    
    
    
    SurfaceMaterialRef SurfaceMaterial::createMatte(const SpectrumTextureRef &reflectance, const FloatTextureRef &sigma) {
        return createShared<DiffuseReflection>(reflectance, sigma);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createMetal(const SpectrumTextureRef &coeffR, const SpectrumTextureRef &eta, const SpectrumTextureRef &k) {
        return createShared<SpecularReflection>(coeffR, eta, k);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createGlass(const SpectrumTextureRef &coeff, const SpectrumTextureRef &etaExt, const SpectrumTextureRef &etaInt) {
        return createShared<SpecularScattering>(coeff, etaExt, etaInt);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createModifiedWardDur(const SpectrumTextureRef &reflectance, const FloatTextureRef &anisoX, const FloatTextureRef &anisoY) {
        return createShared<ModifiedWardDurReflection>(reflectance, anisoX, anisoY);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createAshikhminShirley(const SpectrumTextureRef &Rd, const SpectrumTextureRef &Rs, const FloatTextureRef &nu, const FloatTextureRef &nv) {
        return createShared<AshikhminShirleyReflection>(Rs, Rd, nu, nv);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createMicrofacetMetal(const SpectrumTextureRef &eta, const SpectrumTextureRef &k, const FloatTextureRef &alpha_g) {
        SVMicrofacetDistributionRef dist = createShared<SVGGX>(alpha_g);
        return createShared<MicrofacetReflection>(eta, k, dist);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createMicrofacetGlass(const SpectrumTextureRef &etaExt, const SpectrumTextureRef &etaInt, const FloatTextureRef &alpha_g) {
        SVMicrofacetDistributionRef dist = createShared<SVGGX>(alpha_g);
        return createShared<MicrofacetScattering>(etaExt, etaInt, dist);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createInverseMaterial(const SurfaceMaterialRef &baseMat) {
        return createShared<InverseSurfaceMaterial>(baseMat);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createSummedMaterial(const SurfaceMaterialRef &mat0, const SurfaceMaterialRef &mat1) {
        return createShared<SummedSurfaceMaterial>(mat0, mat1);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createMixedMaterial(const SurfaceMaterialRef &mat0, const SurfaceMaterialRef &mat1, const FloatTextureRef &factor) {
        return createShared<MixedSurfaceMaterial>(mat0, mat1, factor);
    }
    
    EmitterSurfacePropertyRef SurfaceMaterial::createDiffuseEmitter(const SpectrumTextureRef &emittance) {
        return createShared<DiffuseEmission>(emittance);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createEmitterSurfaceMaterial(const SurfaceMaterialRef &mat, const EmitterSurfacePropertyRef &emit) {
        return createShared<EmitterSurfaceMaterial>(mat, emit);
    }
}

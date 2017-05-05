//
//  surface_materials.cpp
//
//  Created by 渡部 心 on 2015/10/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "surface_materials.h"

#include <libSLR/SurfaceMaterial/surface_material_headers.h>
#include "textures.h"
#include "Scene/Scene.h"

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
    
    GGXSVMicrofacetDistribution::GGXSVMicrofacetDistribution(const FloatTextureRef &alpha_gx, const FloatTextureRef &alpha_gy) : m_alpha_gx(alpha_gx), m_alpha_gy(alpha_gy) {
        m_rawData = new SLR::GGXSVMicrofacetDistribution(alpha_gx->getRaw(), alpha_gy->getRaw());
    }
    
    
    
    DiffuseReflectionSurfaceMaterial::DiffuseReflectionSurfaceMaterial(const SpectrumTextureRef &reflectance, const FloatTextureRef &sigma) :
    m_reflectance(reflectance), m_sigma(sigma) {
        m_rawData = new SLR::DiffuseReflectionSurfaceMaterial(reflectance->getRaw(), sigma ? sigma->getRaw() : nullptr);
    }
    
    SpecularReflectionSurfaceMaterial::SpecularReflectionSurfaceMaterial(const SpectrumTextureRef &coeffR, const SpectrumTextureRef &eta, const SpectrumTextureRef &k) :
    m_coeffR(coeffR), m_eta(eta), m_k(k) {
        m_rawData = new SLR::SpecularReflectionSurfaceMaterial(coeffR->getRaw(), eta->getRaw(), k->getRaw());
    }
    
    SpecularScatteringSurfaceMaterial::SpecularScatteringSurfaceMaterial(const SpectrumTextureRef &coeff, const SpectrumTextureRef &etaExt, const SpectrumTextureRef &etaInt) :
    m_coeff(coeff), m_etaExt(etaExt), m_etaInt(etaInt) {
        m_rawData = new SLR::SpecularScatteringSurfaceMaterial(coeff->getRaw(), etaExt->getRaw(), etaInt->getRaw());
    }
    
    FlippedSurfaceMaterial::FlippedSurfaceMaterial(const SurfaceMaterialRef &baseMat) :
    m_baseMat(baseMat) {
        m_rawData = new SLR::FlippedSurfaceMaterial(baseMat->getRaw());
    }
    
    AshikhminShirleyReflectionSurfaceMaterial::AshikhminShirleyReflectionSurfaceMaterial(const SpectrumTextureRef &Rs, const SpectrumTextureRef &Rd, const FloatTextureRef &nu, const FloatTextureRef &nv) :
    m_Rs(Rs), m_Rd(Rd), m_nu(nu), m_nv(nv) {
        m_rawData = new SLR::AshikhminShirleyReflectionSurfaceMaterial(Rs->getRaw(), Rd->getRaw(), nu->getRaw(), nv->getRaw());
    }
    
    ModifiedWardDurReflectionSurfaceMaterial::ModifiedWardDurReflectionSurfaceMaterial(const SpectrumTextureRef &reflectance, const FloatTextureRef &anisoX, const FloatTextureRef &anisoY) :
    m_reflectance(reflectance), m_anisoX(anisoX), m_anisoY(anisoY) {
        m_rawData = new SLR::ModifiedWardDurReflectionSurfaceMaterial(reflectance->getRaw(), anisoX->getRaw(), anisoY->getRaw());
    }
    
    MicrofacetReflectionSurfaceMaterial::MicrofacetReflectionSurfaceMaterial(const SpectrumTextureRef &eta, const SpectrumTextureRef &k, const SVMicrofacetDistributionRef &D) :
    m_eta(eta), m_k(k), m_D(D) {
        m_rawData = new SLR::MicrofacetReflectionSurfaceMaterial(eta->getRaw(), k->getRaw(), D->getRaw());
    }
    
    MicrofacetScatteringSurfaceMaterial::MicrofacetScatteringSurfaceMaterial(const SpectrumTextureRef &etaExt, const SpectrumTextureRef &etaInt, const SVMicrofacetDistributionRef &D) :
    m_etaExt(etaExt), m_etaInt(etaInt), m_D(D) {
        m_rawData = new SLR::MicrofacetScatteringSurfaceMaterial(etaExt->getRaw(), etaInt->getRaw(), D->getRaw());
    }
    
    SummedSurfaceMaterial::SummedSurfaceMaterial(const SurfaceMaterialRef &m0, const SurfaceMaterialRef &m1) :
    m_mat0(m0), m_mat1(m1) {
        m_rawData = new SLR::SummedSurfaceMaterial(m0->getRaw(), m1->getRaw());
    }
    
    MixedSurfaceMaterial::MixedSurfaceMaterial(const SurfaceMaterialRef &m0, const SurfaceMaterialRef &m1, const FloatTextureRef &factor) :
    m_mat0(m0), m_mat1(m1), m_factor(factor) {
        m_rawData = new SLR::MixedSurfaceMaterial(m0->getRaw(), m1->getRaw(), factor->getRaw());
    }
    
    
    
    DiffuseEmitterSurfaceProperty::DiffuseEmitterSurfaceProperty(const SpectrumTextureRef &emittance) :
    m_emittance(emittance) {
        m_rawData = new SLR::DiffuseEmitterSurfaceProperty(emittance->getRaw());
    };
    
    IdealDirectionalEmitterSurfaceProperty::IdealDirectionalEmitterSurfaceProperty(const SpectrumTextureRef &emittance, const SLR::Vector3D &dir) :
    m_emittance(emittance), m_direction(dir) {
        m_rawData = new SLR::IdealDirectionalEmitterSurfaceProperty(emittance->getRaw(), dir);
    };
    
    IBLEmitterSurfaceProperty::IBLEmitterSurfaceProperty(const SceneWRef &scene, const SpectrumTextureRef &coeffM, float scale) :
    m_scene(scene), m_coeffM(coeffM), m_scale(scale) {
        SceneRef sceneRef = m_scene.lock();
        m_rawData = new SLR::IBLEmitterSurfaceProperty(sceneRef->getRaw(), coeffM->getRaw(), scale);
    };
    
    
    
    SurfaceMaterialRef SurfaceMaterial::createMatte(const SpectrumTextureRef &reflectance, const FloatTextureRef &sigma) {
        return createShared<DiffuseReflectionSurfaceMaterial>(reflectance, sigma);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createMetal(const SpectrumTextureRef &coeffR, const SpectrumTextureRef &eta, const SpectrumTextureRef &k) {
        return createShared<SpecularReflectionSurfaceMaterial>(coeffR, eta, k);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createGlass(const SpectrumTextureRef &coeff, const SpectrumTextureRef &etaExt, const SpectrumTextureRef &etaInt) {
        return createShared<SpecularScatteringSurfaceMaterial>(coeff, etaExt, etaInt);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createModifiedWardDur(const SpectrumTextureRef &reflectance, const FloatTextureRef &anisoX, const FloatTextureRef &anisoY) {
        return createShared<ModifiedWardDurReflectionSurfaceMaterial>(reflectance, anisoX, anisoY);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createAshikhminShirley(const SpectrumTextureRef &Rd, const SpectrumTextureRef &Rs, const FloatTextureRef &nu, const FloatTextureRef &nv) {
        return createShared<AshikhminShirleyReflectionSurfaceMaterial>(Rs, Rd, nu, nv);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createMicrofacetMetal(const SpectrumTextureRef &eta, const SpectrumTextureRef &k, const FloatTextureRef &alpha_gx, const FloatTextureRef &alpha_gy) {
        SVMicrofacetDistributionRef dist = createShared<GGXSVMicrofacetDistribution>(alpha_gx, alpha_gy);
        return createShared<MicrofacetReflectionSurfaceMaterial>(eta, k, dist);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createMicrofacetGlass(const SpectrumTextureRef &etaExt, const SpectrumTextureRef &etaInt, const FloatTextureRef &alpha_gx, const FloatTextureRef &alpha_gy) {
        SVMicrofacetDistributionRef dist = createShared<GGXSVMicrofacetDistribution>(alpha_gx, alpha_gy);
        return createShared<MicrofacetScatteringSurfaceMaterial>(etaExt, etaInt, dist);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createFlippedMaterial(const SurfaceMaterialRef &baseMat) {
        return createShared<FlippedSurfaceMaterial>(baseMat);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createSummedMaterial(const SurfaceMaterialRef &mat0, const SurfaceMaterialRef &mat1) {
        return createShared<SummedSurfaceMaterial>(mat0, mat1);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createMixedMaterial(const SurfaceMaterialRef &mat0, const SurfaceMaterialRef &mat1, const FloatTextureRef &factor) {
        return createShared<MixedSurfaceMaterial>(mat0, mat1, factor);
    }
    
    EmitterSurfacePropertyRef SurfaceMaterial::createDiffuseEmitter(const SpectrumTextureRef &emittance) {
        return createShared<DiffuseEmitterSurfaceProperty>(emittance);
    }
    
    EmitterSurfacePropertyRef SurfaceMaterial::createIdealDirectionalEmitter(const SpectrumTextureRef &emittance, const SLR::Vector3D &dir) {
        return createShared<IdealDirectionalEmitterSurfaceProperty>(emittance, dir);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createEmitterSurfaceMaterial(const SurfaceMaterialRef &mat, const EmitterSurfacePropertyRef &emit) {
        return createShared<EmitterSurfaceMaterial>(mat, emit);
    }
}

//
//  surface_materials.cpp
//
//  Created by 渡部 心 on 2015/10/06.
//  Copyright © 2015年 渡部 心. All rights reserved.
//

#include "surface_materials.hpp"
#include <libSLR/SurfaceMaterials/surface_material_headers.h>
#include "textures.hpp"

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
    
    SpatialFresnel::~SpatialFresnel() {
        delete m_rawData;
    }
    
    SpatialFresnelNoOp::SpatialFresnelNoOp() {
        m_rawData = new SLR::SpatialFresnelNoOp();
    }
    
    SpatialFresnelConductor::SpatialFresnelConductor(const SpectrumTextureRef &eta, const SpectrumTextureRef &k) : m_eta(eta), m_k(k) {
        m_rawData = new SLR::SpatialFresnelConductor(eta->getRaw(), k->getRaw());
    }
    
    SpatialFresnelDielectric::SpatialFresnelDielectric(const SpectrumTextureRef &etaExt, const SpectrumTextureRef &etaInt) : m_etaExt(etaExt), m_etaInt(etaInt) {
        m_rawData = new SLR::SpatialFresnelDielectric(etaExt->getRaw(), etaInt->getRaw());
    }
    
    DiffuseReflection::DiffuseReflection(const SpectrumTextureRef &reflectance, const FloatTextureRef &sigma) :
    m_reflectance(reflectance), m_sigma(sigma) {
        m_rawData = new SLR::DiffuseReflection(reflectance->getRaw(), m_sigma ? sigma->getRaw() : nullptr);
    }
    
    SpecularReflection::SpecularReflection(const SpectrumTextureRef &coeffR, const SpatialFresnelRef &fresnel) :
    m_coeffR(coeffR), m_fresnel(fresnel) {
        m_rawData = new SLR::SpecularReflection(coeffR->getRaw(), fresnel->getRaw());
    }
    
    SpecularTransmission::SpecularTransmission(const SpectrumTextureRef &coeffT, const SpectrumTextureRef &etaExt, const SpectrumTextureRef &etaInt) :
    m_coeffT(coeffT), m_etaExt(etaExt), m_etaInt(etaInt) {
        m_rawData = new SLR::SpecularTransmission(coeffT->getRaw(), etaExt->getRaw(), etaInt->getRaw());
    }
    
    AshikhminSpecularReflection::AshikhminSpecularReflection(const SpectrumTextureRef &Rs, const FloatTextureRef &nu, const FloatTextureRef &nv) :
    m_Rs(Rs), m_nu(nu), m_nv(nv) {
        m_rawData = new SLR::AshikhminSpecularReflection(Rs->getRaw(), nu->getRaw(), nv->getRaw());
    }
    
    AshikhminDiffuseReflection::AshikhminDiffuseReflection(const SpectrumTextureRef &Rs, const SpectrumTextureRef &Rd) :
    m_Rs(Rs), m_Rd(Rd) {
        m_rawData = new SLR::AshikhminDiffuseReflection(Rs->getRaw(), Rd->getRaw());
    }
    
    ModifiedWardDurReflection::ModifiedWardDurReflection(const SpectrumTextureRef &reflectance, const FloatTextureRef &anisoX, const FloatTextureRef &anisoY) :
    m_reflectance(reflectance), m_anisoX(anisoX), m_anisoY(anisoY) {
        m_rawData = new SLR::ModifiedWardDurReflection(reflectance->getRaw(), anisoX->getRaw(), anisoY->getRaw());
    }
    
    
    SurfaceMaterialRef SurfaceMaterial::createMatte(const SpectrumTextureRef &reflectance, const FloatTextureRef &sigma) {
        return createShared<DiffuseReflection>(reflectance, sigma);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createMetal(const SpectrumTextureRef &coeffR, const SpectrumTextureRef &eta, const SpectrumTextureRef &k) {
        SpatialFresnelRef fresnel = createShared<SpatialFresnelConductor>(eta, k);
        return createShared<SpecularReflection>(coeffR, fresnel);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createGlass(const SpectrumTextureRef &coeffR, const SpectrumTextureRef &coeffT, const SpectrumTextureRef &etaExt, const SpectrumTextureRef &etaInt) {
        SpatialFresnelRef frDiel = createShared<SpatialFresnelDielectric>(etaExt, etaInt);
        SurfaceMaterialRef r = createShared<SpecularReflection>(coeffR, frDiel);
        SurfaceMaterialRef t = createShared<SpecularTransmission>(coeffT, etaExt, etaInt);
        return createSummedMaterial(r, t);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createModifiedWardDur(const SpectrumTextureRef &reflectance, const FloatTextureRef &anisoX, const FloatTextureRef &anisoY) {
        return createShared<ModifiedWardDurReflection>(reflectance, anisoX, anisoY);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createAshikhminShirley(const SpectrumTextureRef &Rd, const SpectrumTextureRef &Rs, const FloatTextureRef &nu, const FloatTextureRef &nv) {
        SurfaceMaterialRef specular = createShared<AshikhminSpecularReflection>(Rs, nu, nv);
        SurfaceMaterialRef diffuse = createShared<AshikhminDiffuseReflection>(Rs, Rd);
        return createSummedMaterial(specular, diffuse);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createSummedMaterial(const SurfaceMaterialRef &mat0, const SurfaceMaterialRef &mat1) {
        return createShared<SummedSurfaceMaterial>(mat0, mat1);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createMixedMaterial(const SurfaceMaterialRef &mat0, const SurfaceMaterialRef &mat1, const FloatTextureRef &factor) {
        return createShared<MixedSurfaceMaterial>(mat0, mat1, factor);
    }
    
    SurfaceMaterialRef SurfaceMaterial::createEmitterSurfaceMaterial(const SurfaceMaterialRef &mat, const EmitterSurfacePropertyRef &emit) {
        return createShared<EmitterSurfaceMaterial>(mat, emit);
    }
    
    EmitterSurfacePropertyRef SurfaceMaterial::createDiffuseEmitter(const SpectrumTextureRef &emittance) {
        return createShared<DiffuseEmission>(emittance);
    }
}

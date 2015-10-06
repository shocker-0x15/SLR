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
        m_rawData = new SLR::DiffuseReflection(reflectance->getRaw(), sigma->getRaw());
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
}

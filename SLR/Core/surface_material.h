//
//  surface_material.h
//
//  Created by 渡部 心 on 2015/04/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__materials__
#define __SLR__materials__

#include "../defines.h"
#include "../references.h"
#include "../BasicTypes/Spectrum.h"

class SurfaceMaterial {
public:
    SurfaceMaterial() { };
    virtual ~SurfaceMaterial() { };
    
    virtual BSDF* getBSDF(const SurfacePoint &surfPt, ArenaAllocator &mem, float scale = 1.0f) const = 0;
    virtual Spectrum emittance(const SurfacePoint &surfPt) const { return Spectrum::Zero; };
    virtual bool isEmitting() const { return false; };
    virtual EDF* getEDF(const SurfacePoint &surfPt, ArenaAllocator &mem, float scale = 1.0f) const { SLRAssert(false, "Not implemented."); return nullptr; };
    
    static SurfaceMaterialRef createMatte(const SpectrumTextureRef &reflectance, const FloatTextureRef &sigma);
    static SurfaceMaterialRef createMetal(const SpectrumTextureRef &coeffR, const SpectrumTextureRef &eta, const SpectrumTextureRef &k);
    static SurfaceMaterialRef createGlass(const SpectrumTextureRef &coeffR, const SpectrumTextureRef &coeffT, const SpatialFresnelDielectricRef &frDiel);
    static SurfaceMaterialRef createModifiedWardDur(const SpectrumTextureRef &reflectance, const FloatTextureRef &anisoX, const FloatTextureRef &anisoY);
    static SurfaceMaterialRef createAshikhminShirley(const SpectrumTextureRef &Rd, const SpectrumTextureRef &Rs, const FloatTextureRef &nx, const FloatTextureRef &ny);
    static SurfaceMaterialRef createAddedMaterial(const SurfaceMaterialRef &mat0, const SurfaceMaterialRef &mat1);
    static SurfaceMaterialRef createMixedMaterial(const SurfaceMaterialRef &mat0, const SurfaceMaterialRef &mat1, const FloatTextureRef &factor);
    static SurfaceMaterialRef createEmitterSurfaceMaterial(const SurfaceMaterialRef &mat, const EmitterSurfacePropertyRef &emit);
};

class EmitterSurfaceProperty {
public:
    virtual Spectrum emittance(const SurfacePoint &surfPt) const = 0;
    virtual EDF* getEDF(const SurfacePoint &surfPt, ArenaAllocator &mem, float scale = 1.0f) const = 0;
    
    static EmitterSurfacePropertyRef createDiffuseEmitter(const SpectrumTextureRef &emittance);
};

class EmitterSurfaceMaterial : public SurfaceMaterial {
    SurfaceMaterialRef m_mat;
    EmitterSurfacePropertyRef m_emit;
public:
    EmitterSurfaceMaterial(const SurfaceMaterialRef &mat, const EmitterSurfacePropertyRef &emit) :
    m_mat(mat), m_emit(emit) {};
    
    BSDF* getBSDF(const SurfacePoint &surfPt, ArenaAllocator &mem, float scale = 1.0f) const override;
    Spectrum emittance(const SurfacePoint &surfPt) const override;
    bool isEmitting() const { return true; };
    EDF* getEDF(const SurfacePoint &surfPt, ArenaAllocator &mem, float scale) const override;
};


class SpatialFresnel {
public:
    SpatialFresnel() { };
    virtual ~SpatialFresnel() { };
    
    virtual Fresnel* getFresnel(const SurfacePoint &surfPt, ArenaAllocator &mem) const = 0;
};

class SpatialFresnelNoOp : public SpatialFresnel {
public:
    SpatialFresnelNoOp() { };
    
    Fresnel* getFresnel(const SurfacePoint &surfPt, ArenaAllocator &mem) const override;
};

class SpatialFresnelConductor : public SpatialFresnel {
    SpectrumTextureRef m_eta;
    SpectrumTextureRef m_k;
public:
    SpatialFresnelConductor(const SpectrumTextureRef &eta, const SpectrumTextureRef &k) : m_eta(eta), m_k(k) { };
    
    Fresnel* getFresnel(const SurfacePoint &surfPt, ArenaAllocator &mem) const override;
};

class SpatialFresnelDielectric : public SpatialFresnel {
    FloatTextureRef m_etaExt, m_etaInt;
public:
    SpatialFresnelDielectric(const FloatTextureRef &etaExt, const FloatTextureRef &etaInt) : m_etaExt(etaExt), m_etaInt(etaInt) { };
    
    const FloatTextureRef &etaExt() const { return m_etaExt; };
    const FloatTextureRef &etaInt() const { return m_etaInt; };
    
    Fresnel* getFresnel(const SurfacePoint &surfPt, ArenaAllocator &mem) const override;
};

#endif

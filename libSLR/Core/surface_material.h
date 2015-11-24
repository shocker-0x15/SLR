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

namespace SLR {
    class SurfaceMaterial {
    public:
        SurfaceMaterial() { };
        virtual ~SurfaceMaterial() { };
        
        virtual BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const = 0;
        virtual SampledSpectrum emittance(const SurfacePoint &surfPt, const WavelengthSamples &wls) const { return SampledSpectrum::Zero; };
        virtual bool isEmitting() const { return false; };
        virtual EDF* getEDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const { SLRAssert(false, "Not implemented."); return nullptr; };
    };
    
    class EmitterSurfaceProperty {
    public:
        virtual ~EmitterSurfaceProperty() { };
        
        virtual SampledSpectrum emittance(const SurfacePoint &surfPt, const WavelengthSamples &wls) const = 0;
        virtual EDF* getEDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const = 0;
    };
    
    class EmitterSurfaceMaterial : public SurfaceMaterial {
        const SurfaceMaterial* m_mat;
        const EmitterSurfaceProperty* m_emit;
    public:
        EmitterSurfaceMaterial(const SurfaceMaterial* mat, const EmitterSurfaceProperty* emit) :
        m_mat(mat), m_emit(emit) {};
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
        SampledSpectrum emittance(const SurfacePoint &surfPt, const WavelengthSamples &wls) const override;
        bool isEmitting() const override { return true; };
        EDF* getEDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const override;
    };
    
    
    class SpatialFresnel {
    public:
        SpatialFresnel() { };
        virtual ~SpatialFresnel() { };
        
        virtual Fresnel* getFresnel(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const = 0;
    };
    
    class SpatialFresnelNoOp : public SpatialFresnel {
    public:
        SpatialFresnelNoOp() { };
        
        Fresnel* getFresnel(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const override;
    };
    
    class SpatialFresnelConductor : public SpatialFresnel {
        const SpectrumTexture* m_eta;
        const SpectrumTexture* m_k;
    public:
        SpatialFresnelConductor(const SpectrumTexture* eta, const SpectrumTexture* k) : m_eta(eta), m_k(k) { };
        
        Fresnel* getFresnel(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const override;
    };
    
    class SpatialFresnelDielectric : public SpatialFresnel {
        const SpectrumTexture* m_etaExt;
        const SpectrumTexture* m_etaInt;
    public:
        SpatialFresnelDielectric(const SpectrumTexture* etaExt, const SpectrumTexture* etaInt) : m_etaExt(etaExt), m_etaInt(etaInt) { };
        
        const SpectrumTexture* etaExt() const { return m_etaExt; };
        const SpectrumTexture* etaInt() const { return m_etaInt; };
        
        Fresnel* getFresnel(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const override;
    };
}

#endif

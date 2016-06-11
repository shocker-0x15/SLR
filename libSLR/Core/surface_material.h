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
#include "../BasicTypes/RGBTypes.h"
#include "../BasicTypes/SpectrumTypes.h"

namespace SLR {
    class SLR_API SurfaceMaterial {
    public:
        SurfaceMaterial() { }
        virtual ~SurfaceMaterial() { }
        
        virtual BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const = 0;
        virtual SampledSpectrum emittance(const SurfacePoint &surfPt, const WavelengthSamples &wls) const { return SampledSpectrum::Zero; }
        virtual bool isEmitting() const { return false; }
        virtual EDF* getEDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const { SLRAssert(false, "Not implemented."); return nullptr; }
        virtual BSSRDF* getBSSRDF(bool lowerHemisphere, const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const { return nullptr; }
    };
    
    
    
//    class SLR_API AlphaMaskedSurfaceMaterial : public SurfaceMaterial {
//        const FloatTexture* m_alpha;
//    public:
//        AlphaMaskedSurfaceMaterial() { }
//        virtual ~AlphaMaskedSurfaceMaterial() { }
//        
//        virtual BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const = 0;
//        virtual SampledSpectrum emittance(const SurfacePoint &surfPt, const WavelengthSamples &wls) const { return SampledSpectrum::Zero; }
//        virtual bool isEmitting() const { return false; }
//        virtual EDF* getEDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const { SLRAssert(false, "Not implemented."); return nullptr; }
//        virtual BSSRDF* getBSSRDF(bool lowerHemisphere, const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const { return nullptr; }
//    };
    
    
    
    class SLR_API EmitterSurfaceProperty {
    public:
        virtual ~EmitterSurfaceProperty() { }
        
        virtual SampledSpectrum emittance(const SurfacePoint &surfPt, const WavelengthSamples &wls) const = 0;
        virtual EDF* getEDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const = 0;
    };
    
    
    
    class SLR_API EmitterSurfaceMaterial : public SurfaceMaterial {
        const SurfaceMaterial* m_mat;
        const EmitterSurfaceProperty* m_emit;
    public:
        EmitterSurfaceMaterial(const SurfaceMaterial* mat, const EmitterSurfaceProperty* emit) :
        m_mat(mat), m_emit(emit) {
            if (m_mat)
                SLRAssert(!m_mat->isEmitting(), "EmitterSurfaceMaterial can not have an emitting SurfaceMaterial.");
        }
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override { return m_mat->getBSDF(surfPt, wls, mem); }
        SampledSpectrum emittance(const SurfacePoint &surfPt, const WavelengthSamples &wls) const override { return m_emit->emittance(surfPt, wls); }
        bool isEmitting() const override { return true; }
        EDF* getEDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const override { return m_emit->getEDF(surfPt, wls, mem); }
        BSSRDF* getBSSRDF(bool lowerHemisphere, const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const override {
            return m_mat->getBSSRDF(lowerHemisphere, surfPt, wls, mem);
        }
    };
    
    
    
    class SLR_API SubSurfaceScatteringSurfaceMaterial : public SurfaceMaterial {
        SpectrumTexture* m_etaExtTex;
        SpectrumTexture* m_etaIntTex;
        FloatTexture* m_alpha_g;
        SVMicrofacetDistribution* m_D;
        SurfaceMaterial* m_surfMat;
        const InputSpectrum* m_l_sigma_a;
        const InputSpectrum* m_l_sigma_s;
        float m_l_g;
        InputSpectrum* m_innerHHReflectance;
        const InputSpectrum* m_u_sigma_a;
        const InputSpectrum* m_u_sigma_s;
        float m_u_g;
        InputSpectrum* m_outerHHReflectance;
    public:
        SubSurfaceScatteringSurfaceMaterial(const InputSpectrum* etaExt, const InputSpectrum* etaInt, float alpha_g,
                                            const InputSpectrum* l_sigma_a, const InputSpectrum* l_sigma_s, float l_g,
                                            const InputSpectrum* u_sigma_a, const InputSpectrum* u_sigma_s, float u_g);
        ~SubSurfaceScatteringSurfaceMaterial();
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override { return m_surfMat->getBSDF(surfPt, wls, mem); }
        SampledSpectrum emittance(const SurfacePoint &surfPt, const WavelengthSamples &wls) const override { return m_surfMat->emittance(surfPt, wls); }
        bool isEmitting() const override { return false; }
        EDF* getEDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const override { return m_surfMat->getEDF(surfPt, wls, mem); }
        BSSRDF* getBSSRDF(bool lowerHemisphere, const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const override;
    };
    
    
    
    class SLR_API SVFresnel {
    public:
        SVFresnel() { }
        virtual ~SVFresnel() { }
        
        virtual Fresnel* getFresnel(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const = 0;
    };
    
    class SLR_API SVFresnelNoOp : public SVFresnel {
    public:
        SVFresnelNoOp() { }
        
        Fresnel* getFresnel(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const override;
    };
    
    class SLR_API SVFresnelConductor : public SVFresnel {
        const SpectrumTexture* m_eta;
        const SpectrumTexture* m_k;
    public:
        SVFresnelConductor(const SpectrumTexture* eta, const SpectrumTexture* k) : m_eta(eta), m_k(k) { }
        
        Fresnel* getFresnel(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const override;
    };
    
    class SLR_API SVFresnelDielectric : public SVFresnel {
        const SpectrumTexture* m_etaExt;
        const SpectrumTexture* m_etaInt;
    public:
        SVFresnelDielectric(const SpectrumTexture* etaExt, const SpectrumTexture* etaInt) : m_etaExt(etaExt), m_etaInt(etaInt) { };
        
        const SpectrumTexture* etaExt() const { return m_etaExt; }
        const SpectrumTexture* etaInt() const { return m_etaInt; }
        
        Fresnel* getFresnel(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const override;
    };
    
    
    
    class SLR_API SVMicrofacetDistribution {
    public:
        SVMicrofacetDistribution() { }
        virtual ~SVMicrofacetDistribution() { }
        
        virtual MicrofacetDistribution* getMicrofacetDistribution(const SurfacePoint &surfPt, ArenaAllocator &mem) const = 0;
    };
    
    class SLR_API SVGGX : public SVMicrofacetDistribution {
        const FloatTexture* m_alpha_g;
    public:
        SVGGX(const FloatTexture* alpha_g) : m_alpha_g(alpha_g) { }
        
        MicrofacetDistribution* getMicrofacetDistribution(const SurfacePoint &surfPt, ArenaAllocator &mem) const override;
    };
}

#endif

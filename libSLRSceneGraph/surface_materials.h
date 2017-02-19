//
//  surface_materials.h
//
//  Created by 渡部 心 on 2015/10/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLRSceneGraph_surface_materials__
#define __SLRSceneGraph_surface_materials__

#include <libSLR/defines.h>
#include "declarations.h"

namespace SLRSceneGraph {
    class SLR_SCENEGRAPH_API SurfaceMaterial {
    protected:
        SLR::SurfaceMaterial* m_rawData;
    public:
        virtual ~SurfaceMaterial();
        const SLR::SurfaceMaterial* getRaw() const {
            return m_rawData;
        };
        
        static SurfaceMaterialRef createMatte(const SpectrumTextureRef &reflectance, const FloatTextureRef &sigma);
        static SurfaceMaterialRef createMetal(const SpectrumTextureRef &coeffR, const SpectrumTextureRef &eta, const SpectrumTextureRef &k);
        static SurfaceMaterialRef createGlass(const SpectrumTextureRef &coeff, const SpectrumTextureRef &etaExt, const SpectrumTextureRef &etaInt);
        static SurfaceMaterialRef createModifiedWardDur(const SpectrumTextureRef &reflectance, const FloatTextureRef &anisoX, const FloatTextureRef &anisoY);
        static SurfaceMaterialRef createAshikhminShirley(const SpectrumTextureRef &Rd, const SpectrumTextureRef &Rs, const FloatTextureRef &nx, const FloatTextureRef &ny);
        static SurfaceMaterialRef createMicrofacetMetal(const SpectrumTextureRef &eta, const SpectrumTextureRef &k, const FloatTextureRef &alpha_g);
        static SurfaceMaterialRef createMicrofacetGlass(const SpectrumTextureRef &etaExt, const SpectrumTextureRef &etaInt, const FloatTextureRef &alpha_g);
        static SurfaceMaterialRef createFlippedMaterial(const SurfaceMaterialRef &baseMat);
        static SurfaceMaterialRef createSummedMaterial(const SurfaceMaterialRef &mat0, const SurfaceMaterialRef &mat1);
        static SurfaceMaterialRef createMixedMaterial(const SurfaceMaterialRef &mat0, const SurfaceMaterialRef &mat1, const FloatTextureRef &factor);
        static EmitterSurfacePropertyRef createDiffuseEmitter(const SpectrumTextureRef &emittance);
        static EmitterSurfacePropertyRef createIdealDirectionalEmitter(const SpectrumTextureRef &emittance);
        static SurfaceMaterialRef createEmitterSurfaceMaterial(const SurfaceMaterialRef &mat, const EmitterSurfacePropertyRef &emit);
    };
    
    
    
    class SLR_SCENEGRAPH_API EmitterSurfaceProperty {
    protected:
        SLR::EmitterSurfaceProperty* m_rawData;
    public:
        virtual ~EmitterSurfaceProperty();
        const SLR::EmitterSurfaceProperty* getRaw() const {
            return m_rawData;
        };
    };
    
    
    
    class SLR_SCENEGRAPH_API EmitterSurfaceMaterial : public SurfaceMaterial {
        SurfaceMaterialRef m_mat;
        EmitterSurfacePropertyRef m_emit;
    public:
        EmitterSurfaceMaterial(const SurfaceMaterialRef &mat, const EmitterSurfacePropertyRef &emit);
    };
    
    
    
    class SLR_SCENEGRAPH_API SVFresnel {
    protected:
        SLR::SVFresnel* m_rawData;
    public:
        virtual ~SVFresnel();
        const SLR::SVFresnel* getRaw() const {
            return m_rawData;
        };
    };
    
    class SLR_SCENEGRAPH_API SVFresnelNoOp : public SVFresnel {
    public:
        SVFresnelNoOp();
    };
    
    class SLR_SCENEGRAPH_API SVFresnelConductor : public SVFresnel {
        SpectrumTextureRef m_eta;
        SpectrumTextureRef m_k;
    public:
        SVFresnelConductor(const SpectrumTextureRef &eta, const SpectrumTextureRef &k);
    };
    
    class SLR_SCENEGRAPH_API SVFresnelDielectric : public SVFresnel {
        SpectrumTextureRef m_etaExt;
        SpectrumTextureRef m_etaInt;
    public:
        SVFresnelDielectric(const SpectrumTextureRef &etaExt, const SpectrumTextureRef &etaInt);
    };
    
    
    
    class SLR_SCENEGRAPH_API SVMicrofacetDistribution {
    protected:
        SLR::SVMicrofacetDistribution* m_rawData;
    public:
        virtual ~SVMicrofacetDistribution();
        const SLR::SVMicrofacetDistribution* getRaw() const {
            return m_rawData;
        };
    };
    
    class SLR_SCENEGRAPH_API SVGGX : public SVMicrofacetDistribution {
        FloatTextureRef m_alpha_g;
    public:
        SVGGX(const FloatTextureRef &alpha_g);
    };
    
    
    
    class SLR_SCENEGRAPH_API DiffuseReflectionSurfaceMaterial : public SurfaceMaterial {
        SpectrumTextureRef m_reflectance;
        FloatTextureRef m_sigma;
    public:
        DiffuseReflectionSurfaceMaterial(const SpectrumTextureRef &reflectance, const FloatTextureRef &sigma);
    };
    
    class SLR_SCENEGRAPH_API SpecularReflectionSurfaceMaterial : public SurfaceMaterial {
        SpectrumTextureRef m_coeffR;
        SpectrumTextureRef m_eta;
        SpectrumTextureRef m_k;
    public:
        SpecularReflectionSurfaceMaterial(const SpectrumTextureRef &coeffR, const SpectrumTextureRef &eta, const SpectrumTextureRef &k);
    };
    
    class SLR_SCENEGRAPH_API SpecularScatteringSurfaceMaterial : public SurfaceMaterial {
        SpectrumTextureRef m_coeff;
        SpectrumTextureRef m_etaExt;
        SpectrumTextureRef m_etaInt;
    public:
        SpecularScatteringSurfaceMaterial(const SpectrumTextureRef &coeff, const SpectrumTextureRef &etaExt, const SpectrumTextureRef &etaInt);
    };
    
    class SLR_SCENEGRAPH_API FlippedSurfaceMaterial : public SurfaceMaterial {
        SurfaceMaterialRef m_baseMat;
    public:
        FlippedSurfaceMaterial(const SurfaceMaterialRef &baseMat);
    };
    
    class SLR_SCENEGRAPH_API AshikhminShirleyReflectionSurfaceMaterial : public SurfaceMaterial {
        SpectrumTextureRef m_Rs;
        SpectrumTextureRef m_Rd;
        FloatTextureRef m_nu;
        FloatTextureRef m_nv;
    public:
        AshikhminShirleyReflectionSurfaceMaterial(const SpectrumTextureRef &Rs, const SpectrumTextureRef &Rd, const FloatTextureRef &nu, const FloatTextureRef &nv);
    };
    
    class SLR_SCENEGRAPH_API ModifiedWardDurReflectionSurfaceMaterial : public SurfaceMaterial {
        SpectrumTextureRef m_reflectance;
        FloatTextureRef m_anisoX;
        FloatTextureRef m_anisoY;
    public:
        ModifiedWardDurReflectionSurfaceMaterial(const SpectrumTextureRef &reflectance, const FloatTextureRef &anisoX, const FloatTextureRef &anisoY);
    };
    
    class SLR_SCENEGRAPH_API MicrofacetReflectionSurfaceMaterial : public SurfaceMaterial {
        SpectrumTextureRef m_eta;
        SpectrumTextureRef m_k;
        const SVMicrofacetDistributionRef m_D;
    public:
        MicrofacetReflectionSurfaceMaterial(const SpectrumTextureRef &eta, const SpectrumTextureRef &k, const SVMicrofacetDistributionRef &D);
    };
    
    class SLR_SCENEGRAPH_API MicrofacetScatteringSurfaceMaterial : public SurfaceMaterial {
        SpectrumTextureRef m_etaExt;
        SpectrumTextureRef m_etaInt;
        SVMicrofacetDistributionRef m_D;
    public:
        MicrofacetScatteringSurfaceMaterial(const SpectrumTextureRef &etaExt, const SpectrumTextureRef &etaInt, const SVMicrofacetDistributionRef &D);
    };
    
    class SLR_SCENEGRAPH_API SummedSurfaceMaterial : public SurfaceMaterial {
        SurfaceMaterialRef m_mat0;
        SurfaceMaterialRef m_mat1;
    public:
        SummedSurfaceMaterial(const SurfaceMaterialRef &m0, const SurfaceMaterialRef &m1);
    };
    
    class SLR_SCENEGRAPH_API MixedSurfaceMaterial : public SurfaceMaterial {
        SurfaceMaterialRef m_mat0;
        SurfaceMaterialRef m_mat1;
        FloatTextureRef m_factor;
    public:
        MixedSurfaceMaterial(const SurfaceMaterialRef &m0, const SurfaceMaterialRef &m1, const FloatTextureRef &factor);
    };
    
    
    
    class SLR_SCENEGRAPH_API DiffuseEmitterSurfaceProperty : public EmitterSurfaceProperty {
        SpectrumTextureRef m_emittance;
    public:
        DiffuseEmitterSurfaceProperty(const SpectrumTextureRef &emittance);
    };
    
    class SLR_SCENEGRAPH_API IdealDirectionalEmitterSurfaceProperty : public EmitterSurfaceProperty {
        SpectrumTextureRef m_emittance;
    public:
        IdealDirectionalEmitterSurfaceProperty(const SpectrumTextureRef &emittance);
    };
    
    class SLR_SCENEGRAPH_API IBLEmitterSurfaceProperty : public EmitterSurfaceProperty {
        SceneWRef m_scene;
        SpectrumTextureRef m_coeffM;
        float m_scale;
    public:
        IBLEmitterSurfaceProperty(const SceneWRef &scene, const SpectrumTextureRef &coeffM, float scale);
    };
}

#endif /* __SLRSceneGraph_surface_materials__ */

//
//  surface_materials.hpp
//
//  Created by 渡部 心 on 2015/10/06.
//  Copyright © 2015年 渡部 心. All rights reserved.
//

#ifndef surface_materials_hpp
#define surface_materials_hpp

#include <libSLR/defines.h>
#include "references.h"

namespace SLRSceneGraph {
    class SurfaceMaterial {
    protected:
        SLR::SurfaceMaterial* m_rawData;
    public:
        virtual ~SurfaceMaterial();
        const SLR::SurfaceMaterial* getRaw() const {
            return m_rawData;
        };
        
        static SurfaceMaterialRef createMatte(const SpectrumTextureRef &reflectance, const FloatTextureRef &sigma);
        static SurfaceMaterialRef createMetal(const SpectrumTextureRef &coeffR, const SpectrumTextureRef &eta, const SpectrumTextureRef &k);
        static SurfaceMaterialRef createGlass(const SpectrumTextureRef &coeffR, const SpectrumTextureRef &coeffT, const SpectrumTextureRef &etaExt, const SpectrumTextureRef &etaInt);
        static SurfaceMaterialRef createModifiedWardDur(const SpectrumTextureRef &reflectance, const FloatTextureRef &anisoX, const FloatTextureRef &anisoY);
        static SurfaceMaterialRef createAshikhminShirley(const SpectrumTextureRef &Rd, const SpectrumTextureRef &Rs, const FloatTextureRef &nx, const FloatTextureRef &ny);
        static SurfaceMaterialRef createSummedMaterial(const SurfaceMaterialRef &mat0, const SurfaceMaterialRef &mat1);
        static SurfaceMaterialRef createMixedMaterial(const SurfaceMaterialRef &mat0, const SurfaceMaterialRef &mat1, const FloatTextureRef &factor);
        static EmitterSurfacePropertyRef createDiffuseEmitter(const SpectrumTextureRef &emittance);
        static SurfaceMaterialRef createEmitterSurfaceMaterial(const SurfaceMaterialRef &mat, const EmitterSurfacePropertyRef &emit);
    };
    
    class EmitterSurfaceProperty {
    protected:
        SLR::EmitterSurfaceProperty* m_rawData;
    public:
        virtual ~EmitterSurfaceProperty();
        const SLR::EmitterSurfaceProperty* getRaw() const {
            return m_rawData;
        };
    };
    
    class EmitterSurfaceMaterial : public SurfaceMaterial {
        SurfaceMaterialRef m_mat;
        EmitterSurfacePropertyRef m_emit;
    public:
        EmitterSurfaceMaterial(const SurfaceMaterialRef &mat, const EmitterSurfacePropertyRef &emit);
    };
    
    
    class SpatialFresnel {
    protected:
        SLR::SpatialFresnel* m_rawData;
    public:
        virtual ~SpatialFresnel();
        const SLR::SpatialFresnel* getRaw() const {
            return m_rawData;
        };
    };
    
    class SpatialFresnelNoOp : public SpatialFresnel {
    public:
        SpatialFresnelNoOp();
    };
    
    class SpatialFresnelConductor : public SpatialFresnel {
        SpectrumTextureRef m_eta;
        SpectrumTextureRef m_k;
    public:
        SpatialFresnelConductor(const SpectrumTextureRef &eta, const SpectrumTextureRef &k);
    };
    
    class SpatialFresnelDielectric : public SpatialFresnel {
        SpectrumTextureRef m_etaExt;
        SpectrumTextureRef m_etaInt;
    public:
        SpatialFresnelDielectric(const SpectrumTextureRef &etaExt, const SpectrumTextureRef &etaInt);
    };
    
    
    class DiffuseReflection : public SurfaceMaterial {
        SpectrumTextureRef m_reflectance;
        FloatTextureRef m_sigma;
    public:
        DiffuseReflection(const SpectrumTextureRef &reflectance, const FloatTextureRef &sigma);
    };
    
    class SpecularReflection : public SurfaceMaterial {
        SpectrumTextureRef m_coeffR;
        SpatialFresnelRef m_fresnel;
    public:
        SpecularReflection(const SpectrumTextureRef &coeffR, const SpatialFresnelRef &fresnel);
    };
    
    class SpecularTransmission : public SurfaceMaterial {
        SpectrumTextureRef m_coeffT;
        SpectrumTextureRef m_etaExt;
        SpectrumTextureRef m_etaInt;
    public:
        SpecularTransmission(const SpectrumTextureRef &coeffT, const SpectrumTextureRef &etaExt, const SpectrumTextureRef &etaInt);
    };
    
    
    class AshikhminSpecularReflection : public SurfaceMaterial {
        SpectrumTextureRef m_Rs;
        FloatTextureRef m_nu;
        FloatTextureRef m_nv;
    public:
        AshikhminSpecularReflection(const SpectrumTextureRef &Rs, const FloatTextureRef &nu, const FloatTextureRef &nv);
    };
    
    class AshikhminDiffuseReflection : public SurfaceMaterial {
        SpectrumTextureRef m_Rs;
        SpectrumTextureRef m_Rd;
    public:
        AshikhminDiffuseReflection(const SpectrumTextureRef &Rs, const SpectrumTextureRef &Rd);
    };
    
    
    class ModifiedWardDurReflection : public SurfaceMaterial {
        SpectrumTextureRef m_reflectance;
        FloatTextureRef m_anisoX;
        FloatTextureRef m_anisoY;
    public:
        ModifiedWardDurReflection(const SpectrumTextureRef &reflectance, const FloatTextureRef &anisoX, const FloatTextureRef &anisoY);
    };
    
    
    class SummedSurfaceMaterial : public SurfaceMaterial {
        SurfaceMaterialRef m_mat0;
        SurfaceMaterialRef m_mat1;
    public:
        SummedSurfaceMaterial(const SurfaceMaterialRef &m0, const SurfaceMaterialRef &m1);
    };
    
    
    class MixedSurfaceMaterial : public SurfaceMaterial {
        SurfaceMaterialRef m_mat0;
        SurfaceMaterialRef m_mat1;
        FloatTextureRef m_factor;
    public:
        MixedSurfaceMaterial(const SurfaceMaterialRef &m0, const SurfaceMaterialRef &m1, const FloatTextureRef &factor);
    };
    
    
    class DiffuseEmission : public EmitterSurfaceProperty {
        SpectrumTextureRef m_emittance;
    public:
        DiffuseEmission(const SpectrumTextureRef &emittance);
    };
    
    
    class IBLEmission : public EmitterSurfaceProperty {
        SceneWRef m_scene;
        SpectrumTextureRef m_coeffM;
        float m_scale;
    public:
        IBLEmission(const SceneWRef &scene, const SpectrumTextureRef &coeffM, float scale);
    };
}

#endif /* surface_materials_hpp */

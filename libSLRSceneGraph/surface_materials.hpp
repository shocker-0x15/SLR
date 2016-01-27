//
//  surface_materials.hpp
//
//  Created by 渡部 心 on 2015/10/06.
//  Copyright c 2015年 渡部 心. All rights reserved.
//

#ifndef surface_materials_hpp
#define surface_materials_hpp

#include <libSLR/defines.h>
#include "references.h"

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
        static SurfaceMaterialRef createGlass(const SpectrumTextureRef &coeffR, const SpectrumTextureRef &coeffT, const SpectrumTextureRef &etaExt, const SpectrumTextureRef &etaInt);
        static SurfaceMaterialRef createModifiedWardDur(const SpectrumTextureRef &reflectance, const FloatTextureRef &anisoX, const FloatTextureRef &anisoY);
        static SurfaceMaterialRef createAshikhminShirley(const SpectrumTextureRef &Rd, const SpectrumTextureRef &Rs, const FloatTextureRef &nx, const FloatTextureRef &ny);
        static SurfaceMaterialRef createSummedMaterial(const SurfaceMaterialRef &mat0, const SurfaceMaterialRef &mat1);
        static SurfaceMaterialRef createMixedMaterial(const SurfaceMaterialRef &mat0, const SurfaceMaterialRef &mat1, const FloatTextureRef &factor);
        static EmitterSurfacePropertyRef createDiffuseEmitter(const SpectrumTextureRef &emittance);
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
    
    
    class SLR_SCENEGRAPH_API SpatialFresnel {
    protected:
        SLR::SpatialFresnel* m_rawData;
    public:
        virtual ~SpatialFresnel();
        const SLR::SpatialFresnel* getRaw() const {
            return m_rawData;
        };
    };
    
    class SLR_SCENEGRAPH_API SpatialFresnelNoOp : public SpatialFresnel {
    public:
        SpatialFresnelNoOp();
    };
    
    class SLR_SCENEGRAPH_API SpatialFresnelConductor : public SpatialFresnel {
        SpectrumTextureRef m_eta;
        SpectrumTextureRef m_k;
    public:
        SpatialFresnelConductor(const SpectrumTextureRef &eta, const SpectrumTextureRef &k);
    };
    
    class SLR_SCENEGRAPH_API SpatialFresnelDielectric : public SpatialFresnel {
        SpectrumTextureRef m_etaExt;
        SpectrumTextureRef m_etaInt;
    public:
        SpatialFresnelDielectric(const SpectrumTextureRef &etaExt, const SpectrumTextureRef &etaInt);
    };
    
    
    class SLR_SCENEGRAPH_API DiffuseReflection : public SurfaceMaterial {
        SpectrumTextureRef m_reflectance;
        FloatTextureRef m_sigma;
    public:
        DiffuseReflection(const SpectrumTextureRef &reflectance, const FloatTextureRef &sigma);
    };
    
    class SLR_SCENEGRAPH_API SpecularReflection : public SurfaceMaterial {
        SpectrumTextureRef m_coeffR;
        SpatialFresnelRef m_fresnel;
    public:
        SpecularReflection(const SpectrumTextureRef &coeffR, const SpatialFresnelRef &fresnel);
    };
    
    class SLR_SCENEGRAPH_API SpecularTransmission : public SurfaceMaterial {
        SpectrumTextureRef m_coeffT;
        SpectrumTextureRef m_etaExt;
        SpectrumTextureRef m_etaInt;
    public:
        SpecularTransmission(const SpectrumTextureRef &coeffT, const SpectrumTextureRef &etaExt, const SpectrumTextureRef &etaInt);
    };
    
    
    class SLR_SCENEGRAPH_API AshikhminSpecularReflection : public SurfaceMaterial {
        SpectrumTextureRef m_Rs;
        FloatTextureRef m_nu;
        FloatTextureRef m_nv;
    public:
        AshikhminSpecularReflection(const SpectrumTextureRef &Rs, const FloatTextureRef &nu, const FloatTextureRef &nv);
    };
    
    class SLR_SCENEGRAPH_API AshikhminDiffuseReflection : public SurfaceMaterial {
        SpectrumTextureRef m_Rs;
        SpectrumTextureRef m_Rd;
    public:
        AshikhminDiffuseReflection(const SpectrumTextureRef &Rs, const SpectrumTextureRef &Rd);
    };
    
    
    class SLR_SCENEGRAPH_API ModifiedWardDurReflection : public SurfaceMaterial {
        SpectrumTextureRef m_reflectance;
        FloatTextureRef m_anisoX;
        FloatTextureRef m_anisoY;
    public:
        ModifiedWardDurReflection(const SpectrumTextureRef &reflectance, const FloatTextureRef &anisoX, const FloatTextureRef &anisoY);
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
    
    
    class SLR_SCENEGRAPH_API DiffuseEmission : public EmitterSurfaceProperty {
        SpectrumTextureRef m_emittance;
    public:
        DiffuseEmission(const SpectrumTextureRef &emittance);
    };
    
    
    class SLR_SCENEGRAPH_API IBLEmission : public EmitterSurfaceProperty {
        SceneWRef m_scene;
        SpectrumTextureRef m_coeffM;
        float m_scale;
    public:
        IBLEmission(const SceneWRef &scene, const SpectrumTextureRef &coeffM, float scale);
    };
}

#endif /* surface_materials_hpp */

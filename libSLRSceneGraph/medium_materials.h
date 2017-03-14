//
//  medium_materials.h
//
//  Created by 渡部 心 on 2017/02/05.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef __SLRSceneGraph_medium_materials__
#define __SLRSceneGraph_medium_materials__

#include <libSLR/defines.h>
#include "declarations.h"

namespace SLRSceneGraph {
    class SLR_SCENEGRAPH_API MediumMaterial {
    protected:
        SLR::MediumMaterial* m_rawData;
    public:
        virtual ~MediumMaterial();
        
        const SLR::MediumMaterial* getRaw() const {
            return m_rawData;
        };
        
        static MediumMaterialRef createIsotropic();
        static MediumMaterialRef createHenyeyGreenstein(const FloatTextureRef &g);
//        static MediumMaterialRef createSchlick(const FloatTextureRef &g);
    };
    
    
    
    class SLR_SCENEGRAPH_API EmitterMediumProperty {
    protected:
        SLR::EmitterMediumProperty* m_rawData;
    public:
        virtual ~EmitterMediumProperty();
        
        const SLR::EmitterMediumProperty* getRaw() const {
            return m_rawData;
        };
    };
    
    
    
    class SLR_SCENEGRAPH_API EmitterMediumMaterial : public MediumMaterial {
        MediumMaterialRef m_mat;
        EmitterMediumPropertyRef m_emit;
    public:
        EmitterMediumMaterial(const MediumMaterialRef &mat, const EmitterMediumPropertyRef &emit);
    };
    
    
    
    class SLR_SCENEGRAPH_API IsotropicScatteringMediumMaterial : public MediumMaterial {
    public:
        IsotropicScatteringMediumMaterial();
    };
    
    
    
    class SLR_SCENEGRAPH_API HenyeyGreensteinScatteringMediumMaterial : public MediumMaterial {
        FloatTextureRef m_g;
    public:
        HenyeyGreensteinScatteringMediumMaterial(const FloatTextureRef &g);
    };
}

#endif /* __SLRSceneGraph_medium_materials__ */

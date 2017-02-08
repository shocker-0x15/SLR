//
//  medium_materials.hpp
//
//  Created by 渡部 心 on 2017/02/05.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef medium_materials_hpp
#define medium_materials_hpp

#include <libSLR/defines.h>
#include "references.h"

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
//        static MediumMaterialRef createHenyeyGreenstein(const FloatTextureRef &g);
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
    
    
    
    class SLR_SCENEGRAPH_API IsotropicScattering : public MediumMaterial {
    public:
        IsotropicScattering();
    };
}

#endif /* medium_materials_hpp */

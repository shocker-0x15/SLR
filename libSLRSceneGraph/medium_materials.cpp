//
//  medium_materials.cpp
//
//  Created by 渡部 心 on 2017/02/05.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "medium_materials.hpp"
#include <libSLR/MediumMaterials/medium_material_headers.h>
#include "textures.hpp"

namespace SLRSceneGraph {
    MediumMaterial::~MediumMaterial() {
        delete m_rawData;
    }
    
    
    
    EmitterMediumProperty::~EmitterMediumProperty() {
        delete m_rawData;
    }
    
    
    
    EmitterMediumMaterial::EmitterMediumMaterial(const MediumMaterialRef &mat, const EmitterMediumPropertyRef &emit) :
    m_mat(mat), m_emit(emit) {
        m_rawData = new SLR::EmitterMediumMaterial(mat->getRaw(), emit->getRaw());
    }
    
    
    
    IsotropicScatteringMediumMaterial::IsotropicScatteringMediumMaterial() {
        m_rawData = new SLR::IsotropicScatteringMediumMaterial();
    }
    
    
    
    MediumMaterialRef MediumMaterial::createIsotropic() {
        return createShared<IsotropicScatteringMediumMaterial>();
    }
}

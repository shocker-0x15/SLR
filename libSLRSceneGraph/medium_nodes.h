//
//  medium_nodes.h
//
//  Created by 渡部 心 on 2017/01/04.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef __SLRSceneGraph_medium_nodes__
#define __SLRSceneGraph_medium_nodes__

#include <libSLR/defines.h>
#include "references.h"
#include "nodes.h"
#include <libSLR/Core/geometry.h>

namespace SLRSceneGraph {
    class SLR_SCENEGRAPH_API HomogeneousMediumNode : public Node {
        SLR::BoundingBox3D m_region;
        InputSpectrumRef m_sigma_s;
        InputSpectrumRef m_sigma_e;
        MediumMaterialRef m_material;
        
        void allocateRawData() override;
        void setupRawData() override;
        void terminateRawData() override;
    public:
        HomogeneousMediumNode(const SLR::BoundingBox3D &region, const InputSpectrumRef &sigma_s, const InputSpectrumRef &sigma_e, const MediumMaterialRef &material);
        
        NodeRef copy() const override;
        
        void prepareForRendering() override;
    };
    
    
    class SLR_SCENEGRAPH_API GridMediumNode : public Node {
        SLR::BoundingBox3D m_region;
        SLR::InputSpectrum** m_sigma_s;
        SLR::InputSpectrum** m_sigma_e;
        uint32_t m_numX, m_numY, m_numZ;
        MediumMaterialRef m_material;
        
        void allocateRawData() override;
        void setupRawData() override;
        void terminateRawData() override;
    public:
        GridMediumNode(const SLR::BoundingBox3D &region, const SLR::InputSpectrum** sigma_s, const SLR::InputSpectrum** sigma_e,
                       uint32_t elemSize, uint32_t numX, uint32_t numY, uint32_t numZ, const MediumMaterialRef &material);
        ~GridMediumNode();
        
        NodeRef copy() const override;
        
        void prepareForRendering() override;
    };
}

#endif /* __SLRSceneGraph_medium_nodes__ */

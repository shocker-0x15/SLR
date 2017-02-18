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
    class SLR_SCENEGRAPH_API HomogeneousMediumNode : public MediumNode {
        SLR::BoundingBox3D m_region;
        AssetSpectrumRef m_sigma_s;
        AssetSpectrumRef m_sigma_e;
        MediumMaterialRef m_material;
        
        void allocateRawData() override;
        void setupRawData() override;
        void terminateRawData() override;
    public:
        HomogeneousMediumNode(const SLR::BoundingBox3D &region, const AssetSpectrumRef &sigma_s, const AssetSpectrumRef &sigma_e, const MediumMaterialRef &material);
        
        NodeRef copy() const override;
        
        void prepareForRendering() override;
    };
    
    
    
    class SLR_SCENEGRAPH_API GridMediumNode : public MediumNode {
        SLR::BoundingBox3D m_region;
        SLR::AssetSpectrum** m_sigma_s;
        SLR::AssetSpectrum** m_sigma_e;
        uint32_t m_numX, m_numY, m_numZ;
        MediumMaterialRef m_material;
        
        void allocateRawData() override;
        void setupRawData() override;
        void terminateRawData() override;
    public:
        GridMediumNode(const SLR::BoundingBox3D &region, const SLR::AssetSpectrum** sigma_s, const SLR::AssetSpectrum** sigma_e,
                       uint32_t elemSize, uint32_t numX, uint32_t numY, uint32_t numZ, const MediumMaterialRef &material);
        ~GridMediumNode();
        
        NodeRef copy() const override;
        
        void prepareForRendering() override;
    };
    
    
    
    class SLR_SCENEGRAPH_API DensityGridMediumNode : public MediumNode {
        SLR::BoundingBox3D m_region;
        AssetSpectrumRef m_base_sigma_s;
        AssetSpectrumRef m_base_sigma_e;
        std::unique_ptr<float[]> m_density_grid;
        uint32_t m_numX, m_numY, m_numZ;
        MediumMaterialRef m_material;
        
        void allocateRawData() override;
        void setupRawData() override;
        void terminateRawData() override;
    public:
        DensityGridMediumNode(const SLR::BoundingBox3D &region, const AssetSpectrumRef &base_sigma_s, const AssetSpectrumRef &base_sigma_e, std::unique_ptr<float[]> &density_grid,
                              uint32_t numX, uint32_t numY, uint32_t numZ, const MediumMaterialRef &material);
        ~DensityGridMediumNode();
        
        NodeRef copy() const override;
        
        void prepareForRendering() override;
    };
}

#endif /* __SLRSceneGraph_medium_nodes__ */

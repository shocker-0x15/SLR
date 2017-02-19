//
//  medium_nodes.h
//
//  Created by 渡部 心 on 2017/01/08.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef __SLR_medium_nodes__
#define __SLR_medium_nodes__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/geometry.h"
#include "node.h"

namespace SLR {
    class SLR_API HomogeneousMediumNode : public MediumNode {
        HomogeneousMediumDistribution* m_medium;
        const MediumMaterial* m_material;
        
        SingleMediumObject* m_obj;
    public:
        HomogeneousMediumNode(const BoundingBox3D &region, const AssetSpectrum* sigma_s, const AssetSpectrum* sigma_e, const MediumMaterial* material);
        ~HomogeneousMediumNode();
        
        bool isDirectlyTransformable() const override { return false; }
        void createRenderingData(Allocator* mem, const Transform* subTF, RenderingData *data) override;
        void destroyRenderingData(Allocator* mem) override;
    };
    
    
    
    class SLR_API DensityGridMediumNode : public MediumNode {
        DensityGridMediumDistribution* m_medium;
        const MediumMaterial* m_material;
        
        SingleMediumObject* m_obj;
    public:
        DensityGridMediumNode(const BoundingBox3D &region, const AssetSpectrum* base_sigma_s, const AssetSpectrum* base_sigma_e, const float* density_grid,
                              uint32_t numX, uint32_t numY, uint32_t numZ, const MediumMaterial* material);
        ~DensityGridMediumNode();
        
        bool isDirectlyTransformable() const override { return false; }
        void createRenderingData(Allocator* mem, const Transform* subTF, RenderingData *data) override;
        void destroyRenderingData(Allocator* mem) override;
    };
    
    
    
    class SLR_API GridMediumNode : public MediumNode {
        GridMediumDistribution* m_medium;
        const MediumMaterial* m_material;
    public:
        GridMediumNode(const BoundingBox3D &region, const AssetSpectrum** sigma_s_grid, const AssetSpectrum** sigma_e_grid,
                       uint32_t numX, uint32_t numY, uint32_t numZ, const MediumMaterial* material);
        ~GridMediumNode();
        
        bool isDirectlyTransformable() const override { return false; }
        void createRenderingData(Allocator* mem, const Transform* subTF, RenderingData *data) override;
        void destroyRenderingData(Allocator* mem) override;
    };
}

#endif /* __SLR_medium_nodes__ */

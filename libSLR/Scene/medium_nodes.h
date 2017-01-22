//
//  medium_nodes.h
//
//  Created by 渡部 心 on 2017/01/08.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef __SLR_medium_nodes_h__
#define __SLR_medium_nodes_h__

#include "../defines.h"
#include "../references.h"
#include "nodes.h"
#include "geometry.h"

namespace SLR {
    class SLR_API HomogeneousMediumNode : public Node {
        HomogeneousMedium* m_medium;
    public:
        HomogeneousMediumNode(const BoundingBox3D &region, const InputSpectrum* sigma_s, const InputSpectrum* sigma_e);
        ~HomogeneousMediumNode();
        
        bool isDirectlyTransformable() const override { return false; }
        void createRenderingData(Allocator* mem, const Transform* subTF, RenderingData *data) override;
        void destroyRenderingData(Allocator* mem) override;
    };
    
    
    
    class SLR_API GridMediumNode : public Node {
        GridMedium* m_medium;
    public:
        GridMediumNode(const BoundingBox3D &region, const InputSpectrum** sigma_s_grid, const InputSpectrum** sigma_e_grid,
                       uint32_t numX, uint32_t numY, uint32_t numZ);
        ~GridMediumNode();
        
        bool isDirectlyTransformable() const override { return false; }
        void createRenderingData(Allocator* mem, const Transform* subTF, RenderingData *data) override;
        void destroyRenderingData(Allocator* mem) override;
    };
}

#endif /* __SLR_medium_nodes_h__ */

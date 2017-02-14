//
//  medium_nodes.cpp
//
//  Created by 渡部 心 on 2017/01/08.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "medium_nodes.h"
#include "../Memory/ArenaAllocator.h"
#include "../Core/medium_material.h"
#include "../Core/MediumObject.h"
#include "../Medium/HomogeneousMedium.h"
#include "../Medium/GridMedium.h"
#include "../Medium/DensityGridMedium.h"

namespace SLR {
    HomogeneousMediumNode::HomogeneousMediumNode(const BoundingBox3D &region, const InputSpectrum* sigma_s, const InputSpectrum* sigma_e, const MediumMaterial* material) :
    m_material(material) {
        m_medium = new HomogeneousMedium(region, sigma_s, sigma_e);
    }
    
    HomogeneousMediumNode::~HomogeneousMediumNode() {
        delete m_medium;
    }
    
    void HomogeneousMediumNode::createRenderingData(Allocator *mem, const Transform *subTF, RenderingData *data) {
        m_obj = mem->create<SingleMediumObject>(m_medium, m_material);
        data->medObjs.push_back(m_obj);
    }
    
    void HomogeneousMediumNode::destroyRenderingData(Allocator *mem) {
        mem->destroy(m_obj);
    }
    
    
    
    GridMediumNode::GridMediumNode(const BoundingBox3D &region, const InputSpectrum** sigma_s_grid, const InputSpectrum** sigma_e_grid,
                                   uint32_t numX, uint32_t numY, uint32_t numZ, const MediumMaterial* material) :
    m_material(material) {
        float maxExtinctionCoefficient = -INFINITY;
        for (int z = 0; z < numZ; ++z) {
            const InputSpectrum* zSlice = sigma_e_grid[z];
            for (int y = 0; y < numY; ++y) {
                for (int x = 0; x < numX; ++x) {
                    const InputSpectrum* spectrum = zSlice + numX * y + x;
                    maxExtinctionCoefficient = std::max(maxExtinctionCoefficient, spectrum->calcBounds());
                }
            }
        }
        m_medium = new GridMedium(region, sigma_s_grid, sigma_e_grid, numX, numY, numZ, maxExtinctionCoefficient);
    }
    
    GridMediumNode::~GridMediumNode() {
        delete m_medium;
    }
    
    void GridMediumNode::createRenderingData(Allocator *mem, const Transform *subTF, RenderingData *data) {
        
    }
    
    void GridMediumNode::destroyRenderingData(Allocator *mem) {
        
    }
    
    
    
    DensityGridMediumNode::DensityGridMediumNode(const BoundingBox3D &region, const InputSpectrum* base_sigma_s, const InputSpectrum* base_sigma_e, const float* density_grid,
                                                 uint32_t numX, uint32_t numY, uint32_t numZ, const MediumMaterial* material) :
    m_material(material) {
        float maxDensity = -INFINITY;
        for (int z = 0; z < numZ; ++z) {
            for (int y = 0; y < numY; ++y) {
                for (int x = 0; x < numX; ++x) {
                    maxDensity = std::max(maxDensity, density_grid[numX * numY * z + numX * y + x]);
                }
            }
        }
        float maxExtinctionCoefficient = base_sigma_e->calcBounds() * maxDensity;
        m_medium = new DensityGridMedium(region, base_sigma_s, base_sigma_e, density_grid, numX, numY, numZ, maxExtinctionCoefficient);
    }
    
    DensityGridMediumNode::~DensityGridMediumNode() {
        delete m_medium;
    }
    
    void DensityGridMediumNode::createRenderingData(Allocator *mem, const Transform *subTF, RenderingData *data) {
        m_obj = mem->create<SingleMediumObject>(m_medium, m_material);
        data->medObjs.push_back(m_obj);
    }
    
    void DensityGridMediumNode::destroyRenderingData(Allocator *mem) {
        mem->destroy(m_obj);
    }
}

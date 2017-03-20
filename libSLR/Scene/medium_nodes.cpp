//
//  medium_nodes.cpp
//
//  Created by 渡部 心 on 2017/01/08.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "medium_nodes.h"

#include "../MemoryAllocators/ArenaAllocator.h"
#include "../Core/medium_object.h"
#include "../Core/medium_material.h"
#include "../MediumDistribution/HomogeneousMediumDistribution.h"
#include "../MediumDistribution/GridMediumDistribution.h"
#include "../MediumDistribution/DensityGridMediumDistribution.h"
#include "../MediumDistribution/VacuumMediumDistribution.h"

namespace SLR {
    HomogeneousMediumNode::HomogeneousMediumNode(const BoundingBox3D &region, const AssetSpectrum* sigma_s, const AssetSpectrum* sigma_e, const MediumMaterial* material) :
    m_material(material) {
        m_medium = new HomogeneousMediumDistribution(region, sigma_s, sigma_e);
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
    
    
    
    DensityGridMediumNode::DensityGridMediumNode(const BoundingBox3D &region, const AssetSpectrum* base_sigma_s, const AssetSpectrum* base_sigma_e, 
                                                 const std::vector<std::vector<float>> &density_grid, uint32_t numX, uint32_t numY, uint32_t numZ, const MediumMaterial* material) :
    m_material(material) {
        m_medium = new DensityGridMediumDistribution(region, base_sigma_s, base_sigma_e, density_grid, numX, numY, numZ);
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
    
    
    
    GridMediumNode::GridMediumNode(const BoundingBox3D &region, const AssetSpectrum** sigma_s_grid, const AssetSpectrum** sigma_e_grid,
                                   uint32_t numX, uint32_t numY, uint32_t numZ, const MediumMaterial* material) :
    m_material(material) {
        m_medium = new GridMediumDistribution(region, sigma_s_grid, sigma_e_grid, numX, numY, numZ);
    }
    
    GridMediumNode::~GridMediumNode() {
        delete m_medium;
    }
    
    void GridMediumNode::createRenderingData(Allocator *mem, const Transform *subTF, RenderingData *data) {
        
    }
    
    void GridMediumNode::destroyRenderingData(Allocator *mem) {
        
    }
    
    
    
    VacuumMediumNode::VacuumMediumNode(const BoundingBox3D &region) {
        m_medium = new VacuumMediumDistribution(region);
        m_material = new NullMediumMaterial();
    }
    
    VacuumMediumNode::~VacuumMediumNode() {
        delete m_material;
        delete m_medium;
    }
    
    void VacuumMediumNode::createRenderingData(Allocator *mem, const Transform *subTF, RenderingData *data) {
        m_obj = mem->create<SingleMediumObject>(m_medium, m_material);
        data->medObjs.push_back(m_obj);
    }
    
    void VacuumMediumNode::destroyRenderingData(Allocator *mem) {
        mem->destroy(m_obj);
    }
}

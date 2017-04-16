//
//  medium_nodes.cpp
//
//  Created by 渡部 心 on 2017/01/04.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "medium_nodes.h"

#include <libSLR/Core/transform.h>
#include <libSLR/Core/medium_object.h>
#include <libSLR/Scene/medium_nodes.h>
#include <libSLR/MediumDistribution/HomogeneousMediumDistribution.h>
#include <libSLR/MediumDistribution/GridMediumDistribution.h>
#include <libSLR/MediumDistribution/DensityGridMediumDistribution.h>
#include <libSLR/MediumDistribution/VacuumMediumDistribution.h>
#include "../medium_materials.h"

namespace SLRSceneGraph {
    void HomogeneousMediumNode::allocateRawData() {
        m_rawData = (SLR::Node*)malloc(sizeof(SLR::HomogeneousMediumNode));
    }
    
    void HomogeneousMediumNode::setupRawData() {
        new (m_rawData) SLR::HomogeneousMediumNode(m_region, m_sigma_s.get(), m_sigma_e.get(), m_material->getRaw());
        m_setup = true;
    }
    
    void HomogeneousMediumNode::terminateRawData() {
        SLR::HomogeneousMediumNode &raw = *(SLR::HomogeneousMediumNode*)m_rawData;
        if (m_setup)
            raw.~HomogeneousMediumNode();
        m_setup = false;
    }
    
    HomogeneousMediumNode::HomogeneousMediumNode(const SLR::BoundingBox3D &region, const AssetSpectrumRef &sigma_s, const AssetSpectrumRef &sigma_e, const MediumMaterialRef &material) :
    m_region(region), m_sigma_s(sigma_s), m_sigma_e(sigma_e), m_material(material) {
        allocateRawData();
    }
    
    NodeRef HomogeneousMediumNode::copy() const {
        return createShared<HomogeneousMediumNode>(m_region, m_sigma_s, m_sigma_e, m_material);
    }
    
    void HomogeneousMediumNode::prepareForRendering() {
        terminateRawData();
        setupRawData();
    }
    
    
    
    void DensityGridMediumNode::allocateRawData() {
        m_rawData = (SLR::Node*)malloc(sizeof(SLR::DensityGridMediumNode));
    }
    
    void DensityGridMediumNode::setupRawData() {
        new (m_rawData) SLR::DensityGridMediumNode(m_region, m_base_sigma_s.get(), m_base_sigma_e.get(), m_density_grid, m_numX, m_numY, m_numZ, m_material->getRaw());
        m_setup = true;
    }
    
    void DensityGridMediumNode::terminateRawData() {
        SLR::DensityGridMediumNode &raw = *(SLR::DensityGridMediumNode*)m_rawData;
        if (m_setup)
            raw.~DensityGridMediumNode();
        m_setup = false;
    }
    
    DensityGridMediumNode::DensityGridMediumNode(const SLR::BoundingBox3D &region, const AssetSpectrumRef &base_sigma_s, const AssetSpectrumRef &base_sigma_e, 
                                                 std::vector<std::vector<float>> &&density_grid, uint32_t numX, uint32_t numY, uint32_t numZ, const MediumMaterialRef &material) :
    m_region(region), m_base_sigma_s(base_sigma_s), m_base_sigma_e(base_sigma_e),
    m_density_grid(density_grid), m_numX(numX), m_numY(numY), m_numZ(numZ), m_material(material) {
        allocateRawData();
    }
    
    NodeRef DensityGridMediumNode::copy() const {
        SLRAssert_NotImplemented();
        return nullptr;
    }
    
    void DensityGridMediumNode::prepareForRendering() {
        terminateRawData();
        setupRawData();
    }
    
    
    
    void GridMediumNode::allocateRawData() {
        m_rawData = (SLR::Node*)malloc(sizeof(SLR::GridMediumNode));
    }
    
    void GridMediumNode::setupRawData() {
        SLRAssert_NotImplemented();
        //        new (m_rawData) SLR::GridMediumNode(m_region, m_sigma_s, m_sigma_e, m_numX, m_numY, m_numZ);
        m_setup = true;
    }
    
    void GridMediumNode::terminateRawData() {
        SLR::GridMediumNode &raw = *(SLR::GridMediumNode*)m_rawData;
        if (m_setup)
            raw.~GridMediumNode();
        m_setup = false;
    }
    
    GridMediumNode::GridMediumNode(const SLR::BoundingBox3D &region, const SLR::AssetSpectrum** sigma_s, const SLR::AssetSpectrum** sigma_e,
                                   uint32_t elemSize, uint32_t numX, uint32_t numY, uint32_t numZ, const MediumMaterialRef &material) :
    m_region(region),
    m_numX(numX), m_numY(numY), m_numZ(numZ), m_material(material) {
        // allocate 3D grid as 2D array.
        m_sigma_s = (SLR::AssetSpectrum**)malloc(sizeof(SLR::AssetSpectrum*) * m_numZ);
        m_sigma_e = (SLR::AssetSpectrum**)malloc(sizeof(SLR::AssetSpectrum*) * m_numZ);
        for (int z = 0; z < m_numZ; ++z) {
            size_t allocateSize = elemSize * m_numX * m_numY;
            m_sigma_s[z] = (SLR::AssetSpectrum*)malloc(allocateSize);
            m_sigma_e[z] = (SLR::AssetSpectrum*)malloc(allocateSize);
            memcpy((void*)m_sigma_s[z], (void*)sigma_s[z], allocateSize);
            memcpy((void*)m_sigma_e[z], (void*)sigma_e[z], allocateSize);
        }
        
        allocateRawData();
    }
    
    GridMediumNode::~GridMediumNode() {
        for (int z = 0; z < m_numZ; ++z) {
            free(m_sigma_s[z]);
            free(m_sigma_e[z]);
        }
        free(m_sigma_s);
        free(m_sigma_e);
        
        Node::~Node();
    }
    
    NodeRef GridMediumNode::copy() const {
        SLRAssert_NotImplemented();
        return nullptr;
    }
    
    void GridMediumNode::prepareForRendering() {
        terminateRawData();
        setupRawData();
    }
    
    
    
    void VacuumMediumNode::allocateRawData() {
        m_rawData = (SLR::Node*)malloc(sizeof(SLR::VacuumMediumNode));
    }
    
    void VacuumMediumNode::setupRawData() {
        new (m_rawData) SLR::VacuumMediumNode(m_region);
        m_setup = true;
    }
    
    void VacuumMediumNode::terminateRawData() {
        SLR::VacuumMediumNode &raw = *(SLR::VacuumMediumNode*)m_rawData;
        if (m_setup)
            raw.~VacuumMediumNode();
        m_setup = false;
    }
    
    VacuumMediumNode::VacuumMediumNode(const SLR::BoundingBox3D &region) :
    m_region(region) {
        allocateRawData();
    }
    
    VacuumMediumNode::~VacuumMediumNode() {
    }
    
    NodeRef VacuumMediumNode::copy() const {
        return createShared<VacuumMediumNode>(m_region);
    }
    
    void VacuumMediumNode::prepareForRendering() {
        terminateRawData();
        setupRawData();
    }
    
    
    
    void CloudMediumNode::allocateRawData() {
        m_rawData = (SLR::Node*)malloc(sizeof(SLR::DensityGridMediumNode));
    }
    
    void CloudMediumNode::setupRawData() {
        new (m_rawData) SLR::CloudMediumNode(m_region, m_featureScale, m_density, m_rngSeed, m_material->getRaw());
        m_setup = true;
    }
    
    void CloudMediumNode::terminateRawData() {
        SLR::CloudMediumNode &raw = *(SLR::CloudMediumNode*)m_rawData;
        if (m_setup)
            raw.~CloudMediumNode();
        m_setup = false;
    }
    
    CloudMediumNode::CloudMediumNode(const SLR::BoundingBox3D &region, float featureScale, float density, uint32_t rngSeed, const MediumMaterialRef &material) :
    m_region(region), m_featureScale(featureScale), m_density(density), m_rngSeed(rngSeed), m_material(material) {
        allocateRawData();
    }
    
    CloudMediumNode::~CloudMediumNode() {
        Node::~Node();
    }
    
    NodeRef CloudMediumNode::copy() const {
        SLRAssert_NotImplemented();
        return nullptr;
    }
    
    void CloudMediumNode::prepareForRendering() {
        terminateRawData();
        setupRawData();
    }
}

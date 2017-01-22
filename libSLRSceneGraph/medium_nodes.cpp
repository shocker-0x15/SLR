//
//  medium_nodes.cpp
//
//  Created by 渡部 心 on 2017/01/04.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "medium_nodes.h"
#include <libSLR/Core/Transform.h>
#include <libSLR/Core/MediumObject.h>
#include <libSLR/Medium/HomogeneousMedium.h>
#include <libSLR/Medium/GridMedium.h>

namespace SLRSceneGraph {
    void HomogeneousMediumNode::setupRawData() {
        
    }
    
    HomogeneousMediumNode::HomogeneousMediumNode(const SLR::BoundingBox3D &region, const InputSpectrumRef &sigma_s, const InputSpectrumRef &sigma_e) :
    m_region(region), m_sigma_s(sigma_s), m_sigma_e(sigma_e) {
        setupRawData();
    }
    
    NodeRef HomogeneousMediumNode::copy() const {
        SLRAssert_NotImplemented();
        return nullptr;
    }
    
    void HomogeneousMediumNode::prepareForRendering() {
        
    }
    
    
    
    void GridMediumNode::setupRawData() {
        
    }
    
    GridMediumNode::GridMediumNode(const SLR::BoundingBox3D &region, const SLR::InputSpectrum** sigma_s, const SLR::InputSpectrum** sigma_e,
                                   uint32_t elemSize, uint32_t numX, uint32_t numY, uint32_t numZ) :
    m_region(region),
    m_numX(numX), m_numY(numY), m_numZ(numZ) {
        // allocate 3D grid as 2D array.
        m_sigma_s = (SLR::InputSpectrum**)malloc(sizeof(SLR::InputSpectrum*) * m_numZ);
        m_sigma_e = (SLR::InputSpectrum**)malloc(sizeof(SLR::InputSpectrum*) * m_numZ);
        for (int z = 0; z < m_numZ; ++z) {
            size_t allocateSize = elemSize * m_numX * m_numY;
            m_sigma_s[z] = (SLR::InputSpectrum*)malloc(allocateSize);
            m_sigma_e[z] = (SLR::InputSpectrum*)malloc(allocateSize);
            memcpy((void*)m_sigma_s[z], (void*)sigma_s[z], allocateSize);
            memcpy((void*)m_sigma_e[z], (void*)sigma_e[z], allocateSize);
        }
        
        setupRawData();
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
        
    }
}

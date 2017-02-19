//
//  textures.cpp
//
//  Created by 渡部 心 on 2015/10/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "textures.h"

#include <libSLR/Texture/texture_headers.h>

namespace SLRSceneGraph {
    Texture2DMapping::Texture2DMapping() {
        m_rawData = new SLR::Texture2DMapping;
    }
    
    Texture2DMapping::~Texture2DMapping() {
        delete m_rawData;
    }
    
    Texture3DMapping::Texture3DMapping() {
        m_rawData = new SLR::Texture3DMapping;
    }
    
    Texture3DMapping::~Texture3DMapping() {
        delete m_rawData;
    }
    
    OffsetAndScale2DMapping::OffsetAndScale2DMapping(float ox, float oy, float sx, float sy) :
    m_offsetX(ox), m_offsetY(oy), m_scaleX(sx), m_scaleY(sy) {
        m_rawData = new SLR::OffsetAndScale2DMapping(ox, oy, sx, sy);
    }
    
    WorldPosition3DMapping::WorldPosition3DMapping() {
        m_rawData = new SLR::WorldPosition3DMapping();
    }
    
    
    
    SpectrumTexture::~SpectrumTexture() {
        delete m_rawData;
    }
    
    NormalTexture::~NormalTexture() {
        delete m_rawData;
    }
    
    FloatTexture::~FloatTexture() {
        delete m_rawData;
    }
    
    ConstantSpectrumTexture::ConstantSpectrumTexture(const AssetSpectrumRef &value) : m_value(value) {
        m_rawData = new SLR::ConstantSpectrumTexture(value.get());
    }
    
    ConstantFloatTexture::ConstantFloatTexture(float value) : m_value(value) {
        m_rawData = new SLR::ConstantFloatTexture(value);
    }
    
    ImageSpectrumTexture::ImageSpectrumTexture(const Texture2DMappingRef &mapping, const TiledImage2DRef &image) :
    m_mapping(mapping), m_data(image) {
        m_rawData = new SLR::ImageSpectrumTexture(image.get(), mapping->getRaw());
    }
    
    ImageNormalTexture::ImageNormalTexture(const Texture2DMappingRef &mapping, const TiledImage2DRef &image) :
    m_mapping(mapping), m_data(image) {
        m_rawData = new SLR::ImageNormalTexture(image.get(), mapping->getRaw());
    }
    
    ImageFloatTexture::ImageFloatTexture(const Texture2DMappingRef &mapping, const TiledImage2DRef &image) :
    m_mapping(mapping), m_data(image) {
        m_rawData = new SLR::ImageFloatTexture(image.get(), mapping->getRaw());
    }
    
    CheckerBoardSpectrumTexture::CheckerBoardSpectrumTexture(const Texture2DMappingRef &mapping, const AssetSpectrumRef &v0, const AssetSpectrumRef &v1) :
    m_mapping(mapping), m_values{v0, v1} {
        m_rawData = new SLR::CheckerBoardSpectrumTexture(mapping->getRaw(), v0.get(), v1.get());
    }
    
    CheckerBoardNormalTexture::CheckerBoardNormalTexture(const Texture2DMappingRef &mapping, float stepWidth, bool reverse) :
    m_mapping(mapping), m_stepWidth(stepWidth), m_reverse(reverse) {
        m_rawData = new SLR::CheckerBoardNormalTexture(mapping->getRaw(), stepWidth, reverse);
    }
    
    CheckerBoardFloatTexture::CheckerBoardFloatTexture(const Texture2DMappingRef &mapping, float v0, float v1) :
    m_mapping(mapping), m_values{v0, v1} {
        m_rawData = new SLR::CheckerBoardFloatTexture(mapping->getRaw(), v0, v1);
    }
    
    VoronoiSpectrumTexture::VoronoiSpectrumTexture(const Texture3DMappingRef &mapping, float scale, float brightness) :
    m_mapping(mapping), m_scale(scale), m_brightness(brightness) {
        m_rawData = new SLR::VoronoiSpectrumTexture(mapping->getRaw(), scale, brightness);
    }
    
    VoronoiNormalTexture::VoronoiNormalTexture(const Texture3DMappingRef &mapping, float scale, float thetaMax) :
    m_mapping(mapping), m_scale(scale), m_thetaMax(thetaMax) {
        m_rawData = new SLR::VoronoiNormalTexture(mapping->getRaw(), scale, thetaMax);
    }
    
    VoronoiFloatTexture::VoronoiFloatTexture(const Texture3DMappingRef &mapping, float scale, float valueScale, bool flat) :
    m_mapping(mapping), m_scale(scale) {
        m_rawData = new SLR::VoronoiFloatTexture(mapping->getRaw(), scale, valueScale, flat);
    }
}

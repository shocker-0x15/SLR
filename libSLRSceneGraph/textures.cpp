//
//  textures.cpp
//
//  Created by 渡部 心 on 2015/10/06.
//  Copyright c 2015年 渡部 心. All rights reserved.
//

#include "textures.hpp"
#include <libSLR/Textures/texture_headers.h>

namespace SLRSceneGraph {
    SpectrumTexture::~SpectrumTexture() {
        delete m_rawData;
    }
    
    Normal3DTexture::~Normal3DTexture() {
        delete m_rawData;
    }
    
    FloatTexture::~FloatTexture() {
        delete m_rawData;
    }
    
    ConstantSpectrumTexture::ConstantSpectrumTexture(const InputSpectrumRef &value) : m_value(value) {
        m_rawData = new SLR::ConstantSpectrumTexture(value.get());
    }
    
    ConstantFloatTexture::ConstantFloatTexture(float value) : m_value(value) {
        m_rawData = new SLR::ConstantFloatTexture(value);
    }
    
    ImageSpectrumTexture::ImageSpectrumTexture(const TiledImage2DRef &image) : m_data(image) {
        m_rawData = new SLR::ImageSpectrumTexture(image.get());
    }
    
    ImageNormal3DTexture::ImageNormal3DTexture(const TiledImage2DRef &image) : m_data(image) {
        m_rawData = new SLR::ImageNormal3DTexture(image.get());
    }
    
    ImageFloatTexture::ImageFloatTexture(const TiledImage2DRef &image) : m_data(image) {
        m_rawData = new SLR::ImageFloatTexture(image.get());
    }
    
    CheckerBoardSpectrumTexture::CheckerBoardSpectrumTexture(const InputSpectrumRef &v0, const InputSpectrumRef &v1) : m_values{v0, v1} {
        m_rawData = new SLR::CheckerBoardSpectrumTexture(v0.get(), v1.get());
    }
    
    CheckerBoardNormal3DTexture::CheckerBoardNormal3DTexture(float stepWidth, bool reverse) : m_stepWidth(stepWidth), m_reverse(reverse) {
        m_rawData = new SLR::CheckerBoardNormal3DTexture(stepWidth, reverse);
    }
    
    CheckerBoardFloatTexture::CheckerBoardFloatTexture(float v0, float v1) : m_values{v0, v1} {
        m_rawData = new SLR::CheckerBoardFloatTexture(v0, v1);
    }
    
    VoronoiSpectrumTexture::VoronoiSpectrumTexture(float scale) : m_scale(scale) {
        m_rawData = new SLR::VoronoiSpectrumTexture(scale);
    }
    
    VoronoiNormal3DTexture::VoronoiNormal3DTexture(float scale) : m_scale(scale) {
        m_rawData = new SLR::VoronoiNormal3DTexture(scale);
    }
    
    VoronoiFloatTexture::VoronoiFloatTexture(float scale) : m_scale(scale) {
        m_rawData = new SLR::VoronoiFloatTexture(scale);
    }
}

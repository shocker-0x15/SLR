//
//  textures.cpp
//
//  Created by 渡部 心 on 2015/10/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "textures.h"

#include <libSLR/Texture/texture_headers.h>

#include "images.h"

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
    
    bool ConstantSpectrumTexture::generateLuminanceChannel() {
        SLR::ConstantSpectrumTexture &raw = *(SLR::ConstantSpectrumTexture*)m_rawData;
        raw.generateLuminanceChannel();
        return true;
    }
    
    ConstantFloatTexture::ConstantFloatTexture(float value) : m_value(value) {
        m_rawData = new SLR::ConstantFloatTexture(value);
    }
    
    ImageSpectrumTexture::ImageSpectrumTexture(const Texture2DMappingRef &mapping, const Image2DRef &image) :
    m_mapping(mapping), m_data(image) {
        m_rawData = new SLR::ImageSpectrumTexture(image->getRaw(), mapping->getRaw());
    }
    
    ImageNormalTexture::ImageNormalTexture(const Texture2DMappingRef &mapping, const Image2DRef &image) :
    m_mapping(mapping), m_data(image) {
        m_rawData = new SLR::ImageNormalTexture(image->getRaw(), mapping->getRaw());
    }
    
    ImageFloatTexture::ImageFloatTexture(const Texture2DMappingRef &mapping, const Image2DRef &image) :
    m_mapping(mapping), m_data(image) {
        m_rawData = new SLR::ImageFloatTexture(image->getRaw(), mapping->getRaw());
    }
    
    CheckerBoardSpectrumTexture::CheckerBoardSpectrumTexture(const Texture2DMappingRef &mapping, const AssetSpectrumRef &v0, const AssetSpectrumRef &v1) :
    m_mapping(mapping), m_values{v0, v1} {
        m_rawData = new SLR::CheckerBoardSpectrumTexture(mapping->getRaw(), v0.get(), v1.get());
    }
    
    bool CheckerBoardSpectrumTexture::generateLuminanceChannel() {
        SLR::CheckerBoardSpectrumTexture &raw = *(SLR::CheckerBoardSpectrumTexture*)m_rawData;
        raw.generateLuminanceChannel();
        return true;
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
    
    WorleyNoiseFloatTexture::WorleyNoiseFloatTexture(const Texture3DMappingRef &mapping, 
                                                     uint32_t numOctaves, float initialFrequency, float supValueOrInitialAmplitude, bool supSpecified, float clipValue, 
                                                     float frequencyMultiplier, float persistence, uint32_t repeat) : 
    m_mapping(mapping), m_numOctaves(numOctaves), 
    m_initialFrequency(initialFrequency), m_supValueOrInitialAmplitude(supValueOrInitialAmplitude), m_supSpecified(supSpecified), m_clipValue(clipValue), 
    m_frequencyMultiplier(frequencyMultiplier), m_persistence(persistence), m_repeat(repeat) {
        m_rawData = new SLR::WorleyNoiseFloatTexture(mapping->getRaw(), numOctaves, initialFrequency, supValueOrInitialAmplitude, supSpecified, clipValue, 
                                                     frequencyMultiplier, persistence, repeat);
    }
    
    PerlinNoiseNormalTexture::PerlinNoiseNormalTexture(const Texture3DMappingRef &mapping, float thetaMax, 
                                                       uint32_t numOctaves, float initialFrequencyPhi, float initialFrequencyTheta, 
                                                       float frequencyMultiplier, float persistence, uint32_t repeat) : 
    m_mapping(mapping), m_thetaMax(thetaMax),  
    m_numOctaves(numOctaves), m_initialFrequencyPhi(initialFrequencyPhi), m_initialFrequencyTheta(initialFrequencyTheta), 
    m_frequencyMultiplier(frequencyMultiplier), m_persistence(persistence), m_repeat(repeat) {
        m_rawData = new SLR::PerlinNoiseNormalTexture(mapping->getRaw(), thetaMax, 
                                                      numOctaves, initialFrequencyPhi, initialFrequencyTheta, 
                                                      frequencyMultiplier, persistence, repeat);
    }
    
    PerlinNoiseFloatTexture::PerlinNoiseFloatTexture(const Texture3DMappingRef &mapping, 
                                                     uint32_t numOctaves, float initialFrequency, float supValueOrInitialAmplitude, bool supSpecified, 
                                                     float frequencyMultiplier, float persistence, uint32_t repeat) : 
    m_mapping(mapping), m_numOctaves(numOctaves), 
    m_initialFrequency(initialFrequency), m_supValueOrInitialAmplitude(supValueOrInitialAmplitude), m_supSpecified(supSpecified), 
    m_frequencyMultiplier(frequencyMultiplier), m_persistence(persistence), m_repeat(repeat) {
        m_rawData = new SLR::PerlinNoiseFloatTexture(mapping->getRaw(), numOctaves, initialFrequency, supValueOrInitialAmplitude, supSpecified, 
                                                     frequencyMultiplier, persistence, repeat);
    }
    
    AnalyticSkySpectrumTexture::AnalyticSkySpectrumTexture(const Texture2DMappingRef &mapping, 
                                                           float solarRadius, float soloarElevation, float turbidity, const AssetSpectrumRef &groundAlbedo, float extAngleOfHorizon) : 
    m_mapping(mapping), m_solarElevation(soloarElevation), m_turbidity(turbidity), m_groundAlbedo(groundAlbedo), m_extAngleOfHorizon(extAngleOfHorizon) {
        m_rawData = new SLR::AnalyticSkySpectrumTexture(solarRadius, soloarElevation, turbidity, groundAlbedo.get(), extAngleOfHorizon, mapping->getRaw());
    }
}

//
//  Image.h
//
//  Created by 渡部 心 on 2015/07/11.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__Image__
#define __SLR__Image__

#include "../defines.h"
#include "../references.h"
#include "../Memory/Allocator.h"
#include "../BasicTypes/Spectrum.h"
#include <half.h>

namespace SLR {
    enum class ColorFormat {
        RGB8x3 = 0,
        RGB_8x4,
        RGBA8x4,
        RGBA16Fx4,
        Gray8,
        uvs16Fx3,
        uvsA16Fx4,
        Num
    };
    
    struct RGB8x3 { uint8_t r, g, b; };
    struct RGB_8x4 { uint8_t r, g, b, dummy; };
    struct RGBA8x4  { uint8_t r, g, b, a; };
    struct RGBA16Fx4 { half r, g, b, a; };
    struct uvs16Fx3 { half u, v, s; };
    struct uvsA16Fx4 { half u, v, s, a; };
    struct Gray8 { uint8_t v; };
    
    extern const size_t sizesOfColorFormats[(uint32_t)ColorFormat::Num];
    
    class Image2D {
    protected:
        uint32_t m_width, m_height;
        SpectrumType m_spType;
        ColorFormat m_colorFormat;
        
        virtual const void* getInternal(uint32_t x, uint32_t y) const = 0;
        virtual void setInternal(uint32_t x, uint32_t y, const void* data, size_t size) = 0;
    public:
        Image2D() { };
        Image2D(uint32_t w, uint32_t h, ColorFormat fmt, SpectrumType spType) : m_width(w), m_height(h), m_colorFormat(fmt), m_spType(spType) { };
        ~Image2D() { };
        
        template <typename ColFmt>
        const ColFmt &get(uint32_t x, uint32_t y) const { return *(const ColFmt*)getInternal(x, y); };
        template <typename ColFmt>
        void set(uint32_t x, uint32_t y, const ColFmt &data) { setInternal(x, y, &data, sizeof(data)); };
        
        void areaAverage(float xLeft, float xRight, float yTop, float yBottom, void* avg) const;
        
        uint32_t width() const { return m_width; };
        uint32_t height() const { return m_height; };
        SpectrumType spectrumType() const { return m_spType; };
        ColorFormat format() const { return m_colorFormat; };
    };
    
    
    template <uint32_t log2_tileWidth>
    class TiledImage2DTemplate : public Image2D {
        static const size_t tileWidth = 1 << log2_tileWidth;
        static const uint32_t localMask = (1 << log2_tileWidth) - 1;
        size_t m_stride;
        size_t m_numTileX;
        size_t m_allocSize;
        uint8_t* m_data;
        
        const void* getInternal(uint32_t x, uint32_t y) const override {
            uint32_t tx = x >> log2_tileWidth;
            uint32_t ty = y >> log2_tileWidth;
            uint32_t lx = x & localMask;
            uint32_t ly = y & localMask;
            return m_data + m_stride * ((ty * m_numTileX + tx) * tileWidth * tileWidth + ly * tileWidth + lx);
        };
        
        void setInternal(uint32_t x, uint32_t y, const void* data, size_t size) override {
            uint32_t tx = x >> log2_tileWidth;
            uint32_t ty = y >> log2_tileWidth;
            uint32_t lx = x & localMask;
            uint32_t ly = y & localMask;
            std::memcpy(m_data + m_stride * ((ty * m_numTileX + tx) * tileWidth * tileWidth + ly * tileWidth + lx), data, size);
        };
    public:
        ~TiledImage2DTemplate() { };
        
        TiledImage2DTemplate(const void* linearData, uint32_t width, uint32_t height, ColorFormat fmt, Allocator* mem, SpectrumType spType) {
            m_width = width;
            m_height = height;
            m_colorFormat = fmt;
            m_spType = spType;
            
#ifdef Use_Spectral_Representation
            ColorFormat intermFormat;
            switch (m_colorFormat) {
                case ColorFormat::RGB8x3:
                case ColorFormat::RGB_8x4:
                    intermFormat = ColorFormat::uvs16Fx3;
                    break;
                case ColorFormat::RGBA8x4:
                case ColorFormat::RGBA16Fx4:
                    intermFormat = ColorFormat::uvsA16Fx4;
                    break;
                case ColorFormat::Gray8:
                    intermFormat = ColorFormat::Gray8;
                    break;
                default:
                    SLRAssert(false, "Color format is invalid.");
                    break;
            }
            m_stride = sizesOfColorFormats[(uint32_t)intermFormat];
            m_numTileX = (m_width + (tileWidth - 1)) >> log2_tileWidth;
            size_t numTileY = (m_height + (tileWidth - 1)) >> log2_tileWidth;
            size_t tileSize = m_stride * tileWidth * tileWidth;
            m_allocSize = m_numTileX * numTileY * tileSize;
            m_data = (uint8_t*)mem->alloc(m_allocSize, SLR_L1_Cacheline_Size);
            
            for (int i = 0; i < m_height; ++i) {
                for (int j = 0; j < m_width; ++j) {
                    switch (m_colorFormat) {
                        case ColorFormat::RGB8x3: {
                            const RGB8x3 &val = *((RGB8x3*)linearData + m_width * i + j);
                            float RGB[3] = {val.r / 255.0f, val.g / 255.0f, val.b / 255.0f};
                            float uvs[3];
                            Upsampling::sRGB_to_uvs(spType, RGB, uvs);
                            SLRAssert(!std::isnan(uvs[0]) && !std::isinf(uvs[0]) &&
                                      !std::isnan(uvs[1]) && !std::isinf(uvs[1]) &&
                                      !std::isnan(uvs[2]) && !std::isinf(uvs[2]), "Invalid value.");
                            uvs16Fx3 storedVal{(half)uvs[0], (half)uvs[1], (half)uvs[2]};
                            setInternal(j, i, &storedVal, m_stride);
                            break;
                        }
                        case ColorFormat::RGB_8x4: {
                            const RGB_8x4 &val = *((RGB_8x4*)linearData + m_width * i + j);
                            float RGB[3] = {val.r / 255.0f, val.g / 255.0f, val.b / 255.0f};
                            float uvs[3];
                            Upsampling::sRGB_to_uvs(spType, RGB, uvs);
                            SLRAssert(!std::isnan(uvs[0]) && !std::isinf(uvs[0]) &&
                                      !std::isnan(uvs[1]) && !std::isinf(uvs[1]) &&
                                      !std::isnan(uvs[2]) && !std::isinf(uvs[2]), "Invalid value.");
                            uvs16Fx3 storedVal{(half)uvs[0], (half)uvs[1], (half)uvs[2]};
                            setInternal(j, i, &storedVal, m_stride);
                            break;
                        }
                        case ColorFormat::RGBA8x4: {
                            const RGBA8x4 &val = *((RGBA8x4*)linearData + m_width * i + j);
                            float RGB[3] = {val.r / 255.0f, val.g / 255.0f, val.b / 255.0f};
                            float uvs[3];
                            Upsampling::sRGB_to_uvs(spType, RGB, uvs);
                            SLRAssert(!std::isnan(uvs[0]) && !std::isinf(uvs[0]) &&
                                      !std::isnan(uvs[1]) && !std::isinf(uvs[1]) &&
                                      !std::isnan(uvs[2]) && !std::isinf(uvs[2]), "Invalid value.");
                            uvsA16Fx4 storedVal{(half)uvs[0], (half)uvs[1], (half)uvs[2], (half)(val.a / 255.0f)};
                            setInternal(j, i, &storedVal, m_stride);
                            break;
                        }
                        case ColorFormat::RGBA16Fx4: {
                            const RGBA16Fx4 &val = *((RGBA16Fx4*)linearData + m_width * i + j);
                            float RGB[3] = {val.r, val.g, val.b};
                            float uvs[3];
                            Upsampling::sRGB_to_uvs(spType, RGB, uvs);
                            SLRAssert(!std::isnan(uvs[0]) && !std::isinf(uvs[0]) &&
                                      !std::isnan(uvs[1]) && !std::isinf(uvs[1]) &&
                                      !std::isnan(uvs[2]) && !std::isinf(uvs[2]), "Invalid value.");
                            uvsA16Fx4 storedVal{(half)uvs[0], (half)uvs[1], (half)uvs[2], (half)val.a};
                            setInternal(j, i, &storedVal, m_stride);
                            break;
                        }
                        case ColorFormat::Gray8: {
                            const Gray8 &val = *((Gray8*)linearData + m_width * i + j);
                            setInternal(j, i, &val, m_stride);
                            break;
                        }
                        default:
                            break;
                    }
                }
            }
            
            m_colorFormat = intermFormat;
#else
            m_stride = sizesOfColorFormats[(uint32_t)m_colorFormat];
            m_numTileX = (m_width + (tileWidth - 1)) >> log2_tileWidth;
            size_t numTileY = (m_height + (tileWidth - 1)) >> log2_tileWidth;
            size_t tileSize = m_stride * tileWidth * tileWidth;
            m_allocSize = m_numTileX * numTileY * tileSize;
            m_data = (uint8_t*)mem->alloc(m_allocSize, SLR_L1_Cacheline_Size);
            
            for (int i = 0; i < m_height; ++i) {
                for (int j = 0; j < m_width; ++j) {
                    setInternal(j, i, (uint8_t*)linearData + m_stride * (m_width * i + j), m_stride);
                }
            }
#endif
        };
        
        const uint8_t* data() const { return m_data; };
        size_t size() const { return m_allocSize; };
    };    
}

#endif
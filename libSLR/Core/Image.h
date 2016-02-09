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
    
    struct SLR_API RGB8x3 { uint8_t r, g, b; };
    struct SLR_API RGB_8x4 { uint8_t r, g, b, dummy; };
    struct SLR_API RGBA8x4  { uint8_t r, g, b, a; };
    struct SLR_API RGBA16Fx4 { half r, g, b, a; };
    struct SLR_API uvs16Fx3 { half u, v, s; };
    struct SLR_API uvsA16Fx4 { half u, v, s, a; };
    struct SLR_API Gray8 { uint8_t v; };
    
    extern SLR_API const size_t sizesOfColorFormats[(uint32_t)ColorFormat::Num];
    
    class SLR_API Image2D {
    protected:
        uint32_t m_width, m_height;
        ColorFormat m_colorFormat;
        SpectrumType m_spType;
        
        virtual const void* getInternal(uint32_t x, uint32_t y) const = 0;
        virtual void setInternal(uint32_t x, uint32_t y, const void* data, size_t size) = 0;
    public:
        Image2D() { }
        Image2D(uint32_t w, uint32_t h, ColorFormat fmt, SpectrumType spType) : m_width(w), m_height(h), m_colorFormat(fmt), m_spType(spType) { }
        ~Image2D() { }
        
        template <typename ColFmt>
        const ColFmt &get(uint32_t x, uint32_t y) const { return *(const ColFmt*)getInternal(x, y); }
        template <typename ColFmt>
        void set(uint32_t x, uint32_t y, const ColFmt &data) { setInternal(x, y, &data, sizeof(data)); }
        
        void areaAverage(float xLeft, float xRight, float yTop, float yBottom, void* avg) const;
        
        uint32_t width() const { return m_width; }
        uint32_t height() const { return m_height; }
        SpectrumType spectrumType() const { return m_spType; }
        ColorFormat format() const { return m_colorFormat; }
        
        void saveImage(const std::string &filepath, bool gammaCorrection) const;
    };
    
    
    template <uint32_t log2_tileWidth>
    class SLR_API TiledImage2DTemplate : public Image2D {
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
        }
        
        void setInternal(uint32_t x, uint32_t y, const void* data, size_t size) override {
            uint32_t tx = x >> log2_tileWidth;
            uint32_t ty = y >> log2_tileWidth;
            uint32_t lx = x & localMask;
            uint32_t ly = y & localMask;
            std::memcpy(m_data + m_stride * ((ty * m_numTileX + tx) * tileWidth * tileWidth + ly * tileWidth + lx), data, size);
        }
    public:
        ~TiledImage2DTemplate() {
            
        }
        
        TiledImage2DTemplate(uint32_t width, uint32_t height, ColorFormat fmt, Allocator* mem) {
            m_width = width;
            m_height = height;
            m_spType = SpectrumType::Reflectance;
            m_colorFormat = fmt;
            m_stride = sizesOfColorFormats[(uint32_t)m_colorFormat];
            m_numTileX = (m_width + (tileWidth - 1)) >> log2_tileWidth;
            size_t numTileY = (m_height + (tileWidth - 1)) >> log2_tileWidth;
            size_t tileSize = m_stride * tileWidth * tileWidth;
            m_allocSize = m_numTileX * numTileY * tileSize;
            m_data = (uint8_t*)mem->alloc(m_allocSize, SLR_L1_Cacheline_Size);
            
            memset(m_data, 0, m_allocSize);
        }
        
        TiledImage2DTemplate(const void* linearData, uint32_t width, uint32_t height, ColorFormat fmt, Allocator* mem, SpectrumType spType) {
            m_width = width;
            m_height = height;
            m_spType = spType;
            
#ifdef Use_Spectral_Representation
            switch (fmt) {
                case ColorFormat::RGB8x3:
                case ColorFormat::RGB_8x4:
                    m_colorFormat = ColorFormat::uvs16Fx3;
                    break;
                case ColorFormat::RGBA8x4:
                case ColorFormat::RGBA16Fx4:
                    m_colorFormat = ColorFormat::uvsA16Fx4;
                    break;
                case ColorFormat::Gray8:
                    m_colorFormat = ColorFormat::Gray8;
                    break;
                default:
                    SLRAssert(false, "Color format is invalid.");
                    break;
            }
            m_stride = sizesOfColorFormats[(uint32_t)m_colorFormat];
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
                            uvs16Fx3 storedVal{(half)uvs[0], (half)uvs[1], (half)uvs[2]};
                            SLRAssert(!std::isnan((float)storedVal.u) && !std::isinf((float)storedVal.u) &&
                                      !std::isnan((float)storedVal.v) && !std::isinf((float)storedVal.v) &&
                                      !std::isnan((float)storedVal.s) && !std::isinf((float)storedVal.s), "Invalid value.");
                            setInternal(j, i, &storedVal, m_stride);
                            break;
                        }
                        case ColorFormat::RGB_8x4: {
                            const RGB_8x4 &val = *((RGB_8x4*)linearData + m_width * i + j);
                            float RGB[3] = {val.r / 255.0f, val.g / 255.0f, val.b / 255.0f};
                            float uvs[3];
                            Upsampling::sRGB_to_uvs(spType, RGB, uvs);
                            uvs16Fx3 storedVal{(half)uvs[0], (half)uvs[1], (half)uvs[2]};
                            SLRAssert(!std::isnan((float)storedVal.u) && !std::isinf((float)storedVal.u) &&
                                      !std::isnan((float)storedVal.v) && !std::isinf((float)storedVal.v) &&
                                      !std::isnan((float)storedVal.s) && !std::isinf((float)storedVal.s), "Invalid value.");
                            setInternal(j, i, &storedVal, m_stride);
                            break;
                        }
                        case ColorFormat::RGBA8x4: {
                            const RGBA8x4 &val = *((RGBA8x4*)linearData + m_width * i + j);
                            float RGB[3] = {val.r / 255.0f, val.g / 255.0f, val.b / 255.0f};
                            float uvs[3];
                            Upsampling::sRGB_to_uvs(spType, RGB, uvs);
                            uvsA16Fx4 storedVal{(half)uvs[0], (half)uvs[1], (half)uvs[2], (half)(val.a / 255.0f)};
                            SLRAssert(!std::isnan((float)storedVal.u) && !std::isinf((float)storedVal.u) &&
                                      !std::isnan((float)storedVal.v) && !std::isinf((float)storedVal.v) &&
                                      !std::isnan((float)storedVal.s) && !std::isinf((float)storedVal.s) &&
                                      !std::isnan((float)storedVal.a) && !std::isinf((float)storedVal.a) && (float)storedVal.a > 0, "Invalid value.");
                            setInternal(j, i, &storedVal, m_stride);
                            break;
                        }
                        case ColorFormat::RGBA16Fx4: {
                            const RGBA16Fx4 &val = *((RGBA16Fx4*)linearData + m_width * i + j);
                            float RGB[3] = {val.r, val.g, val.b};
                            float uvs[3];
                            Upsampling::sRGB_to_uvs(spType, RGB, uvs);
                            uvsA16Fx4 storedVal{(half)uvs[0], (half)uvs[1], (half)uvs[2], (half)val.a};
                            SLRAssert(!std::isnan((float)storedVal.u) && !std::isinf((float)storedVal.u) &&
                                      !std::isnan((float)storedVal.v) && !std::isinf((float)storedVal.v) &&
                                      !std::isnan((float)storedVal.s) && !std::isinf((float)storedVal.s) &&
                                      !std::isnan((float)storedVal.a) && !std::isinf((float)storedVal.a) && (float)storedVal.a > 0, "Invalid value.");
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
#else
            m_colorFormat = fmt;
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
        }
        
        const uint8_t* data() const { return m_data; }
        size_t size() const { return m_allocSize; }
    };    
}

#endif

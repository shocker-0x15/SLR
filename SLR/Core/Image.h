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
#include "../Helper/image_loader.h"
#include "../BasicTypes/Spectrum.h"

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
    
    static std::map<std::string, Image2DRef> s_database;
};


template <uint32_t log2_tileWidth>
class TiledImage2DTemplate : public Image2D {
    static const size_t tileWidth = 1 << log2_tileWidth;
    static const uint32_t localMask = (1 << log2_tileWidth) - 1;
    size_t m_stride;
    size_t m_numTileX;
    size_t m_allocSize;
    uint8_t* m_data;
    
//    TiledImage2DTemplate(const uint8_t* linearImage, uint32_t w, uint32_t h, ColorFormat fmt, Allocator* mem) : Image2D(w, h, fmt) {
//        m_stride = sizesOfColorFormats[(uint32_t)fmt];
//        m_numTileX = (w + (tileWidth - 1)) & ~(tileWidth - 1);
//        uint32_t numTileY = (h + (tileWidth - 1)) & ~(tileWidth - 1);
//        uint64_t tileSize = m_stride * tileWidth * tileWidth;
//        tileSize = (tileSize + (SLR_L1_Cacheline_Size - 1)) & ~(SLR_L1_Cacheline_Size - 1);
//        m_size = m_numTileX * numTileY * tileSize;
//        m_data = (uint8_t*)mem->alloc(m_size, SLR_L1_Cacheline_Size);
//        
//        for (int i = 0; i < h; ++i) {
//            for (int j = 0; j < w; ++j) {
//                (*this)(j, i) = linearImage[w * i + j];
//            }
//        }
//    };
    
    TiledImage2DTemplate(const std::string &filepath, Allocator* mem, SpectrumType spType, bool gammaCorrection) {
        m_spType = spType;
        
        uint64_t requiredSize;
        bool imgSuccess;
        imgSuccess = getImageInfo(filepath, &m_width, &m_height, &requiredSize, &m_colorFormat);
        SLRAssert(imgSuccess, "Error occured during getting image information.\n%s", filepath.c_str());
        
        void* linearData = malloc(requiredSize);
        imgSuccess = loadImage(filepath, (uint8_t*)linearData, gammaCorrection);
        SLRAssert(imgSuccess, "failed to load the image\n%s", filepath.c_str());
        
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
        free(linearData);
    };
    
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
    ~TiledImage2DTemplate() {
        
    };
    
    static std::shared_ptr<TiledImage2DTemplate> create(const std::string &filepath, Allocator *mem, SpectrumType spType, bool gammaCorrection = false) {
        if (s_database.count(filepath) > 0) {
            return std::static_pointer_cast<TiledImage2DTemplate>(s_database[filepath]);
        }
        else {
            TiledImage2DTemplate* texData = new TiledImage2DTemplate(filepath, mem, spType, gammaCorrection);
            std::shared_ptr<TiledImage2DTemplate> ret = std::shared_ptr<TiledImage2DTemplate>(texData);
            s_database[filepath] = ret;
            return ret;
        }
    };
    
    const uint8_t* data() const { return m_data; };
    size_t size() const { return m_allocSize; };
};

#endif

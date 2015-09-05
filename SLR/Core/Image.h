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

class Image2D {
protected:
    uint32_t m_width, m_height;
    ColorFormat m_colorFormat;
    
    virtual const void* getInternal(uint32_t x, uint32_t y) const = 0;
    virtual void setInternal(uint32_t x, uint32_t y, void* data, size_t size) = 0;
public:
    Image2D() { };
    Image2D(uint32_t w, uint32_t h, ColorFormat fmt) : m_width(w), m_height(h), m_colorFormat(fmt) { };
    
    template <typename ColFmt>
    const ColFmt &get(uint32_t x, uint32_t y) const { return *(const ColFmt*)getInternal(x, y); };
    template <typename ColFmt>
    void set(uint32_t x, uint32_t y, const ColFmt &data) { setInternal(x, y, &data, sizeof(data)); };
    
    void areaAverage(float xLeft, float xRight, float yTop, float yBottom, void* avg) const;
    
    uint32_t width() const { return m_width; };
    uint32_t height() const { return m_height; };
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
    
    TiledImage2DTemplate(const std::string &filepath, Allocator* mem, bool gammaCorrection) {
        uint64_t requiredSize;
        bool imgSuccess;
        imgSuccess = getImageInfo(filepath, &m_width, &m_height, &requiredSize, &m_colorFormat);
        SLRAssert(imgSuccess, "Error occured during getting image information.\n%s", filepath.c_str());
        
        void* linearData = malloc(requiredSize);
        imgSuccess = loadImage(filepath, (uint8_t*)linearData, gammaCorrection);
        SLRAssert(imgSuccess, "failed to load the image\n%s", filepath.c_str());
        
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
        free(linearData);
    };
    
    const void* getInternal(uint32_t x, uint32_t y) const override {
        uint32_t tx = x >> log2_tileWidth;
        uint32_t ty = y >> log2_tileWidth;
        uint32_t lx = x & localMask;
        uint32_t ly = y & localMask;
        return m_data + m_stride * ((ty * m_numTileX + tx) * tileWidth * tileWidth + ly * tileWidth + lx);
    };
    
    void setInternal(uint32_t x, uint32_t y, void* data, size_t size) override {
        uint32_t tx = x >> log2_tileWidth;
        uint32_t ty = y >> log2_tileWidth;
        uint32_t lx = x & localMask;
        uint32_t ly = y & localMask;
        std::memcpy(m_data + m_stride * ((ty * m_numTileX + tx) * tileWidth * tileWidth + ly * tileWidth + lx), data, size);
    };
public:
    static std::shared_ptr<TiledImage2DTemplate> create(const std::string &filepath, Allocator *mem, bool gammaCorrection = false) {
        if (s_database.count(filepath) > 0) {
            return std::static_pointer_cast<TiledImage2DTemplate>(s_database[filepath]);
        }
        else {
            TiledImage2DTemplate* texData = new TiledImage2DTemplate(filepath, mem, gammaCorrection);
            std::shared_ptr<TiledImage2DTemplate> ret = std::shared_ptr<TiledImage2DTemplate>(texData);
            s_database[filepath] = ret;
            return ret;
        }
    };
    
    const uint8_t* data() const { return m_data; };
    size_t size() const { return m_allocSize; };
};

#endif

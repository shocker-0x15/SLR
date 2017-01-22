//
//  ImageSensor.cpp
//
//  Created by 渡部 心 on 2015/08/07.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "ImageSensor.h"
#include "../Helper/bmp_exporter.h"

namespace SLR {
    static const uint32_t s_log2_tileWidth = 3;
    static const uint32_t s_tileWidth = 1 << s_log2_tileWidth;
    static const uint32_t s_localMask = (1 << s_log2_tileWidth) - 1;
    
    ImageSensor::ImageSensor(float sensitivity) :
    m_data(nullptr), m_separatedData(nullptr), m_numSeparated(0), m_sensitivity(sensitivity)
    {}
    
    ImageSensor::ImageSensor(uint32_t width, uint32_t height, float sensitivity) :
    m_data(nullptr), m_separatedData(nullptr), m_numSeparated(0), m_sensitivity(sensitivity) {
        init(width, height);
    }
    
    ImageSensor::~ImageSensor() {
        if (m_data)
            SLR_freealign(m_data);
        if (m_separatedData) {
            for (int i = 0; i < m_numSeparated; ++i)
                SLR_freealign(m_separatedData[i]);
            SLR_freealign(m_separatedData);
        }
    }
    
    void ImageSensor::init(uint32_t width, uint32_t height) {
        m_width = width;
        m_height = height;
        if (m_data)
            SLR_freealign(m_data);
        
        m_numTileX = (width + (s_tileWidth - 1)) >> s_log2_tileWidth;
        m_numTileY = (height + (s_tileWidth - 1)) >> s_log2_tileWidth;
        
        uint64_t tileSize = sizeof(SpectrumStorage) * s_tileWidth * s_tileWidth;
        
        m_allocSize = m_numTileX * m_numTileY * tileSize;
        m_data = (uint8_t*)SLR_memalign(m_allocSize, SLR_L1_Cacheline_Size);
        SLRAssert(m_data, "Failed to allocate frame buffer.");
        
        clear();
    }
    
    void ImageSensor::addSeparatedBuffers(uint32_t numBuffers) {
        uint32_t curIdx = m_numSeparated;
        m_numSeparated = numBuffers;
        m_separatedData = (uint8_t**)SLR_memalign(sizeof(uint8_t*) * numBuffers, SLR_L1_Cacheline_Size);
        for (int i = curIdx; i < m_numSeparated; ++i) {
            m_separatedData[i] = (uint8_t*)SLR_memalign(m_allocSize, SLR_L1_Cacheline_Size);
            SLRAssert(m_separatedData[i], "Failed to allocate a separated buffer.");
        }
        clearSeparatedBuffers();
    }

    uint32_t ImageSensor::tileWidth() const {
        return s_tileWidth;
    }
    
    uint32_t ImageSensor::tileHeight() const {
        return s_tileWidth;
    }
    
    void ImageSensor::clear() {
        for (int i = 0; i < m_allocSize / sizeof(SpectrumStorage); ++i) {
            SpectrumStorage &dst = *((SpectrumStorage*)m_data + i);
            dst = SpectrumStorage(0.0);
        }
    }
    
    void ImageSensor::clearSeparatedBuffers() {
        for (int b = 0; b < m_numSeparated; ++b) {
            for (int i = 0; i < m_allocSize / sizeof(SpectrumStorage); ++i) {
                SpectrumStorage &dst = *((SpectrumStorage*)m_separatedData[b] + i);
                dst = SpectrumStorage(0.0);
            }
        }
    }
    
    DiscretizedSpectrum ImageSensor::pixel(uint32_t x, uint32_t y) const {
        uint32_t tx = x >> s_log2_tileWidth;
        uint32_t ty = y >> s_log2_tileWidth;
        uint32_t lx = x & s_localMask;
        uint32_t ly = y & s_localMask;
        SpectrumStorage &storage = *(SpectrumStorage*)(m_data + sizeof(SpectrumStorage) * ((ty * m_numTileX + tx) * s_tileWidth * s_tileWidth + ly * s_tileWidth + lx));
        return storage.getValue();
    }
    
    SpectrumStorage &ImageSensor::pixel(uint32_t x, uint32_t y) {
        uint32_t tx = x >> s_log2_tileWidth;
        uint32_t ty = y >> s_log2_tileWidth;
        uint32_t lx = x & s_localMask;
        uint32_t ly = y & s_localMask;
        SpectrumStorage &storage = *(SpectrumStorage*)(m_data + sizeof(SpectrumStorage) * ((ty * m_numTileX + tx) * s_tileWidth * s_tileWidth + ly * s_tileWidth + lx));
        return storage;
    }
    
    DiscretizedSpectrum ImageSensor::pixel(uint32_t idx, uint32_t x, uint32_t y) const {
        uint32_t tx = x >> s_log2_tileWidth;
        uint32_t ty = y >> s_log2_tileWidth;
        uint32_t lx = x & s_localMask;
        uint32_t ly = y & s_localMask;
        SpectrumStorage &storage = *(SpectrumStorage*)(m_separatedData[idx] + sizeof(SpectrumStorage) * ((ty * m_numTileX + tx) * s_tileWidth * s_tileWidth + ly * s_tileWidth + lx));
        return storage.getValue();
    }
    
    SpectrumStorage &ImageSensor::pixel(uint32_t idx, uint32_t x, uint32_t y) {
        uint32_t tx = x >> s_log2_tileWidth;
        uint32_t ty = y >> s_log2_tileWidth;
        uint32_t lx = x & s_localMask;
        uint32_t ly = y & s_localMask;
        SpectrumStorage &storage = *(SpectrumStorage*)(m_separatedData[idx] + sizeof(SpectrumStorage) * ((ty * m_numTileX + tx) * s_tileWidth * s_tileWidth + ly * s_tileWidth + lx));
        return storage;
    }
    
    void ImageSensor::add(float px, float py, const WavelengthSamples &wls, const SampledSpectrum &contribution) {
        uint32_t ipx = std::min((uint32_t)px, m_width - 1);
        uint32_t ipy = std::min((uint32_t)py, m_height - 1);
		SLRAssert(contribution.allFinite(), "invalid value: (%u, %u), %s", ipx, ipy, contribution.toString().c_str());
        pixel(ipx, ipy).add(wls, contribution);
    }
    
    void ImageSensor::add(uint32_t idx, float px, float py, const WavelengthSamples &wls, const SampledSpectrum &contribution) {
        uint32_t ipx = std::min((uint32_t)px, m_width - 1);
        uint32_t ipy = std::min((uint32_t)py, m_height - 1);
		SLRAssert(contribution.allFinite(), "invalid value: idx: %u, (%u, %u), %s", idx, ipx, ipy, contribution.toString().c_str());
        pixel(idx, ipx, ipy).add(wls, contribution);
    }
    
    void ImageSensor::saveImage(const std::string &filepath, float scale, float* scaleSeparated) const {
        struct BMP_RGB {
            uint8_t B, G, R;
        };
        
        float* scales = (float*)alloca(sizeof(float) * m_numSeparated);
        float sensitivity = std::isinf(m_sensitivity) ? 1.0f : m_sensitivity;
        for (int i = 0; i < m_numSeparated; ++i)
            scales[i] = (scaleSeparated ? scaleSeparated[i] : scale) * sensitivity;
        scale *= sensitivity;
        
        uint32_t byteWidth = 3 * m_width + m_width % 4;
        uint8_t* bmp = (uint8_t*)malloc(m_height * byteWidth);
        for (int i = 0; i < m_height; ++i) {
            for (int j = 0; j < m_width; ++j) {
                CompensatedSum<DiscretizedSpectrum> pixSum = pixel(j, i) * scale;
                for (int b = 0; b < m_numSeparated; ++b)
                    pixSum += pixel(b, j, i) * scales[b];
                DiscretizedSpectrum pix = pixSum.result;
                if (pix.hasInf())
                    printf("(%u, %u): has an infinite value!\n%s\n", j, i, pix.toString().c_str());
                if (pix.hasNaN())
                    printf("(%u, %u): has NaN!\n%s\n", j, i, pix.toString().c_str());
                if (pix.hasMinus())
                    printf("(%u, %u): has a minus value!\n%s\n", j, i, pix.toString().c_str());
                
                float RGB[3];
                pix.getRGB(RGB);
                RGB[0] = RGB[0] < 0.0f ? 0.0f : RGB[0];
                RGB[1] = RGB[1] < 0.0f ? 0.0f : RGB[1];
                RGB[2] = RGB[2] < 0.0f ? 0.0f : RGB[2];
                
                float Y = 0.222485 * RGB[0] + 0.716905 * RGB[1] + 0.060610 * RGB[2];
                float scaleY = Y != 0 ? (1.0f - std::exp(-Y)) / Y : 0.0f;
                RGB[0] = std::min(scaleY * RGB[0], 1.0f);
                RGB[1] = std::min(scaleY * RGB[1], 1.0f);
                RGB[2] = std::min(scaleY * RGB[2], 1.0f);
                
                uint32_t idx = (m_height - i - 1) * byteWidth + 3 * j;
                BMP_RGB &dst = *(BMP_RGB*)(bmp + idx);
                dst.R = uint8_t(256 * std::min(sRGB_gamma(RGB[0]), 0.999f));
                dst.G = uint8_t(256 * std::min(sRGB_gamma(RGB[1]), 0.999f));
                dst.B = uint8_t(256 * std::min(sRGB_gamma(RGB[2]), 0.999f));
            }
        }
        
        saveBMP(filepath.c_str(), bmp, m_width, m_height);
        free(bmp);
    }
}

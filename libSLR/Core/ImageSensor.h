//
//  ImageSensor.h
//
//  Created by 渡部 心 on 2015/08/07.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_ImageSensor__
#define __SLR_ImageSensor__

#include "../defines.h"
#include "../declarations.h"
#include "../BasicTypes/rgb_types.h"
#include "../BasicTypes/spectrum_types.h"
#include "../BasicTypes/CompensatedSum.h"

namespace SLR {
    class SLR_API ImageSensor {
        uint8_t* m_data;
        uint8_t** m_separatedData;
        uint32_t m_numSeparated;
        uint32_t m_width;
        uint32_t m_height;
        float m_sensitivity;
        
        size_t m_numTileX;
        size_t m_numTileY;
        size_t m_allocSize;
    public:
        ImageSensor(float sensitivity);
        ImageSensor(uint32_t width, uint32_t height, float sensitivity);
        ~ImageSensor();
        
        void init(uint32_t width, uint32_t height);
        void addSeparatedBuffers(uint32_t numBuffers);
        
        void clear();
        void clearSeparatedBuffers();
        
        uint32_t width() const { return m_width; };
        uint32_t height() const { return m_height; };
        uint32_t tileWidth() const;
        uint32_t tileHeight() const;
        uint32_t numTileX() const { return (uint32_t)m_numTileX; };
        uint32_t numTileY() const { return (uint32_t)m_numTileY; };
        
        DiscretizedSpectrum pixel(uint32_t x, uint32_t y) const;
        SpectrumStorage &pixel(uint32_t x, uint32_t y);
        DiscretizedSpectrum pixel(uint32_t idx, uint32_t x, uint32_t y) const;
        SpectrumStorage &pixel(uint32_t idx, uint32_t x, uint32_t y);
        
        void add(float px, float py, const WavelengthSamples &wls, const SampledSpectrum &contribution);
        void add(uint32_t idx, float px, float py, const WavelengthSamples &wls, const SampledSpectrum &contribution);
        
        void saveImage(const std::string &filepath, float scale = 1.0f, float* scaleSeparated = nullptr) const;
    };    
}

#endif /* __SLR_ImageSensor__ */

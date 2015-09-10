//
//  image_texture.cpp
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "image_textures.h"
#include "../Core/Image.h"
#include "../Core/distributions.h"

Spectrum ImageSpectrumTexture::evaluate(const TexCoord2D &tc, const WavelengthSamples &wls) const {
    uint32_t x = std::clamp((uint32_t)std::fmod(m_data->width() * tc.u, m_data->width()), 0u, m_data->width() - 1);
    uint32_t y = std::clamp((uint32_t)std::fmod(m_data->height() * tc.v, m_data->height()), 0u, m_data->height() - 1);
    Spectrum ret;
    switch (m_data->format()) {
        case ColorFormat::RGB8x3: {
            const RGB8x3 &data = m_data->get<RGB8x3>(x, y);
            ret.r = data.r / 255.0f;
            ret.g = data.g / 255.0f;
            ret.b = data.b / 255.0f;
            break;
        }
        case ColorFormat::RGB_8x4: {
            const RGB_8x4 &data = m_data->get<RGB_8x4>(x, y);
            ret.r = data.r / 255.0f;
            ret.g = data.g / 255.0f;
            ret.b = data.b / 255.0f;
            break;
        }
        case ColorFormat::RGBA8x4: {
            const RGBA8x4 &data = m_data->get<RGBA8x4>(x, y);
            ret.r = data.r / 255.0f;
            ret.g = data.g / 255.0f;
            ret.b = data.b / 255.0f;
            break;
        }
        case ColorFormat::RGBA16Fx4: {
            const RGBA16Fx4 &data = m_data->get<RGBA16Fx4>(x, y);
            ret.r = data.r;
            ret.g = data.g;
            ret.b = data.b;
            break;
        }
        case ColorFormat::Gray8: {
            break;
        }
        default:
            break;
    }
    return ret;
}

RegularConstantContinuous2D* ImageSpectrumTexture::createIBLImportanceMap() const {
    uint32_t mapWidth = m_data->width() / 4;
    uint32_t mapHeight = m_data->height() / 4;
    float deltaX = m_data->width() / mapWidth;
    float deltaY = m_data->height() / mapHeight;
    std::function<float(uint32_t, uint32_t)> pickFunc = [this, &deltaX, &deltaY, &mapHeight](uint32_t x, uint32_t y) -> float {
        uint8_t data[16];
        m_data->areaAverage(x * deltaX, (x + 1) * deltaX, y * deltaY, (y + 1) * deltaY, data);
        float luminance;
        switch (m_data->format()) {
            case ColorFormat::RGB8x3: {
                RGB8x3 avg = *(RGB8x3*)data;
                luminance = 0.222485f * avg.r + 0.716905f * avg.g + 0.060610f * avg.b;
                break;
            }
            case ColorFormat::RGB_8x4: {
                RGB_8x4 avg = *(RGB_8x4*)data;
                luminance = 0.222485f * avg.r + 0.716905f * avg.g + 0.060610f * avg.b;
                break;
            }
            case ColorFormat::RGBA8x4: {
                RGBA8x4 avg = *(RGBA8x4*)data;
                luminance = 0.222485f * avg.r + 0.716905f * avg.g + 0.060610f * avg.b;
                break;
            }
            case ColorFormat::RGBA16Fx4: {
                RGBA16Fx4 avg = *(RGBA16Fx4*)data;
                luminance = 0.222485f * avg.r + 0.716905f * avg.g + 0.060610f * avg.b;
                break;
            }
            default:
                return 0.0f;
        }
        return std::sin(M_PI * (y + 0.5f) / mapHeight) * luminance;
    };
    return new RegularConstantContinuous2D(mapWidth, mapHeight, pickFunc);
}

Normal3D ImageNormal3DTexture::evaluate(const TexCoord2D &tc) const {
    uint32_t x = std::clamp((uint32_t)std::fmod(m_data->width() * tc.u, m_data->width()), 0u, m_data->width() - 1);
    uint32_t y = std::clamp((uint32_t)std::fmod(m_data->height() * tc.v, m_data->height()), 0u, m_data->height() - 1);
    Normal3D ret;
    switch (m_data->format()) {
        case ColorFormat::RGB8x3: {
            const RGB8x3 &data = m_data->get<RGB8x3>(x, y);
            ret = normalize(Normal3D(data.r / 255.0f - 0.5f, data.g / 255.0f - 0.5f, data.b / 255.0f - 0.5f));
            break;
        }
        case ColorFormat::RGB_8x4: {
            const RGB_8x4 &data = m_data->get<RGB_8x4>(x, y);
            ret = normalize(Normal3D(data.r / 255.0f - 0.5f, data.g / 255.0f - 0.5f, data.b / 255.0f - 0.5f));
            break;
        }
        case ColorFormat::RGBA8x4: {
            const RGBA8x4 &data = m_data->get<RGBA8x4>(x, y);
            ret = normalize(Normal3D(data.r / 255.0f - 0.5f, data.g / 255.0f - 0.5f, data.b / 255.0f - 0.5f));
            break;
        }
        case ColorFormat::RGBA16Fx4: {
            break;
        }
        case ColorFormat::Gray8: {
            break;
        }
        default:
            break;
    }
    return ret;
}

float ImageFloatTexture::evaluate(const TexCoord2D &tc) const {
    uint32_t x = std::clamp((uint32_t)std::fmod(m_data->width() * tc.u, m_data->width()), 0u, m_data->width() - 1);
    uint32_t y = std::clamp((uint32_t)std::fmod(m_data->height() * tc.v, m_data->height()), 0u, m_data->height() - 1);
    float ret;
    switch (m_data->format()) {
        case ColorFormat::RGB8x3: {
            break;
        }
        case ColorFormat::RGB_8x4: {
            break;
        }
        case ColorFormat::RGBA8x4: {
            break;
        }
        case ColorFormat::RGBA16Fx4: {
            break;
        }
        case ColorFormat::Gray8: {
            const Gray8 &data = m_data->get<Gray8>(x, y);
            ret = data.v / 255.0f;
            break;
        }
        default:
            break;
    }
    return ret;
}

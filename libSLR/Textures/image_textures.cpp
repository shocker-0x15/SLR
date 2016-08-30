//
//  image_texture.cpp
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "image_textures.h"
#include "../Core/Image.h"
#include "../Core/distributions.h"

namespace SLR {
    SampledSpectrum ImageSpectrumTexture::evaluate(const SurfacePoint &surfPt, const WavelengthSamples &wls) const {
        Point3D tc = m_mapping->map(surfPt);
        float u = std::fmod(tc.x, 1.0f);
        float v = std::fmod(tc.y, 1.0f);
        u += u < 0 ? 1.0f : 0.0f;
        v += v < 0 ? 1.0f : 0.0f;
        uint32_t px = std::min((uint32_t)(m_data->width() * u), m_data->width() - 1);
        uint32_t py = std::min((uint32_t)(m_data->height() * v), m_data->height() - 1);
        SampledSpectrum ret;
        switch (m_data->format()) {
#ifdef Use_Spectral_Representation
            case ColorFormat::uvs16Fx3: {
                const uvs16Fx3 &data = m_data->get<uvs16Fx3>(px, py);
                ret = UpsampledContinuousSpectrum(data.u, data.v, data.s / Upsampling::EqualEnergyReflectance).evaluate(wls);
                break;
            }
            case ColorFormat::uvsA16Fx4: {
                const uvsA16Fx4 &data = m_data->get<uvsA16Fx4>(px, py);
                ret = UpsampledContinuousSpectrum(data.u, data.v, data.s / Upsampling::EqualEnergyReflectance).evaluate(wls);
                break;
            }
            case ColorFormat::Gray8: {
                const Gray8 &data = m_data->get<Gray8>(px, py);
                ret = SampledSpectrum(data.v / 255.0f);
                break;
            }
#else
            case ColorFormat::RGB8x3: {
                const RGB8x3 &data = m_data->get<RGB8x3>(px, py);
                ret.r = data.r / 255.0f;
                ret.g = data.g / 255.0f;
                ret.b = data.b / 255.0f;
                break;
            }
            case ColorFormat::RGB_8x4: {
                const RGB_8x4 &data = m_data->get<RGB_8x4>(px, py);
                ret.r = data.r / 255.0f;
                ret.g = data.g / 255.0f;
                ret.b = data.b / 255.0f;
                break;
            }
            case ColorFormat::RGBA8x4: {
                const RGBA8x4 &data = m_data->get<RGBA8x4>(px, py);
                ret.r = data.r / 255.0f;
                ret.g = data.g / 255.0f;
                ret.b = data.b / 255.0f;
                break;
            }
            case ColorFormat::RGBA16Fx4: {
                const RGBA16Fx4 &data = m_data->get<RGBA16Fx4>(px, py);
                ret.r = data.r;
                ret.g = data.g;
                ret.b = data.b;
                break;
            }
            case ColorFormat::Gray8: {
                const Gray8 &data = m_data->get<Gray8>(px, py);
                ret.r = ret.g = ret.b = data.v / 255.0f;
                break;
            }
#endif
            default:
                SLRAssert(false, "Image data format is unknown.");
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
                case ColorFormat::uvs16Fx3: {
                    uvs16Fx3 avg = *(uvs16Fx3*)data;
                    float uvs[3] = {avg.u, avg.v, avg.s};
                    float rgb[3];
                    Upsampling::uvs_to_sRGB(m_data->spectrumType(), uvs, rgb);
                    luminance = 0.222485f * rgb[0] + 0.716905f * rgb[1] + 0.060610f * rgb[2];
                    break;
                }
                case ColorFormat::uvsA16Fx4: {
                    uvsA16Fx4 avg = *(uvsA16Fx4*)data;
                    float uvs[3] = {avg.u, avg.v, avg.s};
                    float rgb[3];
                    Upsampling::uvs_to_sRGB(m_data->spectrumType(), uvs, rgb);
                    luminance = 0.222485f * rgb[0] + 0.716905f * rgb[1] + 0.060610f * rgb[2];
                    break;
                }
                default:
                    return 0.0f;
            }
            SLRAssert(!std::isnan(luminance) && !std::isinf(luminance), "Invalid area average value.");
            return std::sin(M_PI * (y + 0.5f) / mapHeight) * luminance;
        };
        return new RegularConstantContinuous2D(mapWidth, mapHeight, pickFunc);
    }
    
    Normal3D ImageNormal3DTexture::evaluate(const SurfacePoint &surfPt) const {
        Point3D tc = m_mapping->map(surfPt);
        float u = std::fmod(tc.x, 1.0f);
        float v = std::fmod(tc.y, 1.0f);
        u += u < 0 ? 1.0f : 0.0f;
        v += v < 0 ? 1.0f : 0.0f;
        uint32_t px = std::min((uint32_t)(m_data->width() * u), m_data->width() - 1);
        uint32_t py = std::min((uint32_t)(m_data->height() * v), m_data->height() - 1);
        Normal3D ret;
        switch (m_data->format()) {
            case ColorFormat::RGB8x3: {
                const RGB8x3 &data = m_data->get<RGB8x3>(px, py);
                ret = normalize(Normal3D(data.r / 255.0f - 0.5f, data.g / 255.0f - 0.5f, data.b / 255.0f - 0.5f));
                break;
            }
            case ColorFormat::RGB_8x4: {
                const RGB_8x4 &data = m_data->get<RGB_8x4>(px, py);
                ret = normalize(Normal3D(data.r / 255.0f - 0.5f, data.g / 255.0f - 0.5f, data.b / 255.0f - 0.5f));
                break;
            }
            case ColorFormat::RGBA8x4: {
                const RGBA8x4 &data = m_data->get<RGBA8x4>(px, py);
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
    
    float ImageFloatTexture::evaluate(const SurfacePoint &surfPt) const {
        Point3D tc = m_mapping->map(surfPt);
        float u = std::fmod(tc.x, 1.0f);
        float v = std::fmod(tc.y, 1.0f);
        u += u < 0 ? 1.0f : 0.0f;
        v += v < 0 ? 1.0f : 0.0f;
        uint32_t px = std::min((uint32_t)(m_data->width() * u), m_data->width() - 1);
        uint32_t py = std::min((uint32_t)(m_data->height() * v), m_data->height() - 1);
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
                const Gray8 &data = m_data->get<Gray8>(px, py);
                ret = data.v / 255.0f;
                break;
            }
            case ColorFormat::uvsA16Fx4: {
                const uvsA16Fx4 &data = m_data->get<uvsA16Fx4>(px, py);
                ret = data.a;
                break;
            }
            default:
                break;
        }
        return ret;
    }    
}

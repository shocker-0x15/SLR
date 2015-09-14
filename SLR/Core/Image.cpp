//
//  Image.cpp
//
//  Created by 渡部 心 on 2015/07/11.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "Image.h"
#include "../BasicTypes/CompensatedSum.h"

std::map<std::string, Image2DRef> Image2D::s_database;

void Image2D::areaAverage(float xLeft, float xRight, float yTop, float yBottom, void *avg) const {
    uint32_t xLeftPix = (uint32_t)xLeft;
    uint32_t xRightPix = (uint32_t)ceilf(xRight) - 1;
    uint32_t yTopPix = (uint32_t)yTop;
    uint32_t yBottomPix = (uint32_t)ceilf(yBottom) - 1;
    float area = (yBottom - yTop) * (xRight - xLeft);
    
    // UL, UR, LL, LR
    float coeffsCorners[] = {
        (xLeftPix + 1 - xLeft) * (yTopPix + 1 - yTop),
        (xRight - xRightPix) * (yTopPix + 1 - yTop),
        (xLeftPix + 1 - xLeft) * (yBottom - yBottomPix),
        (xRight - xRightPix) * (yBottom - yBottomPix)
    };
    // Top, Left, Right, Bottom
    float coeffsEdges[] = {
        yTopPix + 1 - yTop,
        xLeftPix + 1 - xLeft,
        xRight - xRightPix,
        yBottom - yBottomPix
    };
    
    switch (m_colorFormat) {
        case ColorFormat::RGB8x3: {
            SLRAssert(false, "Not implemented.");
            break;
        }
        case ColorFormat::RGB_8x4: {
            SLRAssert(false, "Not implemented.");
            break;
        }
        case ColorFormat::RGBA8x4: {
            SLRAssert(false, "Not implemented.");
            break;
        }
        case ColorFormat::RGBA16Fx4: {
            FloatSum sumR(0), sumG(0), sumB(0), sumA(0);
            RGBA16Fx4 pix;
            
            uint32_t corners[] = {xLeftPix, yTopPix, xRightPix, yTopPix, xLeftPix, yBottomPix, xRightPix, yBottomPix};
            for (int i = 0; i < 4; ++i) {
                pix = get<RGBA16Fx4>(corners[2 * i + 0], corners[2 * i + 1]);
                sumR += coeffsCorners[i] * float(pix.r);
                sumG += coeffsCorners[i] * float(pix.g);
                sumB += coeffsCorners[i] * float(pix.b);
                sumA += coeffsCorners[i] * float(pix.a);
            }
            
            for (uint32_t x = xLeftPix + 1; x < xRightPix; ++x) {
                pix = get<RGBA16Fx4>(x, yTopPix);
                sumR += coeffsEdges[0] * float(pix.r);
                sumG += coeffsEdges[0] * float(pix.g);
                sumB += coeffsEdges[0] * float(pix.b);
                sumA += coeffsEdges[0] * float(pix.a);
                
                pix = get<RGBA16Fx4>(x, yBottomPix);
                sumR += coeffsEdges[3] * float(pix.r);
                sumG += coeffsEdges[3] * float(pix.g);
                sumB += coeffsEdges[3] * float(pix.b);
                sumA += coeffsEdges[3] * float(pix.a);
            }
            for (uint32_t y = yTopPix + 1; y < yBottomPix; ++y) {
                pix = get<RGBA16Fx4>(xLeftPix, y);
                sumR += coeffsEdges[1] * float(pix.r);
                sumG += coeffsEdges[1] * float(pix.g);
                sumB += coeffsEdges[1] * float(pix.b);
                sumA += coeffsEdges[1] * float(pix.a);
                
                pix = get<RGBA16Fx4>(xRightPix, y);
                sumR += coeffsEdges[2] * float(pix.r);
                sumG += coeffsEdges[2] * float(pix.g);
                sumB += coeffsEdges[2] * float(pix.b);
                sumA += coeffsEdges[2] * float(pix.a);
            }
            
            for (uint32_t y = yTopPix + 1; y < yBottomPix; ++y) {
                for (uint32_t x = xLeftPix + 1; x < xRightPix; ++x) {
                    pix = get<RGBA16Fx4>(x, y);
                    sumR += float(pix.r);
                    sumG += float(pix.g);
                    sumB += float(pix.b);
                    sumA += float(pix.a);
                }
            }
            
            RGBA16Fx4 ret{half(sumR / area), half(sumG / area), half(sumB / area), half(sumA / area)};
            memcpy(avg, &ret, sizeof(ret));
            break;
        }
        case ColorFormat::uvs16Fx3: {
            FloatSum sumR(0), sumG(0), sumB(0);
            uvs16Fx3 pix;
            
            float uvs[3];
            float rgb[3];
            
            uint32_t corners[] = {xLeftPix, yTopPix, xRightPix, yTopPix, xLeftPix, yBottomPix, xRightPix, yBottomPix};
            for (int i = 0; i < 4; ++i) {
                pix = get<uvs16Fx3>(corners[2 * i + 0], corners[2 * i + 1]);
                uvs[0] = pix.u;
                uvs[1] = pix.v;
                uvs[2] = pix.s;
                Upsampling::uvs_to_sRGB(uvs, rgb);
                sumR += coeffsCorners[i] * rgb[0];
                sumG += coeffsCorners[i] * rgb[1];
                sumB += coeffsCorners[i] * rgb[2];
            }
            
            for (uint32_t x = xLeftPix + 1; x < xRightPix; ++x) {
                pix = get<uvs16Fx3>(x, yTopPix);
                uvs[0] = pix.u;
                uvs[1] = pix.v;
                uvs[2] = pix.s;
                Upsampling::uvs_to_sRGB(uvs, rgb);
                sumR += coeffsEdges[0] * rgb[0];
                sumG += coeffsEdges[0] * rgb[1];
                sumB += coeffsEdges[0] * rgb[2];
                
                pix = get<uvs16Fx3>(x, yBottomPix);
                uvs[0] = pix.u;
                uvs[1] = pix.v;
                uvs[2] = pix.s;
                Upsampling::uvs_to_sRGB(uvs, rgb);
                sumR += coeffsEdges[3] * rgb[0];
                sumG += coeffsEdges[3] * rgb[1];
                sumB += coeffsEdges[3] * rgb[2];
            }
            for (uint32_t y = yTopPix + 1; y < yBottomPix; ++y) {
                pix = get<uvs16Fx3>(xLeftPix, y);
                uvs[0] = pix.u;
                uvs[1] = pix.v;
                uvs[2] = pix.s;
                Upsampling::uvs_to_sRGB(uvs, rgb);
                sumR += coeffsEdges[1] * rgb[0];
                sumG += coeffsEdges[1] * rgb[1];
                sumB += coeffsEdges[1] * rgb[2];
                
                pix = get<uvs16Fx3>(xRightPix, y);
                uvs[0] = pix.u;
                uvs[1] = pix.v;
                uvs[2] = pix.s;
                Upsampling::uvs_to_sRGB(uvs, rgb);
                sumR += coeffsEdges[2] * rgb[0];
                sumG += coeffsEdges[2] * rgb[1];
                sumB += coeffsEdges[2] * rgb[2];
            }
            
            for (uint32_t y = yTopPix + 1; y < yBottomPix; ++y) {
                for (uint32_t x = xLeftPix + 1; x < xRightPix; ++x) {
                    pix = get<uvs16Fx3>(x, y);
                    uvs[0] = pix.u;
                    uvs[1] = pix.v;
                    uvs[2] = pix.s;
                    Upsampling::uvs_to_sRGB(uvs, rgb);
                    sumR += rgb[0];
                    sumG += rgb[1];
                    sumB += rgb[2];
                }
            }
            rgb[0] = sumR / area;
            rgb[1] = sumG / area;
            rgb[2] = sumB / area;
            Upsampling::sRGB_to_uvs(rgb, uvs);
            
            uvs16Fx3 ret{half(uvs[0]), half(uvs[1]), half(uvs[2])};
            memcpy(avg, &ret, sizeof(ret));
            break;
        }
        case ColorFormat::uvsA16Fx4: {
            FloatSum sumR(0), sumG(0), sumB(0), sumA(0);
            uvsA16Fx4 pix;
            
            float uvs[3];
            float rgb[3];
            
            uint32_t corners[] = {xLeftPix, yTopPix, xRightPix, yTopPix, xLeftPix, yBottomPix, xRightPix, yBottomPix};
            for (int i = 0; i < 4; ++i) {
                pix = get<uvsA16Fx4>(corners[2 * i + 0], corners[2 * i + 1]);
                uvs[0] = pix.u;
                uvs[1] = pix.v;
                uvs[2] = pix.s;
                Upsampling::uvs_to_sRGB(uvs, rgb);
                sumR += coeffsCorners[i] * rgb[0];
                sumG += coeffsCorners[i] * rgb[1];
                sumB += coeffsCorners[i] * rgb[2];
                sumA += coeffsCorners[i] * pix.a;
            }
            
            for (uint32_t x = xLeftPix + 1; x < xRightPix; ++x) {
                pix = get<uvsA16Fx4>(x, yTopPix);
                uvs[0] = pix.u;
                uvs[1] = pix.v;
                uvs[2] = pix.s;
                Upsampling::uvs_to_sRGB(uvs, rgb);
                sumR += coeffsEdges[0] * rgb[0];
                sumG += coeffsEdges[0] * rgb[1];
                sumB += coeffsEdges[0] * rgb[2];
                sumA += coeffsEdges[0] * pix.a;
                
                pix = get<uvsA16Fx4>(x, yBottomPix);
                uvs[0] = pix.u;
                uvs[1] = pix.v;
                uvs[2] = pix.s;
                Upsampling::uvs_to_sRGB(uvs, rgb);
                sumR += coeffsEdges[3] * rgb[0];
                sumG += coeffsEdges[3] * rgb[1];
                sumB += coeffsEdges[3] * rgb[2];
                sumA += coeffsEdges[3] * pix.a;
            }
            for (uint32_t y = yTopPix + 1; y < yBottomPix; ++y) {
                pix = get<uvsA16Fx4>(xLeftPix, y);
                uvs[0] = pix.u;
                uvs[1] = pix.v;
                uvs[2] = pix.s;
                Upsampling::uvs_to_sRGB(uvs, rgb);
                sumR += coeffsEdges[1] * rgb[0];
                sumG += coeffsEdges[1] * rgb[1];
                sumB += coeffsEdges[1] * rgb[2];
                sumA += coeffsEdges[1] * pix.a;
                
                pix = get<uvsA16Fx4>(xRightPix, y);
                uvs[0] = pix.u;
                uvs[1] = pix.v;
                uvs[2] = pix.s;
                Upsampling::uvs_to_sRGB(uvs, rgb);
                sumR += coeffsEdges[2] * rgb[0];
                sumG += coeffsEdges[2] * rgb[1];
                sumB += coeffsEdges[2] * rgb[2];
                sumA += coeffsEdges[2] * pix.a;
            }
            
            for (uint32_t y = yTopPix + 1; y < yBottomPix; ++y) {
                for (uint32_t x = xLeftPix + 1; x < xRightPix; ++x) {
                    pix = get<uvsA16Fx4>(x, y);
                    uvs[0] = pix.u;
                    uvs[1] = pix.v;
                    uvs[2] = pix.s;
                    Upsampling::uvs_to_sRGB(uvs, rgb);
                    sumR += rgb[0];
                    sumG += rgb[1];
                    sumB += rgb[2];
                    sumA += pix.a;
                }
            }
            rgb[0] = sumR / area;
            rgb[1] = sumG / area;
            rgb[2] = sumB / area;
            Upsampling::sRGB_to_uvs(rgb, uvs);
            
            uvsA16Fx4 ret{half(uvs[0]), half(uvs[1]), half(uvs[2]), half(sumA / area)};
            memcpy(avg, &ret, sizeof(ret));
            break;
        }
        case ColorFormat::Gray8: {
            SLRAssert(false, "Not implemented.");
            break;
        }
        default:
            break;
    }
}

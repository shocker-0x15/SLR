//
//  Image.cpp
//
//  Created by 渡部 心 on 2015/07/11.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "Image.h"
#include "CompensatedSum.h"

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
            break;
        }
        case ColorFormat::RGB_8x4: {
            break;
        }
        case ColorFormat::RGBA8x4: {
            break;
        }
        case ColorFormat::RGBA16Fx4: {
            FloatSum sumR, sumG, sumB, sumA;
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
        }
        case ColorFormat::Gray8: {
            break;
        }
        default:
            break;
    }
}

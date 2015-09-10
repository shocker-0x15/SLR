//
//  checker_board_textures.h
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__checker_board_textures__
#define __SLR__checker_board_textures__

#include "../defines.h"
#include "../references.h"
#include "../Core/textures.h"

class CheckerBoardSpectrumTexture : public SpectrumTexture {
    Spectrum m_values[2];
public:
    CheckerBoardSpectrumTexture(const Spectrum &v0, const Spectrum &v1) : m_values{v0, v1} { };
    
    Spectrum evaluate(const TexCoord2D &tc, const WavelengthSamples &wls) const override { return m_values[((int)(tc.u * 2) + (int)(tc.v * 2)) % 2]; };
    RegularConstantContinuous2D* createIBLImportanceMap() const override;
};

class CheckerBoardNormal3DTexture : public Normal3DTexture {
    float m_stepWidth;
    bool m_reverse;
public:
    CheckerBoardNormal3DTexture(float stepWidth, bool reverse) :
    m_stepWidth(stepWidth), m_reverse(reverse) {
        SLRAssert(stepWidth > 0 && stepWidth <= 1.0f, "stepWidth must be in the range (0, 1].");
    };
    
    Normal3D evaluate(const TexCoord2D &tc) const override;
};

class CheckerBoardFloatTexture : public FloatTexture {
    float m_values[2];
public:
    CheckerBoardFloatTexture(const float &v0, const float &v1) : m_values{v0, v1} { };
    
    float evaluate(const TexCoord2D &tc) const override { return m_values[((int)(tc.u * 2) + (int)(tc.v * 2)) % 2]; };
};

#endif /* defined(__SLR__checker_board_textures__) */

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

namespace SLR {
    class SLR_API CheckerBoardSpectrumTexture : public SpectrumTexture {
        const InputSpectrum* m_values[2];
    public:
        CheckerBoardSpectrumTexture(const InputSpectrum* v0, const InputSpectrum* v1) : m_values{v0, v1} { };
        
        SampledSpectrum evaluate(const TexCoord2D &tc, const WavelengthSamples &wls) const override { return m_values[((int)(tc.u * 2) + (int)(tc.v * 2)) % 2]->evaluate(wls); };
        RegularConstantContinuous2D* createIBLImportanceMap() const override;
    };
    
    class SLR_API CheckerBoardNormal3DTexture : public Normal3DTexture {
        float m_stepWidth;
        bool m_reverse;
    public:
        CheckerBoardNormal3DTexture(float stepWidth, bool reverse) :
        m_stepWidth(stepWidth), m_reverse(reverse) {
            SLRAssert(stepWidth > 0 && stepWidth <= 1.0f, "stepWidth must be in the range (0, 1].");
        };
        
        Normal3D evaluate(const TexCoord2D &tc) const override;
    };
    
    class SLR_API CheckerBoardFloatTexture : public FloatTexture {
        float m_values[2];
    public:
        CheckerBoardFloatTexture(float v0, float v1) : m_values{v0, v1} { };
        
        float evaluate(const TexCoord2D &tc) const override { return m_values[((int)(tc.u * 2) + (int)(tc.v * 2)) % 2]; };
    };
}

#endif /* defined(__SLR__checker_board_textures__) */

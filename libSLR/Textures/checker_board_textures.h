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
        const Texture2DMapping* m_mapping;
        const InputSpectrum* m_values[2];
    public:
        CheckerBoardSpectrumTexture(const Texture2DMapping* mapping, const InputSpectrum* v0, const InputSpectrum* v1) : m_mapping(mapping), m_values{v0, v1} { }
        
        SampledSpectrum evaluate(const SurfacePoint &surfPt, const WavelengthSamples &wls) const override {
            Point3D tc = m_mapping->map(surfPt);
            return m_values[((int)(tc.x * 2) + (int)(tc.y * 2)) % 2]->evaluate(wls);
        }
        RegularConstantContinuous2D* createIBLImportanceMap() const override;
    };
    
    class SLR_API CheckerBoardNormal3DTexture : public Normal3DTexture {
        const Texture2DMapping* m_mapping;
        float m_stepWidth;
        bool m_reverse;
    public:
        CheckerBoardNormal3DTexture(const Texture2DMapping* mapping, float stepWidth, bool reverse) :
        m_mapping(mapping), m_stepWidth(stepWidth), m_reverse(reverse) {
            SLRAssert(stepWidth > 0 && stepWidth <= 1.0f, "stepWidth must be in the range (0, 1].");
        }
        
        Normal3D evaluate(const SurfacePoint &surfPt) const override;
    };
    
    class SLR_API CheckerBoardFloatTexture : public FloatTexture {
        const Texture2DMapping* m_mapping;
        float m_values[2];
    public:
        CheckerBoardFloatTexture(const Texture2DMapping* mapping, float v0, float v1) : m_mapping(mapping), m_values{v0, v1} { }
        
        float evaluate(const SurfacePoint &surfPt) const override {
            Point3D tc = m_mapping->map(surfPt);
            return m_values[((int)(tc.x * 2) + (int)(tc.y * 2)) % 2];
        }
    };
}

#endif /* defined(__SLR__checker_board_textures__) */

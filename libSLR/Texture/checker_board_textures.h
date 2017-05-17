//
//  checker_board_textures.h
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__checker_board_textures__
#define __SLR__checker_board_textures__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/textures.h"

namespace SLR {
    class SLR_API CheckerBoardSpectrumTexture : public SpectrumTexture {
        const Texture2DMapping* m_mapping;
        const AssetSpectrum* m_values[2];
        float m_luminances[2];
    public:
        CheckerBoardSpectrumTexture(const Texture2DMapping* mapping, const AssetSpectrum* v0, const AssetSpectrum* v1) : m_mapping(mapping), m_values{v0, v1} { }
        
        SampledSpectrum evaluate(const Point3D &p, const WavelengthSamples &wls) const {
            return m_values[((int)(p.x * 2) + (int)(p.y * 2)) % 2]->evaluate(wls);
        }
        SampledSpectrum evaluate(const SurfacePoint &surfPt, const WavelengthSamples &wls) const override {
            return evaluate(m_mapping->map(surfPt), wls);
        }
        SampledSpectrum evaluate(const MediumPoint &medPt, const WavelengthSamples &wls) const override {
            return evaluate(m_mapping->map(medPt), wls);
        }
        float evaluateLuminance(const Point3D &p) const {
            return m_luminances[((int)(p.x * 2) + (int)(p.y * 2)) % 2];
        }
        float evaluateLuminance(const SurfacePoint &surfPt) const override {
            return evaluateLuminance(m_mapping->map(surfPt));
        }
        float evaluateLuminance(const MediumPoint &medPt) const override {
            return evaluateLuminance(m_mapping->map(medPt));
        }
        const ContinuousDistribution2D* createIBLImportanceMap() const override;
        
        void generateLuminanceChannel();
    };
    
    
    
    class SLR_API CheckerBoardNormalTexture : public NormalTexture {
        const Texture2DMapping* m_mapping;
        float m_stepWidth;
        bool m_reverse;
    public:
        CheckerBoardNormalTexture(const Texture2DMapping* mapping, float stepWidth, bool reverse) :
        m_mapping(mapping), m_stepWidth(stepWidth), m_reverse(reverse) {
            SLRAssert(stepWidth > 0 && stepWidth <= 1.0f, "stepWidth must be in the range (0, 1].");
        }
        
        Normal3D evaluate(const Point3D &p) const;
        Normal3D evaluate(const SurfacePoint &surfPt) const override {
            return evaluate(m_mapping->map(surfPt));
        }
        Normal3D evaluate(const MediumPoint &medPt) const override {
            return evaluate(m_mapping->map(medPt));
        }
    };
    
    
    
    class SLR_API CheckerBoardFloatTexture : public FloatTexture {
        const Texture2DMapping* m_mapping;
        float m_values[2];
    public:
        CheckerBoardFloatTexture(const Texture2DMapping* mapping, float v0, float v1) : m_mapping(mapping), m_values{v0, v1} { }
        
        float evaluate(const Point3D &p) const {
            return m_values[((int)(p.x * 2) + (int)(p.y * 2)) % 2];
        }
        float evaluate(const SurfacePoint &surfPt) const override {
            return evaluate(m_mapping->map(surfPt));
        }
        float evaluate(const MediumPoint &medPt) const override {
            return evaluate(m_mapping->map(medPt));
        }
    };
}

#endif /* defined(__SLR__checker_board_textures__) */

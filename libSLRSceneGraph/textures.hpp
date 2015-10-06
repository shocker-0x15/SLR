//
//  textures.hpp
//
//  Created by 渡部 心 on 2015/10/06.
//  Copyright © 2015年 渡部 心. All rights reserved.
//

#ifndef textures_hpp
#define textures_hpp

#include <libSLR/defines.h>
#include "references.h"

namespace SLRSceneGraph {
    class SpectrumTexture {
    protected:
        SLR::SpectrumTexture* m_rawData;
    public:
        virtual ~SpectrumTexture();
        const SLR::SpectrumTexture* getRaw() const {
            return m_rawData;
        };
    };
    
    class Normal3DTexture {
    protected:
        SLR::Normal3DTexture* m_rawData;
    public:
        virtual ~Normal3DTexture();
        const SLR::Normal3DTexture* getRaw() const {
            return m_rawData;
        };
    };
    
    class FloatTexture {
    protected:
        SLR::FloatTexture* m_rawData;
    public:
        virtual ~FloatTexture();
        const SLR::FloatTexture* getRaw() const {
            return m_rawData;
        };
    };
    
    
    class ConstantSpectrumTexture : public SpectrumTexture {
        InputSpectrumRef m_value;
    public:
        ConstantSpectrumTexture(const InputSpectrumRef &value);
    };
    
    class ConstantFloatTexture : public FloatTexture {
        float m_value;
    public:
        ConstantFloatTexture(float value);
    };
    
    
    class ImageSpectrumTexture : public SpectrumTexture {
        TiledImage2DRef m_data;
    public:
        ImageSpectrumTexture(const TiledImage2DRef &image);
    };
    
    class ImageNormal3DTexture : public Normal3DTexture {
        TiledImage2DRef m_data;
    public:
        ImageNormal3DTexture(const TiledImage2DRef &image);
    };
    
    class ImageFloatTexture : public FloatTexture {
        TiledImage2DRef m_data;
    public:
        ImageFloatTexture(const TiledImage2DRef &image);
    };
    
    
    class CheckerBoardSpectrumTexture : public SpectrumTexture {
        InputSpectrumRef m_values[2];
    public:
        CheckerBoardSpectrumTexture(const InputSpectrumRef &v0, const InputSpectrumRef &v1);
    };
    
    class CheckerBoardNormal3DTexture : public Normal3DTexture {
        float m_stepWidth;
        bool m_reverse;
    public:
        CheckerBoardNormal3DTexture(float stepWidth, bool reverse);
    };
    
    class CheckerBoardFloatTexture : public FloatTexture {
        float m_values[2];
    public:
        CheckerBoardFloatTexture(float v0, float v1);
    };
}

#endif /* textures_hpp */

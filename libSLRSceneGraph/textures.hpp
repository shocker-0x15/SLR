//
//  textures.hpp
//
//  Created by 渡部 心 on 2015/10/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef textures_hpp
#define textures_hpp

#include <libSLR/defines.h>
#include "references.h"

namespace SLRSceneGraph {
    class SLR_SCENEGRAPH_API Texture2DMapping {
    protected:
        SLR::Texture2DMapping* m_rawData;
    public:
        Texture2DMapping();
        virtual ~Texture2DMapping();
        const SLR::Texture2DMapping* getRaw() const {
            return m_rawData;
        };
        
        static Texture2DMappingRef &sharedInstanceRef() {
            static Texture2DMappingRef s_sharedInstanceRef = createShared<Texture2DMapping>();
            return s_sharedInstanceRef;
        };
    };
    
    class SLR_SCENEGRAPH_API Texture3DMapping {
    protected:
        SLR::Texture3DMapping* m_rawData;
    public:
        Texture3DMapping();
        virtual ~Texture3DMapping();
        const SLR::Texture3DMapping* getRaw() const {
            return m_rawData;
        };
        
        static Texture3DMappingRef &sharedInstanceRef() {
            static Texture3DMappingRef s_sharedInstanceRef = createShared<Texture3DMapping>();
            return s_sharedInstanceRef;
        };
    };
    
    class SLR_SCENEGRAPH_API OffsetAndScale2DMapping : public Texture2DMapping {
        float m_offsetX, m_offsetY;
        float m_scaleX, m_scaleY;
    public:
        OffsetAndScale2DMapping(float ox, float oy, float sx, float sy);
    };
    
    class SLR_SCENEGRAPH_API WorldPosition3DMapping : public Texture3DMapping {
    public:
        WorldPosition3DMapping();
        
        static Texture3DMappingRef &sharedInstanceRef() {
            static Texture3DMappingRef s_sharedInstanceRef = createShared<WorldPosition3DMapping>();
            return s_sharedInstanceRef;
        };
    };
    
    
    
    class SLR_SCENEGRAPH_API SpectrumTexture {
    protected:
        SLR::SpectrumTexture* m_rawData;
    public:
        virtual ~SpectrumTexture();
        const SLR::SpectrumTexture* getRaw() const {
            return m_rawData;
        };
    };
    
    class SLR_SCENEGRAPH_API Normal3DTexture {
    protected:
        SLR::Normal3DTexture* m_rawData;
    public:
        virtual ~Normal3DTexture();
        const SLR::Normal3DTexture* getRaw() const {
            return m_rawData;
        };
    };
    
    class SLR_SCENEGRAPH_API FloatTexture {
    protected:
        SLR::FloatTexture* m_rawData;
    public:
        virtual ~FloatTexture();
        const SLR::FloatTexture* getRaw() const {
            return m_rawData;
        };
    };
    
    
    class SLR_SCENEGRAPH_API ConstantSpectrumTexture : public SpectrumTexture {
        InputSpectrumRef m_value;
    public:
        ConstantSpectrumTexture(const InputSpectrumRef &value);
    };
    
    class SLR_SCENEGRAPH_API ConstantFloatTexture : public FloatTexture {
        float m_value;
    public:
        ConstantFloatTexture(float value);
    };
    
    
    class SLR_SCENEGRAPH_API ImageSpectrumTexture : public SpectrumTexture {
        Texture2DMappingRef m_mapping;
        TiledImage2DRef m_data;
    public:
        ImageSpectrumTexture(const Texture2DMappingRef &mapping, const TiledImage2DRef &image);
    };
    
    class SLR_SCENEGRAPH_API ImageNormal3DTexture : public Normal3DTexture {
        Texture2DMappingRef m_mapping;
        TiledImage2DRef m_data;
    public:
        ImageNormal3DTexture(const Texture2DMappingRef &mapping, const TiledImage2DRef &image);
    };
    
    class SLR_SCENEGRAPH_API ImageFloatTexture : public FloatTexture {
        Texture2DMappingRef m_mapping;
        TiledImage2DRef m_data;
    public:
        ImageFloatTexture(const Texture2DMappingRef &mapping, const TiledImage2DRef &image);
    };
    
    
    class SLR_SCENEGRAPH_API CheckerBoardSpectrumTexture : public SpectrumTexture {
        Texture2DMappingRef m_mapping;
        InputSpectrumRef m_values[2];
    public:
        CheckerBoardSpectrumTexture(const Texture2DMappingRef &mapping, const InputSpectrumRef &v0, const InputSpectrumRef &v1);
    };
    
    class SLR_SCENEGRAPH_API CheckerBoardNormal3DTexture : public Normal3DTexture {
        Texture2DMappingRef m_mapping;
        float m_stepWidth;
        bool m_reverse;
    public:
        CheckerBoardNormal3DTexture(const Texture2DMappingRef &mapping, float stepWidth, bool reverse);
    };
    
    class SLR_SCENEGRAPH_API CheckerBoardFloatTexture : public FloatTexture {
        Texture2DMappingRef m_mapping;
        float m_values[2];
    public:
        CheckerBoardFloatTexture(const Texture2DMappingRef &mapping, float v0, float v1);
    };
    
    
    class SLR_SCENEGRAPH_API VoronoiSpectrumTexture : public SpectrumTexture {
        Texture3DMappingRef m_mapping;
        float m_scale;
        float m_brightness;
    public:
        VoronoiSpectrumTexture(const Texture3DMappingRef &mapping, float scale, float brightness);
    };
    
    class SLR_SCENEGRAPH_API VoronoiNormal3DTexture : public Normal3DTexture {
        Texture3DMappingRef m_mapping;
        float m_scale;
        float m_thetaMax;
    public:
        VoronoiNormal3DTexture(const Texture3DMappingRef &mapping, float scale, float thetaMax);
    };
    
    class SLR_SCENEGRAPH_API VoronoiFloatTexture : public FloatTexture {
        Texture3DMappingRef m_mapping;
        float m_scale;
        float m_valueScale;
        bool m_flat;
    public:
        VoronoiFloatTexture(const Texture3DMappingRef &mapping, float scale, float valueScale, bool flat);
    };
}

#endif /* textures_hpp */

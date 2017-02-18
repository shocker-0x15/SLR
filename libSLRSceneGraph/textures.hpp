//
//  textures.hpp
//
//  Created by 渡部 心 on 2015/10/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLRSceneGraph_textures__
#define __SLRSceneGraph_textures__

#include <libSLR/defines.h>
#include "declarations.h"

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
    
    class SLR_SCENEGRAPH_API NormalTexture {
    protected:
        SLR::NormalTexture* m_rawData;
    public:
        virtual ~NormalTexture();
        const SLR::NormalTexture* getRaw() const {
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
        AssetSpectrumRef m_value;
    public:
        ConstantSpectrumTexture(const AssetSpectrumRef &value);
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
    
    class SLR_SCENEGRAPH_API ImageNormalTexture : public NormalTexture {
        Texture2DMappingRef m_mapping;
        TiledImage2DRef m_data;
    public:
        ImageNormalTexture(const Texture2DMappingRef &mapping, const TiledImage2DRef &image);
    };
    
    class SLR_SCENEGRAPH_API ImageFloatTexture : public FloatTexture {
        Texture2DMappingRef m_mapping;
        TiledImage2DRef m_data;
    public:
        ImageFloatTexture(const Texture2DMappingRef &mapping, const TiledImage2DRef &image);
    };
    
    
    class SLR_SCENEGRAPH_API CheckerBoardSpectrumTexture : public SpectrumTexture {
        Texture2DMappingRef m_mapping;
        AssetSpectrumRef m_values[2];
    public:
        CheckerBoardSpectrumTexture(const Texture2DMappingRef &mapping, const AssetSpectrumRef &v0, const AssetSpectrumRef &v1);
    };
    
    class SLR_SCENEGRAPH_API CheckerBoardNormalTexture : public NormalTexture {
        Texture2DMappingRef m_mapping;
        float m_stepWidth;
        bool m_reverse;
    public:
        CheckerBoardNormalTexture(const Texture2DMappingRef &mapping, float stepWidth, bool reverse);
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
    
    class SLR_SCENEGRAPH_API VoronoiNormalTexture : public NormalTexture {
        Texture3DMappingRef m_mapping;
        float m_scale;
        float m_thetaMax;
    public:
        VoronoiNormalTexture(const Texture3DMappingRef &mapping, float scale, float thetaMax);
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

#endif /* __SLRSceneGraph_textures__ */

//
//  textures.h
//
//  Created by 渡部 心 on 2015/04/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__textures__
#define __SLR__textures__

#include "../defines.h"
#include "../references.h"
#include "geometry.h"

namespace SLR {
    class SLR_API Texture2DMapping {
    public:
        virtual ~Texture2DMapping() {}
        virtual Point3D map(const SurfacePoint &surfPt) const {
            TexCoord2D texCoord = surfPt.getTextureCoordinate();
            return Point3D(texCoord.u, texCoord.v, 0.0f);
        }
    };
    
    class SLR_API Texture3DMapping {
    public:
        virtual ~Texture3DMapping() {}
        virtual Point3D map(const SurfacePoint &surfPt) const {
            TexCoord2D texCoord = surfPt.getTextureCoordinate();
            return Point3D(texCoord.u, texCoord.v, 0.0f);
        }
    };
    
    class SLR_API OffsetAndScale2DMapping : public Texture2DMapping {
        float m_offsetX, m_offsetY;
        float m_scaleX, m_scaleY;
    public:
        OffsetAndScale2DMapping(float ox, float oy, float sx, float sy) : m_offsetX(ox), m_offsetY(oy), m_scaleX(sx), m_scaleY(sy) {}
        Point3D map(const SurfacePoint &surfPt) const override {
            TexCoord2D texCoord = surfPt.getTextureCoordinate();
            return Point3D((texCoord.u + m_offsetX) * m_scaleX,
                           (texCoord.v + m_offsetY) * m_scaleY,
                           0.0f);
        }
    };
    
    class SLR_API WorldPosition3DMapping : public Texture3DMapping {
    public:
        WorldPosition3DMapping() {}
        Point3D map(const SurfacePoint &surfPt) const override {
            return surfPt.getPosition();
        }
    };
    
    
    
    class SLR_API SpectrumTexture {
    public:
        virtual ~SpectrumTexture() { }
        
        virtual SampledSpectrum evaluate(const SurfacePoint &surfPt, const WavelengthSamples &wls) const = 0;
        virtual RegularConstantContinuous2D* createIBLImportanceMap() const = 0;
    };
    
    class SLR_API Normal3DTexture {
    public:
        virtual ~Normal3DTexture() { }
        
        virtual Normal3D evaluate(const SurfacePoint &surfPt) const = 0;
        Normal3D evaluate(const TexCoord2D &tc) const {
            SurfacePoint surfPt;
            surfPt.setTextureCoordinate(tc);
            return evaluate(surfPt);
        }
    };
    
    class SLR_API FloatTexture {
    public:
        virtual ~FloatTexture() { }
        
        virtual float evaluate(const SurfacePoint &surfPt) const = 0;
        float evaluate(const TexCoord2D &tc) const {
            SurfacePoint surfPt;
            surfPt.setTextureCoordinate(tc);
            return evaluate(surfPt);
        }
    };
}

#endif

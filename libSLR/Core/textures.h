//
//  textures.h
//
//  Created by 渡部 心 on 2015/04/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_textures__
#define __SLR_textures__

#include "../defines.h"
#include "../declarations.h"
#include "geometry.h"

namespace SLR {
    class SLR_API Texture2DMapping {
    public:
        virtual ~Texture2DMapping() {}
        
        virtual Point3D map(const SurfacePoint &surfPt) const {
            TexCoord2D texCoord = surfPt.getTextureCoordinate();
            return Point3D(texCoord.u, texCoord.v, 0.0f);
        }
        virtual Point3D map(const MediumPoint &medPt) const {
            return medPt.getPosition();
        }
    };
    
    class SLR_API Texture3DMapping {
    public:
        virtual ~Texture3DMapping() {}
        
        virtual Point3D map(const SurfacePoint &surfPt) const {
            TexCoord2D texCoord = surfPt.getTextureCoordinate();
            return Point3D(texCoord.u, texCoord.v, 0.0f);
        }
        virtual Point3D map(const MediumPoint &medPt) const {
            return medPt.getPosition();
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
        virtual Point3D map(const MediumPoint &medPt) const override {
            Point3D pos = medPt.getPosition();
            return Point3D((pos.x + m_offsetX) * m_scaleX,
                           (pos.y + m_offsetY) * m_scaleY,
                           0.0f);
        }
    };
    
    class SLR_API OffsetAndScale3DMapping : public Texture3DMapping {
        float m_offsetX, m_offsetY, m_offsetZ;
        float m_scaleX, m_scaleY, m_scaleZ;
    public:
        OffsetAndScale3DMapping(float ox, float oy, float oz, float sx, float sy, float sz) :
        m_offsetX(ox), m_offsetY(oy), m_offsetZ(oz), m_scaleX(sx), m_scaleY(sy), m_scaleZ(sz) { }
        
        Point3D map(const SurfacePoint &surfPt) const override {
            TexCoord2D texCoord = surfPt.getTextureCoordinate();
            return Point3D((texCoord.u + m_offsetX) * m_scaleX,
                           (texCoord.v + m_offsetY) * m_scaleY,
                           0.0f);
        }
        virtual Point3D map(const MediumPoint &medPt) const override {
            Point3D pos = medPt.getPosition();
            return Point3D((pos.x + m_offsetX) * m_scaleX,
                           (pos.y + m_offsetY) * m_scaleY,
                           (pos.z + m_offsetZ) * m_scaleZ);
        }
    };
    
    class SLR_API WorldPosition3DMapping : public Texture3DMapping {
    public:
        WorldPosition3DMapping() {}
        
        Point3D map(const SurfacePoint &surfPt) const override {
            return surfPt.getPosition();
        }
        Point3D map(const MediumPoint &medPt) const override {
            return medPt.getPosition();
        }
    };
    
    
    
    class SLR_API SpectrumTexture {
    public:
        virtual ~SpectrumTexture() { }
        
        virtual SampledSpectrum evaluate(const SurfacePoint &surfPt, const WavelengthSamples &wls) const = 0;
        virtual SampledSpectrum evaluate(const MediumPoint &medPt, const WavelengthSamples &wls) const = 0;
        virtual RegularConstantContinuousDistribution2D* createIBLImportanceMap() const = 0;
    };
    
    class SLR_API NormalTexture {
    public:
        virtual ~NormalTexture() { }
        
        virtual Normal3D evaluate(const SurfacePoint &surfPt) const = 0;
        Normal3D evaluate(const TexCoord2D &tc) const {
            SurfacePoint surfPt;
            surfPt.setTextureCoordinate(tc);
            return evaluate(surfPt);
        }
        virtual Normal3D evaluate(const MediumPoint &medPt) const = 0;
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
        virtual float evaluate(const MediumPoint &medPt) const = 0;
    };
}

#endif /* __SLR_textures__ */

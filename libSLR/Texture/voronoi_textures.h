//
//  voronoi_textures.h
//
//  Created by 渡部 心 on 2016/06/17.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef __SLR_voronoi_textures__
#define __SLR_voronoi_textures__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/textures.h"

namespace SLR {
    // References
    // A Cellular Texture Basis Function
    
    class SLR_API VoronoiSpectrumTexture : public SpectrumTexture {
        const Texture3DMapping* m_mapping;
        float m_scale;
        float m_brightness;
    public:
        VoronoiSpectrumTexture(const Texture3DMapping* mapping, float scale, float brightness) :
        m_mapping(mapping), m_scale(scale), m_brightness(brightness) { }
        
        SampledSpectrum evaluate(const Point3D &p, const WavelengthSamples &wls) const;
        SampledSpectrum evaluate(const SurfacePoint &surfPt, const WavelengthSamples &wls) const override {
            return evaluate(m_mapping->map(surfPt) / m_scale, wls);
        }
        SampledSpectrum evaluate(const MediumPoint &medPt, const WavelengthSamples &wls) const override {
            return evaluate(m_mapping->map(medPt) / m_scale, wls);
        }
        const ContinuousDistribution2D* createIBLImportanceMap() const override;
    };
    
    
    
    class SLR_API VoronoiNormalTexture : public NormalTexture {
        const Texture3DMapping* m_mapping;
        float m_scale;
        float m_cosThetaMax;
    public:
        VoronoiNormalTexture(const Texture3DMapping* mapping, float scale, float thetaMax) :
        m_mapping(mapping), m_scale(scale), m_cosThetaMax(std::cos(thetaMax)) { }
        
        Normal3D evaluate(const Point3D &p) const;
        Normal3D evaluate(const SurfacePoint &surfPt) const override {
            return evaluate(m_mapping->map(surfPt) / m_scale);
        }
        Normal3D evaluate(const MediumPoint &medPt) const override {
            return evaluate(m_mapping->map(medPt) / m_scale);
        }
    };
    
    
    
    class SLR_API VoronoiFloatTexture : public FloatTexture {
        const Texture3DMapping* m_mapping;
        float m_scale;
        float m_valueScale;
        bool m_flat;
    public:
        VoronoiFloatTexture(const Texture3DMapping* mapping, float scale, float valueScale, bool flat) :
        m_mapping(mapping), m_scale(scale), m_valueScale(valueScale), m_flat(flat) { }
        
        float evaluate(const Point3D &p) const;
        float evaluate(const SurfacePoint &surfPt) const override {
            return evaluate(m_mapping->map(surfPt) / m_scale);
        }
        float evaluate(const MediumPoint &medPt) const override {
            return evaluate(m_mapping->map(medPt) / m_scale);
        }
    };
}

#endif /* __SLR_voronoi_textures__ */

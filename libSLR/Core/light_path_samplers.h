//
//  light_path_samplers.h
//
//  Created by 渡部 心 on 2016/09/06.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef __SLR_light_path_samplers__
#define __SLR_light_path_samplers__

#include "../defines.h"
#include "../references.h"
#include "cameras.h"
#include "directional_distribution_functions.h"
#include "SurfaceObject.h"
#include "MediumObject.h"

#include "../RNGs/XORShiftRNG.h"

namespace SLR {
    struct SLR_API PixelPosition {
        float x, y;
        PixelPosition(float xx, float yy) : x(xx), y(yy) { }
    };
    
    class SLR_API LightPathSampler {
    public:
        virtual ~LightPathSampler() { }
        
        virtual float getTimeSample(float timeBegin, float timeEnd) = 0;
        virtual PixelPosition getPixelPositionSample(uint32_t baseX, uint32_t baseY) = 0;
        virtual float getWavelengthSample() = 0;
        virtual float getWLSelectionSample() = 0;
        virtual LensPosSample getLensPosSample() = 0;
        virtual IDFSample getIDFSample() = 0;
        virtual BSDFSample getBSDFSample() = 0;
        virtual PFSample getPFSample() = 0;
        virtual FreePathSampler &getFreePathSampler() = 0;
        virtual float getPathTerminationSample() = 0;
        virtual float getLightSelectionSample() = 0;
        virtual SurfaceLightPosSample getSurfaceLightPosSample() = 0;
        virtual VolumetricLightPosSample getVolumetricLightPosSample() = 0;
        virtual EDFSample getEDFSample() = 0;
    };
    
    class SLR_API FreePathSampler {
        XORShiftRNG &m_rng;
    public:
        FreePathSampler(XORShiftRNG& rng) : m_rng(rng) { }
        
        float getSample() { return m_rng.getFloat0cTo1o(); }
    };
    
    
    
    class SLR_API IndependentLightPathSampler : public LightPathSampler {
        XORShiftRNG m_rng;
        FreePathSampler m_freePathSampler;
    public:
        IndependentLightPathSampler() :m_rng() , m_freePathSampler(m_rng) { }
        IndependentLightPathSampler(uint32_t seed) : m_rng(seed), m_freePathSampler(m_rng) { }
        
        float getTimeSample(float timeBegin, float timeEnd) override { float v = m_rng.getFloat0cTo1o(); return timeBegin * (1 - v) + timeEnd * v; }
        PixelPosition getPixelPositionSample(uint32_t baseX, uint32_t baseY) override { return PixelPosition(baseX + m_rng.getFloat0cTo1o(), baseY + m_rng.getFloat0cTo1o()); }
        float getWavelengthSample() override { return m_rng.getFloat0cTo1o(); }
        float getWLSelectionSample() override { return m_rng.getFloat0cTo1o(); }
        LensPosSample getLensPosSample() override { return LensPosSample(m_rng.getFloat0cTo1o(), m_rng.getFloat0cTo1o()); }
        IDFSample getIDFSample() override { return IDFSample(m_rng.getFloat0cTo1o(), m_rng.getFloat0cTo1o()); }
        BSDFSample getBSDFSample() override { return BSDFSample(m_rng.getFloat0cTo1o(), m_rng.getFloat0cTo1o(), m_rng.getFloat0cTo1o()); }
        PFSample getPFSample() override { return PFSample(m_rng.getFloat0cTo1o(), m_rng.getFloat0cTo1o()); }
        FreePathSampler &getFreePathSampler() override { return m_freePathSampler; }
        float getPathTerminationSample() override { return m_rng.getFloat0cTo1o(); }
        float getLightSelectionSample() override { return m_rng.getFloat0cTo1o(); }
        SurfaceLightPosSample getSurfaceLightPosSample() override { return SurfaceLightPosSample(m_rng.getFloat0cTo1o(), m_rng.getFloat0cTo1o()); }
        VolumetricLightPosSample getVolumetricLightPosSample() override { return VolumetricLightPosSample(m_rng.getFloat0cTo1o(), m_rng.getFloat0cTo1o(), m_rng.getFloat0cTo1o()); }
        EDFSample getEDFSample() override { return EDFSample(m_rng.getFloat0cTo1o(), m_rng.getFloat0cTo1o(), m_rng.getFloat0cTo1o()); }
    };
}

#endif /* __SLR_light_path_samplers__ */

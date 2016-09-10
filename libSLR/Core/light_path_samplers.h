//
//  light_path_samplers.h
//
//  Created by 渡部 心 on 2016/09/06.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef light_path_samplers_h
#define light_path_samplers_h

#include "../defines.h"
#include "../references.h"
#include "cameras.h"
#include "directional_distribution_functions.h"
#include "SurfaceObject.h"

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
        virtual BSSRDFSample getBSSRDFSample() = 0;
        virtual float getPathTerminationSample() = 0;
        virtual float getLightSelectionSample() = 0;
        virtual LightPosSample getLightPosSample() = 0;
        virtual EDFSample getEDFSample() = 0;
    };
    
    
    
    class SLR_API IndependentLightPathSampler : public LightPathSampler {
        XORShiftRNG m_rng;
    public:
        IndependentLightPathSampler() { }
        IndependentLightPathSampler(uint32_t seed) : m_rng(seed) { }
        
        float getTimeSample(float timeBegin, float timeEnd) override { float v = m_rng.getFloat0cTo1o(); return timeBegin * (1 - v) + timeEnd * v; }
        PixelPosition getPixelPositionSample(uint32_t baseX, uint32_t baseY) override { return PixelPosition(baseX + m_rng.getFloat0cTo1o(), baseY + m_rng.getFloat0cTo1o()); }
        float getWavelengthSample() override { return m_rng.getFloat0cTo1o(); }
        float getWLSelectionSample() override { return m_rng.getFloat0cTo1o(); }
        LensPosSample getLensPosSample() override { return LensPosSample(m_rng.getFloat0cTo1o(), m_rng.getFloat0cTo1o()); }
        IDFSample getIDFSample() override { return IDFSample(m_rng.getFloat0cTo1o(), m_rng.getFloat0cTo1o()); }
        BSDFSample getBSDFSample() override { return BSDFSample(m_rng.getFloat0cTo1o(), m_rng.getFloat0cTo1o(), m_rng.getFloat0cTo1o()); }
        BSSRDFSample getBSSRDFSample() override { return BSSRDFSample(m_rng.getFloat0cTo1o(), m_rng.getFloat0cTo1o(), m_rng.getFloat0cTo1o(), m_rng.getFloat0cTo1o(), m_rng.getFloat0cTo1o()); };
        float getPathTerminationSample() override { return m_rng.getFloat0cTo1o(); }
        float getLightSelectionSample() override { return m_rng.getFloat0cTo1o(); }
        LightPosSample getLightPosSample() override { return LightPosSample(m_rng.getFloat0cTo1o(), m_rng.getFloat0cTo1o()); }
        EDFSample getEDFSample() override { return EDFSample(m_rng.getFloat0cTo1o(), m_rng.getFloat0cTo1o(), m_rng.getFloat0cTo1o()); }
    };
}

#endif /* light_path_samplers_hpp */

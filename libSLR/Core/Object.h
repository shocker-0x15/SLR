//
//  Object.h
//
//  Created by 渡部 心 on 2016/09/19.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef __SLR_Object__
#define __SLR_Object__

#include "../defines.h"
#include "../references.h"
#include "geometry.h"
#include "directional_distribution_functions.h"

namespace SLR {
    struct SLR_API LightPosQuery {
        float time;
        WavelengthSamples wls;
        LightPosQuery(float t, const WavelengthSamples &lambdas) : time(t), wls(lambdas) { }
    };
    
    struct SLR_API LightPosQueryResult {
        DirectionType posType;
        
        virtual InteractionPoint* getInteractionPoint() = 0;
        virtual float spatialPDF() const = 0;
        DirectionType sampledPositionType() const { return posType; }
    };
    
    class SLR_API Light {
    protected:
        StaticTransform m_appliedTransform;
    public:
        void applyTransformFromLeft(const StaticTransform &transform) {
            m_appliedTransform = transform * m_appliedTransform;
        }
        
        virtual SampledSpectrum sample(const LightPosQuery &query, LightPathSampler &pathSampler, ArenaAllocator &mem, LightPosQueryResult** lpResult) const = 0;
        virtual Ray sampleRay(const LightPosQuery &lightPosQuery, LightPathSampler &pathSampler, const EDFQuery &edfQuery, ArenaAllocator &mem,
                              LightPosQueryResult** lightPosResult, SampledSpectrum* Le0, EDF** edf,
                              EDFQueryResult* edfResult, SampledSpectrum* Le1) const= 0;
    };
    
    
    
    class SLR_API Object {
    public:
        virtual ~Object() { }
        
        virtual BoundingBox3D bounds() const = 0;
        virtual BoundingBox3D choppedBounds(BoundingBox3D::Axis chopAxis, float minChopPos, float maxChopPos) const {
            BoundingBox3D baseBBox = bounds();
            if (maxChopPos < baseBBox.minP[chopAxis])
                return BoundingBox3D();
            if (minChopPos > baseBBox.maxP[chopAxis])
                return BoundingBox3D();
            if (minChopPos < baseBBox.minP[chopAxis] && maxChopPos > baseBBox.maxP[chopAxis])
                return baseBBox;
            BoundingBox3D ret = baseBBox;
            ret.minP[chopAxis] = std::max(minChopPos, ret.minP[chopAxis]);
            ret.maxP[chopAxis] = std::min(maxChopPos, ret.maxP[chopAxis]);
            return ret;
        }
        virtual void splitBounds(BoundingBox3D::Axis splitAxis, float splitPos, BoundingBox3D* bbox0, BoundingBox3D* bbox1) const {
            BoundingBox3D baseBBox = bounds();
            if (splitPos < baseBBox.minP[splitAxis]) {
                *bbox0 = BoundingBox3D();
                *bbox1 = baseBBox;
                return;
            }
            if (splitPos > baseBBox.maxP[splitAxis]) {
                *bbox0 = baseBBox;
                *bbox1 = BoundingBox3D();
                return;
            }
            *bbox0 = baseBBox;
            bbox0->maxP[splitAxis] = std::min(bbox0->maxP[splitAxis], splitPos);
            *bbox1 = baseBBox;
            bbox1->minP[splitAxis] = std::max(bbox1->minP[splitAxis], splitPos);
        }
    };
}

#endif /* __SLR_Object__ */

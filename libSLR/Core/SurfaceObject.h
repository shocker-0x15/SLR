//
//  SurfaceObject.h
//
//  Created by 渡部 心 on 2015/07/15.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__SurfaceObject__
#define __SLR__SurfaceObject__

#include "../defines.h"
#include "../references.h"
#include "Object.h"

namespace SLR {
    class SLR_API SurfaceObject : public Object {
    public:
        SurfaceObject() { }
        virtual ~SurfaceObject() { }
        
        virtual float costForIntersect() const = 0;
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
        virtual bool intersect(Ray &ray, Intersection* isect) const = 0;
        virtual Point3D getIntersectionPoint(const Intersection &isect) const { return isect.obj.back()->getIntersectionPoint(isect); }
        virtual const SurfaceMaterial* getSurfaceMaterial() const { SLRAssert_NotImplemented(); return nullptr; }
        virtual void getSurfacePoint(const Intersection &isect, SurfacePoint* surfPt) const { isect.obj.back()->getSurfacePoint(isect, surfPt); }
        
        bool intersect(Ray &ray, SurfacePoint* surfPt) const;
        bool testVisibility(const SurfacePoint &shdP, const SurfacePoint &lightP, float time) const;
    };
    
    
    
    class SLR_API SingleSurfaceObject : public SurfaceObject {
    protected:
        const Surface* m_surface;
        const SurfaceMaterial* m_material;
    public:
        SingleSurfaceObject() { }
        SingleSurfaceObject(const Surface* surf, const SurfaceMaterial* mat) : m_surface(surf), m_material(mat) { }
        virtual ~SingleSurfaceObject() { }
        
        float costForIntersect() const override { return m_surface->costForIntersect(); }
        BoundingBox3D bounds() const override { return m_surface->bounds(); }
        BoundingBox3D choppedBounds(BoundingBox3D::Axis chopAxis, float minChopPos, float maxChopPos) const override {
            return m_surface->choppedBounds(chopAxis, minChopPos, maxChopPos);
        }
        void splitBounds(BoundingBox3D::Axis chopAxis, float splitPos, BoundingBox3D* bbox0, BoundingBox3D* bbox1) const override {
            m_surface->splitBounds(chopAxis, splitPos, bbox0, bbox1);
        }
        bool intersect(Ray &ray, Intersection* isect) const override;
        Point3D getIntersectionPoint(const Intersection &isect) const override { return isect.p; }
        const SurfaceMaterial* getSurfaceMaterial() const override { return m_material; }
        void getSurfacePoint(const Intersection &isect, SurfacePoint* surfPt) const override;
        
        bool isEmitting() const override;
        float importance() const override;
        void selectLight(float u, Light* light, float* prob) const override;
        float evaluateProb(const Light &light) const override;
        
        SampledSpectrum sample(const Light &light, const LightPosQuery &query, const LightPosSample &smp, LightPosQueryResult* result) const override;
        Ray sampleRay(const Light &light,
                      const LightPosQuery &lightPosQuery, const LightPosSample &lightPosSample, LightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                      const EDFQuery &edfQuery, const EDFSample &edfSample, EDFQueryResult* edfResult, SampledSpectrum* Le1,
                      ArenaAllocator &mem) const override;
        
        virtual BSDF* createBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const;
        
        virtual float evaluateAreaPDF(const SurfacePoint& surfPt) const;
        virtual SampledSpectrum emittance(const SurfacePoint& surfPt, const WavelengthSamples &wls) const;
        virtual EDF* createEDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const;
    };
    
    class SLR_API BumpSingleSurfaceObject : public SingleSurfaceObject {
        const Normal3DTexture* m_normalMap;
    public:
        BumpSingleSurfaceObject(const Surface* surf, const SurfaceMaterial* mat, const Normal3DTexture* normalMap) :
        SingleSurfaceObject(surf, mat), m_normalMap(normalMap) { }
        
        void getSurfacePoint(const Intersection &isect, SurfacePoint* surfPt) const override;
    };
    
    class SLR_API InfiniteSphereSurfaceObject : public SingleSurfaceObject {
        const Scene* m_scene;
        const RegularConstantContinuous2D* m_dist;
    public:
        InfiniteSphereSurfaceObject(const Scene* scene, const IBLEmission* emitter);
        ~InfiniteSphereSurfaceObject();
        
        bool isEmitting() const override;
        float importance() const override;
        
        SampledSpectrum sample(const Light &light, const LightPosQuery &query, const LightPosSample &smp, LightPosQueryResult* result) const override;
        Ray sampleRay(const Light &light,
                      const LightPosQuery &lightPosQuery, const LightPosSample &lightPosSample, LightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                      const EDFQuery &edfQuery, const EDFSample &edfSample, EDFQueryResult* edfResult, SampledSpectrum* Le1,
                      ArenaAllocator &mem) const override;
        
        BSDF* createBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const override;
        
        float evaluateAreaPDF(const SurfacePoint& surfPt) const override;
    };
    
    
    
    class SLR_API SurfaceObjectAggregate : public SurfaceObject {
        Accelerator* m_accelerator;
        const Object** m_lightList;
        RegularConstantDiscrete1D* m_lightDist1D;
        std::map<const Object*, uint32_t> m_revMap;
    public:
        SurfaceObjectAggregate(std::vector<SurfaceObject*> &objs);
        ~SurfaceObjectAggregate();
        
        float costForIntersect() const override;
        BoundingBox3D bounds() const override;
        bool intersect(Ray &ray, Intersection* isect) const override;
        
        bool isEmitting() const override;
        float importance() const override;
        void selectLight(float u, Light* light, float* prob) const override;
        float evaluateProb(const Light &light) const override;
        
        SampledSpectrum sample(const Light &light, const LightPosQuery &query, const LightPosSample &smp, LightPosQueryResult* result) const override {
            SLRAssert(false, "SurfaceObjectAggregate::sample() should not be called.");
            return SampledSpectrum::Zero;
        }
        Ray sampleRay(const Light &light,
                      const LightPosQuery &lightPosQuery, const LightPosSample &lightPosSample, LightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                      const EDFQuery &edfQuery, const EDFSample &edfSample, EDFQueryResult* edfResult, SampledSpectrum* Le1,
                      ArenaAllocator &mem) const override {
            SLRAssert(false, "SurfaceObjectAggregate::sampleRay() should not be called.");
            return Ray();
        }
    };
    
    
    
    class SLR_API TransformedSurfaceObject : public SurfaceObject {
        const SurfaceObject* m_surfObj;
        const Transform* m_transform;
        friend class Light;
    public:
        TransformedSurfaceObject(const SurfaceObject* surfObj, const Transform* transform) : m_surfObj(surfObj), m_transform(transform) { }
        
        float costForIntersect() const override { return m_surfObj->costForIntersect(); }
        BoundingBox3D bounds() const override;
        bool intersect(Ray &ray, Intersection* isect) const override;
        Point3D getIntersectionPoint(const Intersection &isect) const override;
        void getSurfacePoint(const Intersection &isect, SurfacePoint* surfPt) const override;
        
        bool isEmitting() const override;
        float importance() const override;
        void selectLight(float u, Light* light, float* prob) const override;
        float evaluateProb(const Light &light) const override;
        
        SampledSpectrum sample(const Light &light, const LightPosQuery &query, const LightPosSample &smp, LightPosQueryResult* result) const override;
        Ray sampleRay(const Light &light,
                      const LightPosQuery &lightPosQuery, const LightPosSample &lightPosSample, LightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                      const EDFQuery &edfQuery, const EDFSample &edfSample, EDFQueryResult* edfResult, SampledSpectrum* Le1,
                      ArenaAllocator &mem) const override;
        
        const SurfaceMaterial* getSurfaceMaterial() const override { return m_surfObj->getSurfaceMaterial(); }
        
        void setTransform(const Transform* t) { m_transform = t; }
    };
}

#endif

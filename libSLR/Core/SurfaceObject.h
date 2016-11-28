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
        
        //----------------------------------------------------------------
        // Object's methods
        SampledSpectrum sample(const Light &light, const LightPosQuery &query, const LightPosSample &smp, LightPosQueryResult* result) const override {
            SLRAssert(light.top() == this, "Object stored in Intersection does not match this object.");
            ScopedPop<const Object*> sp = light.scopedPop();
            return light.top()->sample(light, query, smp, result);
        }
        Ray sampleRay(const Light &light,
                      const LightPosQuery &lightPosQuery, const LightPosSample &lightPosSample, LightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                      const EDFQuery &edfQuery, const EDFSample &edfSample, EDFQueryResult* edfResult, SampledSpectrum* Le1,
                      ArenaAllocator &mem) const override {
            SLRAssert(light.top() == this, "Object stored in Intersection does not match this object.");
            ScopedPop<const Object*> sp = light.scopedPop();
            return light.top()->sampleRay(light, lightPosQuery, lightPosSample, lightPosResult, Le0, edf, edfQuery, edfSample, edfResult, Le1, mem);
        }
        //----------------------------------------------------------------
        
        virtual float costForIntersect() const = 0;
        virtual bool intersect(Ray &ray, Intersection* isect) const = 0;
        virtual const SurfaceMaterial* getSurfaceMaterial() const { SLRAssert_NotImplemented(); return nullptr; }
        virtual void getSurfacePoint(const Intersection &isect, SurfacePoint* surfPt) const { SLRAssert_NotImplemented(); }
        
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
        
        //----------------------------------------------------------------
        // Object's methods
        BoundingBox3D bounds() const override { return m_surface->bounds(); }
        BoundingBox3D choppedBounds(BoundingBox3D::Axis chopAxis, float minChopPos, float maxChopPos) const override {
            return m_surface->choppedBounds(chopAxis, minChopPos, maxChopPos);
        }
        void splitBounds(BoundingBox3D::Axis chopAxis, float splitPos, BoundingBox3D* bbox0, BoundingBox3D* bbox1) const override {
            m_surface->splitBounds(chopAxis, splitPos, bbox0, bbox1);
        }
        
        bool isEmitting() const override;
        float importance() const override;
        void selectLight(float u, Light* light, float* prob) const override;
        float evaluateProb(const Light &light) const override;
        
        SampledSpectrum sample(const Light &light, const LightPosQuery &query, const LightPosSample &smp, LightPosQueryResult* result) const override;
        Ray sampleRay(const Light &light,
                      const LightPosQuery &lightPosQuery, const LightPosSample &lightPosSample, LightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                      const EDFQuery &edfQuery, const EDFSample &edfSample, EDFQueryResult* edfResult, SampledSpectrum* Le1,
                      ArenaAllocator &mem) const override;
        //----------------------------------------------------------------
        
        //----------------------------------------------------------------
        // SurfaceObject's methods
        float costForIntersect() const override { return m_surface->costForIntersect(); }
        bool intersect(Ray &ray, Intersection* isect) const override;
        const SurfaceMaterial* getSurfaceMaterial() const override { return m_material; }
        void getSurfacePoint(const Intersection &isect, SurfacePoint* surfPt) const override;
        //----------------------------------------------------------------
        
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
        
        //----------------------------------------------------------------
        // SurfaceObject's methods
        void getSurfacePoint(const Intersection &isect, SurfacePoint* surfPt) const override;
        //----------------------------------------------------------------
    };
    
    class SLR_API InfiniteSphereSurfaceObject : public SingleSurfaceObject {
        const Scene* m_scene;
        const RegularConstantContinuous2D* m_dist;
    public:
        InfiniteSphereSurfaceObject(const Scene* scene, const IBLEmission* emitter);
        ~InfiniteSphereSurfaceObject();
        
        //----------------------------------------------------------------
        // Object's methods
        bool isEmitting() const override;
        float importance() const override;
        
        SampledSpectrum sample(const Light &light, const LightPosQuery &query, const LightPosSample &smp, LightPosQueryResult* result) const override;
        Ray sampleRay(const Light &light,
                      const LightPosQuery &lightPosQuery, const LightPosSample &lightPosSample, LightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                      const EDFQuery &edfQuery, const EDFSample &edfSample, EDFQueryResult* edfResult, SampledSpectrum* Le1,
                      ArenaAllocator &mem) const override;
        //----------------------------------------------------------------
        
        //----------------------------------------------------------------
        // SingleSurfaceObject's methods
        BSDF* createBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const override;
        
        float evaluateAreaPDF(const SurfacePoint& surfPt) const override;
        //----------------------------------------------------------------
    };
    
    
    
    class SLR_API SurfaceObjectAggregate : public SurfaceObject {
        Accelerator* m_accelerator;
        const Object** m_lightList;
        RegularConstantDiscrete1D* m_lightDist1D;
        std::map<const Object*, uint32_t> m_revMap;
    public:
        SurfaceObjectAggregate(std::vector<SurfaceObject*> &objs);
        ~SurfaceObjectAggregate();
        
        //----------------------------------------------------------------
        // Object's methods
        BoundingBox3D bounds() const override;
        
        bool isEmitting() const override;
        float importance() const override;
        void selectLight(float u, Light* light, float* prob) const override;
        float evaluateProb(const Light &light) const override;
        //----------------------------------------------------------------
        
        //----------------------------------------------------------------
        // SurfaceObject's methods
        float costForIntersect() const override;
        bool intersect(Ray &ray, Intersection* isect) const override;
        void getSurfacePoint(const Intersection &isect, SurfacePoint* surfPt) const override {
            SLRAssert(isect.top() == this, "Object stored in Intersection does not match this object.");
            ScopedPop<const SurfaceObject*> sp = isect.scopedPop();
            isect.top()->getSurfacePoint(isect, surfPt);
        }
        //----------------------------------------------------------------
    };
    
    
    
    class SLR_API TransformedSurfaceObject : public SurfaceObject {
        const SurfaceObject* m_surfObj;
        const Transform* m_transform;
        friend class Light;
    public:
        TransformedSurfaceObject(const SurfaceObject* surfObj, const Transform* transform) : m_surfObj(surfObj), m_transform(transform) { }
        
        //----------------------------------------------------------------
        // Object's methods
        BoundingBox3D bounds() const override;
        
        bool isEmitting() const override;
        float importance() const override;
        void selectLight(float u, Light* light, float* prob) const override;
        float evaluateProb(const Light &light) const override;
        
        SampledSpectrum sample(const Light &light, const LightPosQuery &query, const LightPosSample &smp, LightPosQueryResult* result) const override;
        Ray sampleRay(const Light &light,
                      const LightPosQuery &lightPosQuery, const LightPosSample &lightPosSample, LightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                      const EDFQuery &edfQuery, const EDFSample &edfSample, EDFQueryResult* edfResult, SampledSpectrum* Le1,
                      ArenaAllocator &mem) const override;
        //----------------------------------------------------------------
        
        //----------------------------------------------------------------
        // SurfaceObject's methods
        float costForIntersect() const override { return m_surfObj->costForIntersect(); }
        bool intersect(Ray &ray, Intersection* isect) const override;
        const SurfaceMaterial* getSurfaceMaterial() const override { return m_surfObj->getSurfaceMaterial(); }
        void getSurfacePoint(const Intersection &isect, SurfacePoint* surfPt) const override;
        //----------------------------------------------------------------
        
        void setTransform(const Transform* t) { m_transform = t; }
    };
}

#endif

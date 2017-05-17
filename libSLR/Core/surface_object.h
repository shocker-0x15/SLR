//
//  surface_object.h
//
//  Created by 渡部 心 on 2015/07/15.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_surface_object__
#define __SLR_surface_object__

#include "../defines.h"
#include "../declarations.h"
#include "object.h"

namespace SLR {
    struct SLR_API SurfaceLightPosSample {
        float uPos[2];
        SurfaceLightPosSample(float up0, float up1) : uPos{up0, up1} { }
    };
    
    struct SLR_API SurfaceLightPosQueryResult : public LightPosQueryResult {
        SurfacePoint surfPt;
        float areaPDF;
        
        InteractionPoint* getInteractionPoint() override { return &surfPt; }
        float spatialPDF() const override { return areaPDF; }
    };
    
    class SLR_API SurfaceLight : public Light {
        const SingleSurfaceObject* m_obj;
    public:
        SurfaceLight() { }
        
        void setObject(const SingleSurfaceObject* obj) { m_obj = obj; }
        
        SampledSpectrum sample(const LightPosQuery &query, const SurfaceLightPosSample &smp, SurfaceLightPosQueryResult* result) const;
        void sampleRay(const LightPosQuery &lightPosQuery, const SurfaceLightPosSample &lightPosSample,
                       const EDFQuery &edfQuery, const EDFSample &edfSample,
                       ArenaAllocator &mem,
                       SurfaceLightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                       EDFQueryResult* edfResult, SampledSpectrum* Le1, Ray* ray, float* epsilon) const;
        
        // ----------------------------------------------------------------
        // Light's methods
        
        SampledSpectrum sample(const LightPosQuery &query, LightPathSampler &pathSampler, ArenaAllocator &mem, LightPosQueryResult** lpResult) const override;
        void sampleRay(const LightPosQuery &lightPosQuery, LightPathSampler &pathSampler, const EDFQuery &edfQuery, ArenaAllocator &mem,
                       LightPosQueryResult** lightPosResult, SampledSpectrum* Le0, EDF** edf,
                       EDFQueryResult* edfResult, SampledSpectrum* Le1, Ray* ray, float* epsilon) const override;
        
        // END: Light's methods
        // ----------------------------------------------------------------
    };
    
    
    
    class SLR_API SurfaceObject : public Object {
    public:
        virtual bool isEmitting() const = 0;
        virtual float importance() const = 0;
        virtual void selectLight(float u, float time, SurfaceLight* light, float* prob) const = 0;
        
        virtual SampledSpectrum sample(const StaticTransform &transform,
                                       const LightPosQuery &query, const SurfaceLightPosSample &smp, SurfaceLightPosQueryResult* result) const {
            SLRAssert_ShouldNotBeCalled();
            return SampledSpectrum::Zero;
        }
        virtual void sampleRay(const StaticTransform &transform,
                               const LightPosQuery &lightPosQuery, const SurfaceLightPosSample &lightPosSample,
                               const EDFQuery &edfQuery, const EDFSample &edfSample,
                               ArenaAllocator &mem,
                               SurfaceLightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                               EDFQueryResult* edfResult, SampledSpectrum* Le1, Ray* ray, float* epsilon) const {
            SLRAssert_ShouldNotBeCalled();
        }
        
        virtual float costForIntersect() const = 0;
        virtual bool contains(const Point3D &p, float time) const { return false; }
        virtual bool intersect(const Ray &ray, const RaySegment &segment, SurfaceInteraction* si) const = 0;
        virtual void calculateSurfacePoint(const SurfaceInteraction &si, SurfacePoint* surfPt) const {
            SLRAssert_ShouldNotBeCalled();
        }
        
        bool testVisibility(const SurfacePoint &shdP, const SurfacePoint &lightP, float time) const;
    };
    
    
    
    class SLR_API SingleSurfaceObject : public SurfaceObject {
    protected:
        const SurfaceShape* m_surface;
        const SurfaceMaterial* m_material;
    public:
        SingleSurfaceObject() { }
        SingleSurfaceObject(const SurfaceShape* surf, const SurfaceMaterial* mat) : m_surface(surf), m_material(mat) { }
        virtual ~SingleSurfaceObject() { }
        
        virtual BSDF* createBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const;
        
        virtual float evaluateAreaPDF(const SurfacePoint& surfPt) const;
        virtual SampledSpectrum emittance(const SurfacePoint& surfPt, const WavelengthSamples &wls) const;
        virtual EDF* createEDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const;
        
        // ----------------------------------------------------------------
        // Object's methods
        
        BoundingBox3D bounds() const override { return m_surface->bounds(); }
        BoundingBox3D choppedBounds(BoundingBox3D::Axis chopAxis, float minChopPos, float maxChopPos) const override {
            return m_surface->choppedBounds(chopAxis, minChopPos, maxChopPos);
        }
        void splitBounds(BoundingBox3D::Axis chopAxis, float splitPos, BoundingBox3D* bbox0, BoundingBox3D* bbox1) const override {
            m_surface->splitBounds(chopAxis, splitPos, bbox0, bbox1);
        }
        
        // END: Object's methods
        // ----------------------------------------------------------------
        
        // ----------------------------------------------------------------
        // SurfaceObject's methods
        
        bool isEmitting() const override;
        float importance() const override;
        void selectLight(float u, float time, SurfaceLight* light, float* prob) const override;
        
        SampledSpectrum sample(const StaticTransform &transform,
                               const LightPosQuery &query, const SurfaceLightPosSample &smp, SurfaceLightPosQueryResult* result) const override;
        void sampleRay(const StaticTransform &transform,
                       const LightPosQuery &lightPosQuery, const SurfaceLightPosSample &lightPosSample,
                       const EDFQuery &edfQuery, const EDFSample &edfSample,
                       ArenaAllocator &mem,
                       SurfaceLightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                       EDFQueryResult* edfResult, SampledSpectrum* Le1, Ray* ray, float* epsilon) const override;
        
        float costForIntersect() const override { return m_surface->costForIntersect(); }
        bool intersect(const Ray &ray, const RaySegment &segment, SurfaceInteraction* si) const override;
        void calculateSurfacePoint(const SurfaceInteraction &si, SurfacePoint* surfPt) const override;
        
        // END: SurfaceObject's methods
        // ----------------------------------------------------------------
    };
    
    
    
    class SLR_API BumpSingleSurfaceObject : public SingleSurfaceObject {
        const NormalTexture* m_normalMap;
    public:
        BumpSingleSurfaceObject(const SurfaceShape* surf, const SurfaceMaterial* mat, const NormalTexture* normalMap) :
        SingleSurfaceObject(surf, mat), m_normalMap(normalMap) { }
        
        // ----------------------------------------------------------------
        // SurfaceObject's methods
        
        void calculateSurfacePoint(const SurfaceInteraction &si, SurfacePoint* surfPt) const override;
        
        // END: SurfaceObject's methods
        // ----------------------------------------------------------------
    };
    
    
    
    class SLR_API InfiniteSphereSurfaceObject : public SingleSurfaceObject {
        const Scene* m_scene;
        const ContinuousDistribution2D* m_dist;
    public:
        InfiniteSphereSurfaceObject(const Scene* scene, const IBLEmitterSurfaceProperty* emitter);
        ~InfiniteSphereSurfaceObject();
        
        // ----------------------------------------------------------------
        // SurfaceObject's methods
        
        bool isEmitting() const override;
        float importance() const override;
        
        SampledSpectrum sample(const StaticTransform &transform,
                               const LightPosQuery &query, const SurfaceLightPosSample &smp, SurfaceLightPosQueryResult* result) const override;
        void sampleRay(const StaticTransform &transform,
                       const LightPosQuery &lightPosQuery, const SurfaceLightPosSample &lightPosSample,
                       const EDFQuery &edfQuery, const EDFSample &edfSample,
                       ArenaAllocator &mem,
                       SurfaceLightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                       EDFQueryResult* edfResult, SampledSpectrum* Le1, Ray* ray, float* epsilon) const override;
        
        // END: SurfaceObject's methods
        // ----------------------------------------------------------------
        
        // ----------------------------------------------------------------
        // SingleSurfaceObject's methods
        
        BSDF* createBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const override;
        
        float evaluateAreaPDF(const SurfacePoint& surfPt) const override;
        
        // END: SingleSurfaceObject's methods
        // ----------------------------------------------------------------
    };
    
    
    
    class SLR_API TransformedSurfaceObject : public SurfaceObject {
        const SurfaceObject* m_surfObj;
        const Transform* m_transform;
        friend class Light;
    public:
        TransformedSurfaceObject(const SurfaceObject* surfObj, const Transform* transform) : m_surfObj(surfObj), m_transform(transform) { }
        
        void setTransform(const Transform* t) { m_transform = t; }
        
        // ----------------------------------------------------------------
        // Object's methods
        
        BoundingBox3D bounds() const override;
        
        // END: Object's methods
        // ----------------------------------------------------------------
        
        // ----------------------------------------------------------------
        // SurfaceObject's methods
        
        bool isEmitting() const override;
        float importance() const override;
        void selectLight(float u, float time, SurfaceLight* light, float* prob) const override;
        
        float costForIntersect() const override { return m_surfObj->costForIntersect(); }
        bool contains(const Point3D &p, float time) const override;
        bool intersect(const Ray &ray, const RaySegment &segment, SurfaceInteraction* si) const override;
        
        // END: SurfaceObject's methods
        // ----------------------------------------------------------------
    };
    
    
    
    class SLR_API SurfaceObjectAggregate : public SurfaceObject {
        Accelerator* m_accelerator;
        const SurfaceObject** m_lightList;
        std::map<uint32_t, uint32_t> m_objToLightMap;
        uint32_t m_numLights;
        DiscreteDistribution1D* m_lightDist1D;
    public:
        SurfaceObjectAggregate(std::vector<SurfaceObject*> &objs);
        ~SurfaceObjectAggregate();
        
        // ----------------------------------------------------------------
        // Object's methods
        
        BoundingBox3D bounds() const override;
        
        // END: Object's methods
        // ----------------------------------------------------------------
        
        // ----------------------------------------------------------------
        // SurfaceObject's methods
        
        bool isEmitting() const override;
        float importance() const override;
        void selectLight(float u, float time, SurfaceLight* light, float* prob) const override;
        
        float costForIntersect() const override;
        bool contains(const Point3D &p, float time) const override;
        bool intersect(const Ray &ray, const RaySegment &segment, SurfaceInteraction* si) const override;
        
        // END: SurfaceObject's methods
        // ----------------------------------------------------------------
    };
}

#endif /* __SLR_surface_object__ */

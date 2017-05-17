//
//  medium_object.h
//
//  Created by 渡部 心 on 2016/09/19.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef __SLR_medium_object__
#define __SLR_medium_object__

#include "../defines.h"
#include "../declarations.h"
#include "object.h"

namespace SLR {
    struct SLR_API VolumetricLightPosSample {
        float uPos[3];
        VolumetricLightPosSample(float up0, float up1, float up2) : uPos{up0, up1, up2} { }
    };
    
    struct SLR_API VolumetricLightPosQueryResult : public LightPosQueryResult {
        MediumPoint medPt;
        float volumePDF;
        
        InteractionPoint* getInteractionPoint() override { return &medPt; }
        float spatialPDF() const override { return volumePDF; }
    };
    
    class SLR_API VolumetricLight : public Light {
        const SingleMediumObject* m_obj;
    public:
        VolumetricLight() { }
        
        void setObject(const SingleMediumObject* obj) { m_obj = obj; }
        
        SampledSpectrum sample(const LightPosQuery &query, const VolumetricLightPosSample &smp, VolumetricLightPosQueryResult* result) const;
        void sampleRay(const LightPosQuery &lightPosQuery, const VolumetricLightPosSample &lightPosSample,
                       const EDFQuery &edfQuery, const EDFSample &edfSample,
                       ArenaAllocator &mem,
                       VolumetricLightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
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
    
    
    
    class SLR_API MediumObject : public Object {
    public:
        virtual bool isEmitting() const = 0;
        virtual float importance() const = 0;
        virtual void selectLight(float u, float time, VolumetricLight* light, float* prob) const = 0;
        
        virtual SampledSpectrum sample(const StaticTransform &transform,
                                       const LightPosQuery &query, const VolumetricLightPosSample &smp, VolumetricLightPosQueryResult* result) const {
            SLRAssert_ShouldNotBeCalled();
            return SampledSpectrum::Zero;
        }
        virtual void sampleRay(const StaticTransform &transform,
                               const LightPosQuery &lightPosQuery, const VolumetricLightPosSample &lightPosSample,
                               const EDFQuery &edfQuery, const EDFSample &edfSample,
                               ArenaAllocator &mem,
                               VolumetricLightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                               EDFQueryResult* edfResult, SampledSpectrum* Le1, Ray* ray, float* epsilon) const {
            SLRAssert_ShouldNotBeCalled();
        }
        
        virtual bool contains(const Point3D &p, float time) const = 0;
        virtual bool intersectBoundary(const Ray &ray, const RaySegment &segment, float* distToBoundary, bool* enter) const = 0;
        virtual bool interact(const Ray &ray, const RaySegment &segment, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                              MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const = 0;
        virtual SampledSpectrum evaluateTransmittance(const Ray &ray, const RaySegment &segment, const WavelengthSamples &wls, LightPathSampler &pathSampler, bool* singleWavelength) const = 0;
        virtual void calculateMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const {
            SLRAssert_ShouldNotBeCalled();
        }
    };
    
    
    
    class SLR_API SingleMediumObject : public MediumObject {
        const MediumDistribution* m_medium;
        const MediumMaterial* m_material;
    public:
        SingleMediumObject(const MediumDistribution* medium, const MediumMaterial* material) :
        m_medium(medium), m_material(material) {}
        
        virtual SampledSpectrum evaluateExtinctionCoefficient(const MediumPoint &medPt, const WavelengthSamples &wls) const;
        virtual AbstractBDF* createAbstractBDF(const MediumPoint &medPt, const WavelengthSamples &wls, ArenaAllocator &mem) const;
        
        // ----------------------------------------------------------------
        // Object's methods
        
        BoundingBox3D bounds() const override { return m_medium->bounds(); }
        
        // END: Object's methods
        // ----------------------------------------------------------------
        
        // ----------------------------------------------------------------
        // MediumObject's methods
        
        bool isEmitting() const override;
        float importance() const override;
        void selectLight(float u, float time, VolumetricLight* light, float* prob) const override;
        
        SampledSpectrum sample(const StaticTransform &transform,
                               const LightPosQuery &query, const VolumetricLightPosSample &smp, VolumetricLightPosQueryResult* result) const override;
        void sampleRay(const StaticTransform &transform,
                       const LightPosQuery &lightPosQuery, const VolumetricLightPosSample &lightPosSample,
                       const EDFQuery &edfQuery, const EDFSample &edfSample,
                       ArenaAllocator &mem,
                       VolumetricLightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                       EDFQueryResult* edfResult, SampledSpectrum* Le1, Ray* ray, float* epsilon) const override;

        bool contains(const Point3D &p, float time) const override { return m_medium->contains(p); }
        bool intersectBoundary(const Ray &ray, const RaySegment &segment, float* distToBoundary, bool* enter) const override {
            return m_medium->intersectBoundary(ray, segment, distToBoundary, enter);
        }
        bool interact(const Ray &ray, const RaySegment &segment, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                      MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const override;
        SampledSpectrum evaluateTransmittance(const Ray &ray, const RaySegment &segment, const WavelengthSamples &wls, LightPathSampler &pathSampler, bool* singleWavelength) const override;
        void calculateMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const override;
        
        // END: MediumObject's methods
        // ----------------------------------------------------------------
    };
    
    
    
    class SLR_API TransformedMediumObject : public MediumObject {
        const MediumObject* m_medObj;
        const Transform* m_transform;
        friend class Light;
    public:
        TransformedMediumObject(const MediumObject* medObj, const Transform* transform) : m_medObj(medObj), m_transform(transform) { }
        
        // ----------------------------------------------------------------
        // Object's methods
        
        BoundingBox3D bounds() const override;
        
        // END: Object's methods
        // ----------------------------------------------------------------
        
        // ----------------------------------------------------------------
        // MediumObject's methods
        
        bool isEmitting() const override;
        float importance() const override;
        void selectLight(float u, float time, VolumetricLight* light, float* prob) const override;
        
        bool contains(const Point3D &p, float time) const override;
        bool intersectBoundary(const Ray &ray, const RaySegment &segment, float* distToBoundary, bool* enter) const override;
        bool interact(const Ray &ray, const RaySegment &segment, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                      MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const override;
        SampledSpectrum evaluateTransmittance(const Ray &ray, const RaySegment &segment, const WavelengthSamples &wls, LightPathSampler &pathSampler, bool* singleWavelength) const override;
        
        // END: MediumObject's methods
        // ----------------------------------------------------------------
        
        void setTransform(const Transform* t) { m_transform = t; }
    };
    
    
    
    class SLR_API EnclosedMediumObject : public MediumObject {
        const MediumObject* m_medObj;
        const SurfaceObject* m_boundary;
        const StaticTransform m_medToSurfTF;
    public:
        EnclosedMediumObject(const MediumObject* medObj, const SurfaceObject* boundary, const StaticTransform medToSurfTF) :
        m_medObj(medObj), m_boundary(boundary), m_medToSurfTF(medToSurfTF) { }
        
        // ----------------------------------------------------------------
        // Object's methods
        
        BoundingBox3D bounds() const override;
        
        // END: Object's methods
        // ----------------------------------------------------------------
        
        // ----------------------------------------------------------------
        // MediumObject's methods
        
        bool isEmitting() const override { return m_medObj->isEmitting(); }
        float importance() const override { return m_medObj->importance(); }
        void selectLight(float u, float time, VolumetricLight* light, float* prob) const override {
            m_medObj->selectLight(u, time, light, prob);
        }
        
        bool contains(const Point3D &p, float time) const override;
        bool intersectBoundary(const Ray &ray, const RaySegment &segment, float* distToBoundary, bool* enter) const override;
        bool interact(const Ray &ray, const RaySegment &segment, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                      MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const override;
        SampledSpectrum evaluateTransmittance(const Ray &ray, const RaySegment &segment, const WavelengthSamples &wls, LightPathSampler &pathSampler, bool* singleWavelength) const override;
        
        // END: MediumObject's methods
        // ----------------------------------------------------------------
    };
    
    
    
    class SLR_API MediumObjectAggregate : public MediumObject {
//        Accelerator* m_accelerator;
        BoundingBox3D m_bounds;
        std::vector<const MediumObject*> m_objLists;
        const MediumObject** m_lightList;
        std::map<uint32_t, uint32_t> m_objToLightMap;
        uint32_t m_numLights;
        DiscreteDistribution1D* m_lightDist1D;
    public:
        MediumObjectAggregate(const std::vector<MediumObject*> &objs);
        ~MediumObjectAggregate();
        
        // ----------------------------------------------------------------
        // Object's methods
        
        BoundingBox3D bounds() const override;
        
        // END: Object's methods
        // ----------------------------------------------------------------
        
        // ----------------------------------------------------------------
        // MediumObject's methods
        
        bool isEmitting() const override;
        float importance() const override;
        void selectLight(float u, float time, VolumetricLight* light, float* prob) const override;
        
        bool contains(const Point3D &p, float time) const override {
            return m_bounds.contains(p);
        }
        bool intersectBoundary(const Ray &ray, const RaySegment &segment, float* distToBoundary, bool* enter) const override {
            return m_bounds.intersectBoundary(ray, segment, distToBoundary, enter);
        }
        bool interact(const Ray &ray, const RaySegment &segment, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                      MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const override;
        SampledSpectrum evaluateTransmittance(const Ray &ray, const RaySegment &segment, const WavelengthSamples &wls, LightPathSampler &pathSampler, bool* singleWavelength) const override;
        
        // END: MediumObject's methods
        // ----------------------------------------------------------------
    };
}

#endif /* __SLR_medium_object__ */

//
//  MediumObject.h
//
//  Created by 渡部 心 on 2016/09/19.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef __SLR_MediumObject_h__
#define __SLR_MediumObject_h__

#include "../defines.h"
#include "../references.h"
#include "Object.h"

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
        mutable FixedStack<const MediumObject*> m_hierarchy;
    public:
        VolumetricLight() { }
        VolumetricLight(const FixedStack<const MediumObject*> &hierarchy) : m_hierarchy(hierarchy) { }
        VolumetricLight(const MediumInteraction &mi);
        
        void push(const MediumObject* obj) const { m_hierarchy.push(obj); }
        ScopedPop<const MediumObject*> scopedPop() const { return ScopedPop<const MediumObject*>(m_hierarchy); }
        const MediumObject* top() const { return m_hierarchy.top(); }
        
        SampledSpectrum sample(const LightPosQuery &query, const VolumetricLightPosSample &smp, VolumetricLightPosQueryResult* result) const;
        Ray sampleRay(const LightPosQuery &lightPosQuery, const VolumetricLightPosSample &lightPosSample,
                      const EDFQuery &edfQuery, const EDFSample &edfSample,
                      ArenaAllocator &mem,
                      VolumetricLightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                      EDFQueryResult* edfResult, SampledSpectrum* Le1) const;
        
        SampledSpectrum sample(const LightPosQuery &query, LightPathSampler &pathSampler, ArenaAllocator &mem, LightPosQueryResult** lpResult) const override;
    };
    
    
    
    class SLR_API MediumObject : public Object {
    public:
        virtual bool isEmitting() const = 0;
        virtual float importance() const = 0;
        virtual void selectLight(float u, VolumetricLight* light, float* prob) const = 0;
        virtual float evaluateProb(const VolumetricLight &light) const = 0;
        
        virtual SampledSpectrum sample(const VolumetricLight &light,
                                       const LightPosQuery &query, const VolumetricLightPosSample &smp, VolumetricLightPosQueryResult* result) const {
            SLRAssert(light.top() == this, "Object stored in Intersection does not match this object.");
            ScopedPop<const MediumObject*> sp = light.scopedPop();
            return light.top()->sample(light, query, smp, result);
        }
        virtual Ray sampleRay(const VolumetricLight &light,
                              const LightPosQuery &lightPosQuery, const VolumetricLightPosSample &lightPosSample,
                              const EDFQuery &edfQuery, const EDFSample &edfSample,
                              ArenaAllocator &mem,
                              VolumetricLightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                              EDFQueryResult* edfResult, SampledSpectrum* Le1) const {
            SLRAssert(light.top() == this, "Object stored in Intersection does not match this object.");
            ScopedPop<const MediumObject*> sp = light.scopedPop();
            return light.top()->sampleRay(light, lightPosQuery, lightPosSample, edfQuery, edfSample, mem, lightPosResult, Le0, edf, edfResult, Le1);
        }
        
        virtual bool getMediumContaining(const Point3D &p, const SingleMediumObject** curMedium) const = 0;
        virtual bool intersectBoundary(const Ray &ray, float* distToBoundary, bool* enter, const SingleMediumObject** boundaryMedium) const = 0;
        // The term "majorant" comes from the paper of Residual Ratio Tracking.
        virtual float majorantExtinctionCoefficient() const = 0;
        virtual SampledSpectrum extinctionCoefficient(const Point3D &p, const WavelengthSamples &wls) const = 0;
        virtual bool interact(Ray &ray, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                              MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const = 0;
        virtual SampledSpectrum evaluateTransmittance(Ray &ray, const WavelengthSamples &wls, LightPathSampler &pathSampler, bool* singleWavelength) const = 0;
        virtual void getMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const = 0;
        
        bool contains(const Light* light) const override {
            return ((const VolumetricLight*)light)->top() == this;
        }
        float evaluateProbability(const Light* light) const override {
            return evaluateProb(*(const VolumetricLight*)light);
        }
    };
    
    
    
    class SLR_API SingleMediumObject : public MediumObject {
        const Medium* m_medium;
        const MediumMaterial* m_material;
    public:
        SingleMediumObject(const Medium* medium, const MediumMaterial* material) :
        m_medium(medium), m_material(material) {}
        
        //----------------------------------------------------------------
        // Object's methods
        BoundingBox3D bounds() const override { return m_medium->bounds(); }
        //----------------------------------------------------------------
        
        //----------------------------------------------------------------
        // MediumObject's methods
        bool isEmitting() const override;
        float importance() const override;
        void selectLight(float u, VolumetricLight* light, float* prob) const override;
        float evaluateProb(const VolumetricLight &light) const override;
        
        SampledSpectrum sample(const VolumetricLight &light,
                               const LightPosQuery &query, const VolumetricLightPosSample &smp, VolumetricLightPosQueryResult* result) const override;
        Ray sampleRay(const VolumetricLight &light,
                      const LightPosQuery &lightPosQuery, const VolumetricLightPosSample &lightPosSample,
                      const EDFQuery &edfQuery, const EDFSample &edfSample,
                      ArenaAllocator &mem,
                      VolumetricLightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                      EDFQueryResult* edfResult, SampledSpectrum* Le1) const override;

        bool getMediumContaining(const Point3D &p, const SingleMediumObject** curMedium) const override;
        bool intersectBoundary(const Ray &ray, float* distToBoundary, bool* enter, const SingleMediumObject** boundaryMedium) const override;
        float majorantExtinctionCoefficient() const override;
        SampledSpectrum extinctionCoefficient(const Point3D &p, const WavelengthSamples &wls) const override;
        bool interact(Ray &ray, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                      MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const override;
        SampledSpectrum evaluateTransmittance(Ray &ray, const WavelengthSamples &wls, LightPathSampler &pathSampler, bool* singleWavelength) const override;
        void getMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const override;
        //----------------------------------------------------------------
        
        virtual AbstractBDF* createAbstractBDF(const MediumPoint &medPt, const WavelengthSamples &wls, ArenaAllocator &mem) const;
    };
    
    
    
    class SLR_API MediumObjectAggregate : public MediumObject {
//        Accelerator* m_accelerator;
        BoundingBox3D m_bounds;
        std::vector<const MediumObject*> m_objLists;
        const MediumObject** m_lightList;
        RegularConstantDiscrete1D* m_lightDist1D;
        std::map<const MediumObject*, uint32_t> m_revMap;
    public:
        MediumObjectAggregate(const std::vector<MediumObject*> &objs);
        ~MediumObjectAggregate();
        
        //----------------------------------------------------------------
        // Object's methods
        BoundingBox3D bounds() const override;
        //----------------------------------------------------------------
        
        //----------------------------------------------------------------
        // MediumObject's methods
        bool isEmitting() const override;
        float importance() const override;
        void selectLight(float u, VolumetricLight* light, float* prob) const override;
        float evaluateProb(const VolumetricLight &light) const override;
        
        bool getMediumContaining(const Point3D &p, const SingleMediumObject** curMedium) const override;
        bool intersectBoundary(const Ray &ray, float* distToBoundary, bool* enter, const SingleMediumObject** boundaryMedium) const override;
        float majorantExtinctionCoefficient() const override;
        SampledSpectrum extinctionCoefficient(const Point3D &p, const WavelengthSamples &wls) const override;
        bool interact(Ray &ray, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                      MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const override;
        SampledSpectrum evaluateTransmittance(Ray &ray, const WavelengthSamples &wls, LightPathSampler &pathSampler, bool* singleWavelength) const override;
        void getMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const override {
            SLRAssert(mi.top() == this, "Object stored in MediumInteraction does not match this object.");
            ScopedPop<const MediumObject*> sp = mi.scopedPop();
            mi.getMediumPoint(medPt);
        }
        //----------------------------------------------------------------
    };
}

#endif /* MediumObject_h */

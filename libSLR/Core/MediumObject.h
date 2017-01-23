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
        const SingleMediumObject* m_obj;
    public:
        VolumetricLight() { }
        
        void setObject(const SingleMediumObject* obj) { m_obj = obj; }
        
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
        virtual void selectLight(float u, float time, VolumetricLight* light, float* prob) const = 0;
        
        virtual SampledSpectrum sample(const StaticTransform &transform,
                                       const LightPosQuery &query, const VolumetricLightPosSample &smp, VolumetricLightPosQueryResult* result) const {
            SLRAssert_ShouldNotBeCalled();
            return SampledSpectrum::Zero;
        }
        virtual Ray sampleRay(const StaticTransform &transform,
                              const LightPosQuery &lightPosQuery, const VolumetricLightPosSample &lightPosSample,
                              const EDFQuery &edfQuery, const EDFSample &edfSample,
                              ArenaAllocator &mem,
                              VolumetricLightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                              EDFQueryResult* edfResult, SampledSpectrum* Le1) const {
            SLRAssert_ShouldNotBeCalled();
            return Ray();
        }
        
        virtual bool getMediumContaining(const Point3D &p, const SingleMediumObject** curMedium) const = 0;
        virtual bool intersectBoundary(const Ray &ray, float* distToBoundary, bool* enter, const SingleMediumObject** boundaryMedium) const = 0;
        // The term "majorant" comes from the paper of Residual Ratio Tracking.
        virtual float majorantExtinctionCoefficient() const = 0;
        virtual SampledSpectrum extinctionCoefficient(const Point3D &p, const WavelengthSamples &wls) const = 0;
        virtual bool interact(Ray &ray, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                              MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const = 0;
        virtual SampledSpectrum evaluateTransmittance(Ray &ray, const WavelengthSamples &wls, LightPathSampler &pathSampler, bool* singleWavelength) const = 0;
        virtual void getMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const {
            SLRAssert_ShouldNotBeCalled();
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
        void selectLight(float u, float time, VolumetricLight* light, float* prob) const override;
        
        SampledSpectrum sample(const StaticTransform &transform,
                               const LightPosQuery &query, const VolumetricLightPosSample &smp, VolumetricLightPosQueryResult* result) const override;
        Ray sampleRay(const StaticTransform &transform,
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
        std::map<uint32_t, uint32_t> m_objToLightMap;
        uint32_t m_numLights;
        RegularConstantDiscrete1D* m_lightDist1D;
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
        void selectLight(float u, float time, VolumetricLight* light, float* prob) const override;
        
        bool getMediumContaining(const Point3D &p, const SingleMediumObject** curMedium) const override;
        bool intersectBoundary(const Ray &ray, float* distToBoundary, bool* enter, const SingleMediumObject** boundaryMedium) const override;
        float majorantExtinctionCoefficient() const override;
        SampledSpectrum extinctionCoefficient(const Point3D &p, const WavelengthSamples &wls) const override;
        bool interact(Ray &ray, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                      MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const override;
        SampledSpectrum evaluateTransmittance(Ray &ray, const WavelengthSamples &wls, LightPathSampler &pathSampler, bool* singleWavelength) const override;
        //----------------------------------------------------------------
    };
}

#endif /* MediumObject_h */

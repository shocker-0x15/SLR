//
//  medium_object.cpp
//
//  Created by 渡部 心 on 2016/09/19.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "medium_object.h"

#include "../MemoryAllocators/ArenaAllocator.h"
#include "distributions.h"
#include "surface_object.h"
#include "light_path_sampler.h"
#include "textures.h"
#include "medium_material.h"
#include "../Accelerator/StandardBVH.h"
#include "../Accelerator/SBVH.h"
#include "../Accelerator/QBVH.h"
#include "../Scene/Scene.h"

namespace SLR {
    SampledSpectrum VolumetricLight::sample(const LightPosQuery &query, const VolumetricLightPosSample &smp, VolumetricLightPosQueryResult* result) const {
        return m_obj->sample(m_appliedTransform, query, smp, result);
    }
    
    Ray VolumetricLight::sampleRay(const LightPosQuery &lightPosQuery, const VolumetricLightPosSample &lightPosSample,
                                   const EDFQuery &edfQuery, const EDFSample &edfSample,
                                   ArenaAllocator &mem,
                                   VolumetricLightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                                   EDFQueryResult* edfResult, SampledSpectrum* Le1) const {
        return m_obj->sampleRay(m_appliedTransform, lightPosQuery, lightPosSample, edfQuery, edfSample, mem, lightPosResult, Le0, edf, edfResult, Le1);
    }
    
    SampledSpectrum VolumetricLight::sample(const LightPosQuery &query, LightPathSampler &pathSampler, ArenaAllocator &mem, LightPosQueryResult** lpResult) const {
        VolumetricLightPosQueryResult* result = mem.create<VolumetricLightPosQueryResult>();
        SampledSpectrum ret = sample(query, pathSampler.getVolumetricLightPosSample(), result);
        *lpResult = result;
        return ret;
    }
    
    Ray VolumetricLight::sampleRay(const LightPosQuery &lightPosQuery, LightPathSampler &pathSampler, const EDFQuery &edfQuery, ArenaAllocator &mem,
                                   LightPosQueryResult** lightPosResult, SampledSpectrum* Le0, EDF** edf,
                                   EDFQueryResult* edfResult, SampledSpectrum* Le1) const {
        *lightPosResult = mem.create<VolumetricLightPosQueryResult>();
        return sampleRay(lightPosQuery, pathSampler.getVolumetricLightPosSample(), edfQuery, pathSampler.getEDFSample(), mem,
                         (VolumetricLightPosQueryResult*)*lightPosResult, Le0, edf, edfResult, Le1);
    }
    
    
    
    bool SingleMediumObject::isEmitting() const {
        return m_material->isEmitting();
    }
    
    float SingleMediumObject::importance() const {
        return 1.0f;
    }
    
    void SingleMediumObject::selectLight(float u, float time, VolumetricLight* light, float* prob) const {
        light->setObject(this);
        *prob = 1.0f;
    }
    
    SampledSpectrum SingleMediumObject::sample(const StaticTransform &transform,
                                               const LightPosQuery &query, const VolumetricLightPosSample &smp, VolumetricLightPosQueryResult* result) const {
        SLRAssert_NotImplemented();
        return SampledSpectrum::Zero;
    }
    
    Ray SingleMediumObject::sampleRay(const StaticTransform &transform,
                                      const LightPosQuery &lightPosQuery, const VolumetricLightPosSample &lightPosSample,
                                      const EDFQuery &edfQuery, const EDFSample &edfSample,
                                      ArenaAllocator &mem,
                                      VolumetricLightPosQueryResult* lightPosResult, SampledSpectrum* Le0, SLR::EDF** edf,
                                      EDFQueryResult* edfResult, SampledSpectrum* Le1) const {
        SLRAssert_NotImplemented();
        return Ray();
    }
    
    bool SingleMediumObject::interact(Ray &ray, float distanceLimit, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                                      MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const {
        if (!m_medium->interact(ray, distanceLimit, wls, pathSampler, mi, medThroughput, singleWavelength))
            return false;
        mi->setObject(this);
        mi->setLightProb(isEmitting() ? 1.0f : 0.0f);
        
        return true;
    }
    
    SampledSpectrum SingleMediumObject::evaluateTransmittance(Ray &ray, float distanceLimit, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                                                              bool* singleWavelength) const {
        return m_medium->evaluateTransmittance(ray, distanceLimit, wls, pathSampler, singleWavelength);
    }
    
    void SingleMediumObject::getMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const {
        m_medium->getMediumPoint(mi, medPt);
        medPt->setObject(this);
        medPt->applyTransform(mi.getAppliedTransform());
    }
    
    SampledSpectrum SingleMediumObject::extinctionCoefficient(const MediumPoint &medPt, const WavelengthSamples &wls) const {
        float u, v, w;
        medPt.getMediumParameter(&u, &v, &w);
        return m_medium->getExtinctionCoefficient(Point3D(u, v, w), wls);
    }
    
    AbstractBDF* SingleMediumObject::createAbstractBDF(const MediumPoint &medPt, const WavelengthSamples &wls, ArenaAllocator &mem) const {
        PhaseFunction* pf = m_material->getPhaseFunction(medPt, wls, mem);
        float u, v, w;
        medPt.getMediumParameter(&u, &v, &w);
        SampledSpectrum albedo = m_medium->getAlbedo(Point3D(u, v, w), wls);
        return mem.create<VolumetricBSDF>(albedo, pf);
    }
    
    
    
    BoundingBox3D TransformedMediumObject::bounds() const {
        return m_transform->motionBounds(m_medObj->bounds());
    }
    
    bool TransformedMediumObject::isEmitting() const {
        return m_medObj->isEmitting();
    }
    
    float TransformedMediumObject::importance() const {
        return m_medObj->importance();
    }
    
    void TransformedMediumObject::selectLight(float u, float time, VolumetricLight* light, float* prob) const {
        m_medObj->selectLight(u, time, light, prob);
        StaticTransform tf;
        m_transform->sample(time, &tf);
        light->applyTransformFromLeft(tf);
    }
    
    bool TransformedMediumObject::contains(const Point3D &p, float time) const {
        StaticTransform tf;
        m_transform->sample(time, &tf);
        Point3D localP = invert(tf) * p;
        return m_medObj->contains(localP, time);
    }
    
    bool TransformedMediumObject::intersectBoundary(const Ray &ray, float* distToBoundary, bool* enter) const {
        StaticTransform tf;
        m_transform->sample(ray.time, &tf);
        Ray localRay = invert(tf) * ray;
        return m_medObj->intersectBoundary(localRay, distToBoundary, enter);
    }
    
    bool TransformedMediumObject::interact(Ray &ray, float distanceLimit, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                                           MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const {
        StaticTransform tf;
        m_transform->sample(ray.time, &tf);
        Ray localRay = invert(tf) * ray;
        bool hit = m_medObj->interact(localRay, distanceLimit, wls, pathSampler, mi, medThroughput, singleWavelength);
        ray.distMin = distanceLimit;
        if (hit)
            mi->applyTransformFromLeft(tf);
        return hit;
    }
    
    SampledSpectrum TransformedMediumObject::evaluateTransmittance(Ray &ray, float distanceLimit, const WavelengthSamples &wls, LightPathSampler &pathSampler, bool* singleWavelength) const {
        StaticTransform tf;
        m_transform->sample(ray.time, &tf);
        Ray localRay = invert(tf) * ray;
        SampledSpectrum ret = m_medObj->evaluateTransmittance(localRay, distanceLimit, wls, pathSampler, singleWavelength);
        ray.distMin = distanceLimit;
        return ret;
    }
    
    
    
    BoundingBox3D EnclosedMediumObject::bounds() const {
        return intersection(m_boundary->bounds(), m_medToSurfTF * m_medObj->bounds());
    }
    
    bool EnclosedMediumObject::contains(const Point3D &p, float time) const {
        if (!m_medObj->contains(p, time))
            return false;
        
        Point3D pInSurf = m_medToSurfTF * p;
        return m_boundary->contains(pInSurf, time);
    }
    
    bool EnclosedMediumObject::intersectBoundary(const Ray &ray, float* distToBoundary, bool* enter) const {
        SurfaceInteraction si;
        Ray r = m_medToSurfTF * ray;
        bool hit = m_boundary->intersect(r, &si);
        if (!hit)
            return false;
        float distToSurf = si.getDistance();
        bool enterToSurf = dot(si.getGeometricNormal(), r.dir) < 0.0f;
        
        float distToRawMed;
        bool enterToRawMed;
        hit = m_medObj->intersectBoundary(ray, &distToRawMed, &enterToRawMed);
        if (!hit)
            return false;
        
        if (enterToSurf && enterToRawMed) {
            *distToBoundary = std::max(distToSurf, distToRawMed);
            *enter = true;
        }
        else if (enterToSurf && !enterToRawMed) {
            *distToBoundary = std::min(distToSurf, distToRawMed);
            *enter = true;
        }
        else if (!enterToSurf && enterToRawMed) {
            *distToBoundary = std::min(distToSurf, distToRawMed);
            *enter = true;
        }
        else {
            *distToBoundary = std::min(distToSurf, distToRawMed);
            *enter = false;
        }
        
        return true;
    }
    
    bool EnclosedMediumObject::interact(Ray &ray, float distanceLimit, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                                        MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const {
        return m_medObj->interact(ray, distanceLimit, wls, pathSampler, mi, medThroughput, singleWavelength);
    }
    
    SampledSpectrum EnclosedMediumObject::evaluateTransmittance(Ray &ray, float distanceLimit, const WavelengthSamples &wls, LightPathSampler &pathSampler, bool* singleWavelength) const {
        return m_medObj->evaluateTransmittance(ray, distanceLimit, wls, pathSampler, singleWavelength);
    }
    
    
    
    MediumObjectAggregate::MediumObjectAggregate(const std::vector<MediumObject*> &objs) {
        BoundingBox3D bbox;
        for (int i = 0; i < objs.size(); ++i)
            bbox.unify(objs[i]->bounds());
        m_bounds = bbox;
        
        for (int i = 0; i < objs.size(); ++i)
            m_objLists.push_back(objs[i]);
        
        std::vector<uint32_t> lightIndices;
        std::vector<float> lightImportances;
        for (int i = 0; i < objs.size(); ++i) {
            const MediumObject* obj = objs[i];
            if (obj->isEmitting()) {
                lightIndices.push_back(i);
                lightImportances.push_back(obj->importance());
            }
        }
        
        m_numLights = (uint32_t)lightImportances.size();
        m_lightList = new const MediumObject*[m_numLights];
        m_lightDist1D = new DiscreteDistribution1D(lightImportances);
        
        for (int i = 0; i < m_numLights; ++i) {
            uint32_t objIdx = lightIndices[i];
            const MediumObject* light = objs[objIdx];
            m_lightList[i] = light;
            m_objToLightMap[objIdx] = i;
        }
    }
    
    MediumObjectAggregate::~MediumObjectAggregate() {
        delete m_lightDist1D;
        delete[] m_lightList;
    }
    
    BoundingBox3D MediumObjectAggregate::bounds() const {
        return m_bounds;
    }
    
    bool MediumObjectAggregate::isEmitting() const {
        return m_numLights > 0;
    }
    
    float MediumObjectAggregate::importance() const {
        return m_lightDist1D->integral();
    }
    
    void MediumObjectAggregate::selectLight(float u, float time, VolumetricLight* light, float* prob) const {
        uint32_t lIdx = m_lightDist1D->sample(u, prob, &u);
        const MediumObject* obj = m_lightList[lIdx];
        float cProb;
        obj->selectLight(u, time, light, &cProb);
        *prob *= cProb;
    }
    
    bool MediumObjectAggregate::interact(Ray &ray, float distanceLimit, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                                         MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const {        
        *medThroughput = SampledSpectrum::One;
        *singleWavelength = false;
        
        Point3D currentPoint = ray.org + ray.distMin * ray.dir;
        // The following process is logically the same as:
        // curMedia = m_accelerator->queryCurrentMedia(currentPoint);
        const MediumObject* curMedium = nullptr;
        for (int i = 0; i < m_objLists.size(); ++i) {
            if (m_objLists[i]->contains(currentPoint, ray.time))
                curMedium = m_objLists[i];
        }
        
        while (true) {
            const MediumObject* nextMedium = nullptr;
            float distToNextBoundary = INFINITY;
            // The following process is logically the same as:
            // nextMedium = m_accelerator->queryNextMedia(ray, &distToNextBoundary);
            uint32_t objIdx = -1;
            for (int i = 0; i < m_objLists.size(); ++i) {
                float distToBoundary = INFINITY;
                bool enter;
                if (m_objLists[i]->intersectBoundary(ray, &distToBoundary, &enter)) {
                    if (distToBoundary < distToNextBoundary) {
                        distToNextBoundary = distToBoundary;
                        objIdx = enter ? i : -1;
                        nextMedium = enter ? m_objLists[i] : nullptr;
                    }
                }
            }
            distToNextBoundary = std::min(distToNextBoundary, distanceLimit);
            if (curMedium && std::isinf(distToNextBoundary))
                return false;
            
            bool hit = false;
            if (curMedium) {
                SampledSpectrum curMedThroughput;
                bool curSingleWavelength;
                hit = curMedium->interact(ray, distToNextBoundary, wls, pathSampler, mi, &curMedThroughput, &curSingleWavelength);
                *medThroughput *= curMedThroughput;
                *singleWavelength |= curSingleWavelength;
            }
            if (hit) {
                if (m_objToLightMap.count(objIdx) > 0) {
                    uint32_t lightIdx = m_objToLightMap.at(objIdx);
                    mi->setLightProb(m_lightDist1D->evaluatePMF(lightIdx) * mi->getLightProb());
                }
                return true;
            }
            
            if (distToNextBoundary == distanceLimit)
                return false;
            
            ray.distMin = distToNextBoundary * (1.0f + Ray::Epsilon);
            curMedium = nextMedium;
        }
        
        SLRAssert(false, "This code path should never be executed.");
        return true;
    }
    
    SampledSpectrum MediumObjectAggregate::evaluateTransmittance(Ray &ray, float distanceLimit, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                                                                 bool* singleWavelength) const {
        SampledSpectrum transmittance = SampledSpectrum::One;
        *singleWavelength = false;
        
        Point3D currentPoint = ray.org + ray.distMin * ray.dir;
        // The following process is logically the same as:
        // curMedia = m_accelerator->queryCurrentMedia(currentPoint);
        const MediumObject* curMedium = nullptr;
        for (int i = 0; i < m_objLists.size(); ++i) {
            if (m_objLists[i]->contains(currentPoint, ray.time))
                curMedium = m_objLists[i];
        }
        
        while (true) {
            const MediumObject* nextMedium = nullptr;
            float distToNextBoundary = INFINITY;
            // The following process is logically the same as:
            // nextMedium = m_accelerator->queryNextMedia(ray, &distToNextBoundary);
            for (int i = 0; i < m_objLists.size(); ++i) {
                float distToBoundary = INFINITY;
                bool enter;
                if (m_objLists[i]->intersectBoundary(ray, &distToBoundary, &enter)) {
                    if (distToBoundary < distToNextBoundary) {
                        distToNextBoundary = distToBoundary;
                        nextMedium = enter ? m_objLists[i] : nullptr;
                    }
                }
            }
            distToNextBoundary = std::min(distToNextBoundary, distanceLimit);
            if (curMedium && std::isinf(distToNextBoundary))
                return false;
            
            if (curMedium) {
                SampledSpectrum curTransmittance;
                bool curSingleWavelength;
                curTransmittance = curMedium->evaluateTransmittance(ray, distToNextBoundary, wls, pathSampler, &curSingleWavelength);
                transmittance *= curTransmittance;
                *singleWavelength |= curSingleWavelength;
            }
            
            if (distToNextBoundary == distanceLimit)
                break;
            
            ray.distMin = distToNextBoundary * (1.0f + Ray::Epsilon);
            curMedium = nextMedium;
        }
        
        return transmittance;
    }
}

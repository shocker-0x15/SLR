//
//  MediumObject.cpp
//
//  Created by 渡部 心 on 2016/09/19.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "MediumObject.h"
#include "distributions.h"
#include "textures.h"
#include "../Memory/ArenaAllocator.h"
#include "../Core/medium_material.h"
#include "../Core/light_path_samplers.h"
#include "../Accelerator/StandardBVH.h"
#include "../Accelerator/SBVH.h"
#include "../Accelerator/QBVH.h"
#include "../Core/light_path_samplers.h"
#include "../Scene/Scene.h"

namespace SLR {
    VolumetricLight::VolumetricLight(const MediumInteraction &mi) : m_hierarchy(mi.getHierarchy()) {
        SLRAssert(m_hierarchy.top()->isEmitting(), "This is not a light.");
    }
    
    SampledSpectrum VolumetricLight::sample(const LightPosQuery &query, const VolumetricLightPosSample &smp, VolumetricLightPosQueryResult* result) const {
        return m_hierarchy.top()->sample(*this, query, smp, result);
    }
    
    Ray VolumetricLight::sampleRay(const LightPosQuery &lightPosQuery, const VolumetricLightPosSample &lightPosSample,
                                   const EDFQuery &edfQuery, const EDFSample &edfSample,
                                   ArenaAllocator &mem,
                                   VolumetricLightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                                   EDFQueryResult* edfResult, SampledSpectrum* Le1) const {
        return m_hierarchy.top()->sampleRay(*this, lightPosQuery, lightPosSample, edfQuery, edfSample, mem, lightPosResult, Le0, edf, edfResult, Le1);
    }
    
    SampledSpectrum VolumetricLight::sample(const LightPosQuery &query, LightPathSampler &pathSampler, ArenaAllocator &mem, LightPosQueryResult** lpResult) const {
        VolumetricLightPosQueryResult* result = mem.create<VolumetricLightPosQueryResult>();
        SampledSpectrum ret = sample(query, pathSampler.getVolumetricLightPosSample(), result);
        *lpResult = result;
        return ret;
    }
    
    
    
    bool SingleMediumObject::isEmitting() const {
        return m_material->isEmitting();
    }
    
    float SingleMediumObject::importance() const {
        return 1.0f;
    }
    
    void SingleMediumObject::selectLight(float u, VolumetricLight* light, float* prob) const {
        light->push(this);
        *prob = 1.0f;
    }
    
    float SingleMediumObject::evaluateProb(const VolumetricLight &light) const {
        SLRAssert(light.top() == this, "Object stored in Light does not match this object.");
        return 1.0f;
    }
    
    SampledSpectrum SingleMediumObject::sample(const VolumetricLight &light,
                                               const LightPosQuery &query, const VolumetricLightPosSample &smp, VolumetricLightPosQueryResult* result) const {
        SLRAssert_NotImplemented();
        return SampledSpectrum::Zero;
    }
    
    Ray SingleMediumObject::sampleRay(const VolumetricLight &light,
                                      const LightPosQuery &lightPosQuery, const VolumetricLightPosSample &lightPosSample,
                                      const EDFQuery &edfQuery, const EDFSample &edfSample,
                                      ArenaAllocator &mem,
                                      VolumetricLightPosQueryResult* lightPosResult, SampledSpectrum* Le0, SLR::EDF** edf,
                                      EDFQueryResult* edfResult, SampledSpectrum* Le1) const {
        SLRAssert_NotImplemented();
        return Ray();
    }
    
    bool SingleMediumObject::getMediumContaining(const Point3D &p, const SingleMediumObject** curMedium) const {
        if (m_medium->contains(p)) {
            *curMedium = this;
            return true;
        }
        return false;
    }
    
    bool SingleMediumObject::intersectBoundary(const Ray &ray, float* distToBoundary, bool* enter, const SingleMediumObject** boundaryMedium) const {
        if (m_medium->intersectBoundary(ray, distToBoundary, enter)) {
            *boundaryMedium = this;
            return true;
        }
        return false;
    }
    
    float SingleMediumObject::majorantExtinctionCoefficient() const {
        return m_medium->majorantExtinctionCoefficient();
    }
    
    SampledSpectrum SingleMediumObject::extinctionCoefficient(const Point3D &p, const WavelengthSamples &wls) const {
        return m_medium->extinctionCoefficient(p, wls);
    }
    
    bool SingleMediumObject::interact(Ray &ray, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                                      MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const {
        SLRAssert_NotImplemented();
        return true;
    }
    
    SampledSpectrum SingleMediumObject::evaluateTransmittance(Ray &ray, const WavelengthSamples &wls, LightPathSampler &pathSampler, bool* singleWavelength) const {
        SLRAssert_NotImplemented();
        return SampledSpectrum::Zero;
    }
    
    void SingleMediumObject::getMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const {
        SLRAssert_NotImplemented();
    }
    
    AbstractBDF* SingleMediumObject::createAbstractBDF(const MediumPoint &medPt, const WavelengthSamples &wls, ArenaAllocator &mem) const {
        PhaseFunction* pf = m_material->getPhaseFunction(medPt, wls, mem);
        SampledSpectrum sigma_s, sigma_e;
        m_medium->queryCoefficients(medPt.getPosition(), wls, &sigma_s, &sigma_e);
        return mem.create<VolumetricBSDF>(sigma_s / sigma_e, pf);
    }
    
    
    
    MediumObjectAggregate::MediumObjectAggregate(const std::vector<MediumObject*> &objs) {
        BoundingBox3D bbox;
        for (int i = 0; i < objs.size(); ++i)
            bbox.unify(objs[i]->bounds());
        m_bounds = bbox;
        
        for (int i = 0; i < objs.size(); ++i)
            m_objLists.push_back(objs[i]);
        
        std::vector<const MediumObject*> lights;
        std::vector<float> lightImportances;
        for (int i = 0; i < objs.size(); ++i) {
            const MediumObject* obj = objs[i];
            if (obj->isEmitting()) {
                lights.push_back(obj);
                lightImportances.push_back(obj->importance());
            }
        }
        
        m_lightList = new const MediumObject*[lightImportances.size()];
        m_lightDist1D = new RegularConstantDiscrete1D(lightImportances);
        
        for (int i = 0; i < lights.size(); ++i) {
            const MediumObject* light = lights[i];
            m_lightList[i] = light;
            m_revMap[light] = i;
        }
    }
    
    MediumObjectAggregate::~MediumObjectAggregate() {
        delete m_lightList;
        delete m_lightDist1D;
    }
    
    BoundingBox3D MediumObjectAggregate::bounds() const {
        return m_bounds;
    }
    
    bool MediumObjectAggregate::isEmitting() const {
        return m_revMap.size() > 0;
    }
    
    float MediumObjectAggregate::importance() const {
        return m_lightDist1D->integral();
    }
    
    void MediumObjectAggregate::selectLight(float u, VolumetricLight* light, float* prob) const {
        uint32_t lIdx = m_lightDist1D->sample(u, prob, &u);
        const MediumObject* obj = m_lightList[lIdx];
        float cProb;
        obj->selectLight(u, light, &cProb);
        *prob *= cProb;
        light->push(this);
    }
    
    float MediumObjectAggregate::evaluateProb(const VolumetricLight &light) const {
        SLRAssert(light.top() == this, "Object stored in Intersection does not match this object.");
        
        ScopedPop<const MediumObject*> sp = light.scopedPop();
        SLRAssert(m_revMap.count(light.top()) > 0, "Specified light is not found.");
        float prob = m_lightDist1D->evaluatePMF(m_revMap.at(light.top()));
        return prob * light.top()->evaluateProb(light);
    }
    
    bool MediumObjectAggregate::getMediumContaining(const Point3D &p, const SingleMediumObject** curMedium) const {
        for (int i = 0; i < m_objLists.size(); ++i) {
            if (m_objLists[i]->getMediumContaining(p, curMedium))
                return true;
        }
        return false;
    }
    
    bool MediumObjectAggregate::intersectBoundary(const Ray &ray, float* distToBoundary, bool* enter, const SingleMediumObject** boundaryMedium) const {
        *distToBoundary = INFINITY;
        for (int i = 0; i < m_objLists.size(); ++i) {
            float dist = INFINITY;
            bool isEntering;
            const SingleMediumObject* medium;
            if (m_objLists[i]->intersectBoundary(ray, &dist, &isEntering, &medium)) {
                if (dist < *distToBoundary) {
                    *distToBoundary = dist;
                    *enter = isEntering;
                    *boundaryMedium = medium;
                }
            }
        }
        return !std::isinf(*distToBoundary);
    }
    
    float MediumObjectAggregate::majorantExtinctionCoefficient() const {
        float maxMajorantExtinctionCoefficient = 0;
        for (int i = 0; i < m_objLists.size(); ++i)
            maxMajorantExtinctionCoefficient = std::max(maxMajorantExtinctionCoefficient, m_objLists[i]->majorantExtinctionCoefficient());
        return maxMajorantExtinctionCoefficient;
    }
    
    SampledSpectrum MediumObjectAggregate::extinctionCoefficient(const Point3D &p, const WavelengthSamples &wls) const {
        const SingleMediumObject* medium;
        if (getMediumContaining(p, &medium))
            return medium->extinctionCoefficient(p, wls);
        return SampledSpectrum::Zero;
    }
    
    bool MediumObjectAggregate::interact(Ray &ray, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                                         MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const {
        FreePathSampler sampler = pathSampler.getFreePathSampler();
        FloatSum sampledDistance = ray.distMin;
        
        // Current implementation does not allow overlapped media.
        // If allowing overlapping, it might require to get medium list containing queryPoint.
        Point3D queryPoint = ray.org + ray.distMin * ray.dir;
        const SingleMediumObject* curMedium = nullptr;
        for (int i = 0; i < m_objLists.size(); ++i) {
            m_objLists[i]->getMediumContaining(queryPoint, &curMedium);
            if (curMedium)
                break;
        }
        
        while (true) {
            // find the next boundary distance and medium
            ray.distMin = sampledDistance;
            bool enter;
            float distToNextBoundary = ray.distMax;
            const SingleMediumObject* boundaryMedium = nullptr;
            if (curMedium) {
                ray.distMin *= 1 + Ray::Epsilon;
                float distToBoundary = INFINITY;
                curMedium->intersectBoundary(ray, &distToBoundary, &enter, &boundaryMedium);
            }
            else {
                for (int i = 0; i < m_objLists.size(); ++i) {
                    float distToBoundary = INFINITY;
                    if (m_objLists[i]->intersectBoundary(ray, &distToBoundary, &enter, &boundaryMedium)) {
                        // Because overlapped media are not allowed, exitting in this code path is invalid.
                        if (enter == false)
                            continue;
                        distToNextBoundary = std::min(distToNextBoundary, distToBoundary);
                    }
                }
            }
            
            float majorant = 0.0f;
            if (curMedium)
                majorant = curMedium->majorantExtinctionCoefficient();
            
            sampledDistance += -std::log(sampler.getSample()) / majorant;
            while (sampledDistance < distToNextBoundary) {
                queryPoint = ray.org + sampledDistance * ray.dir;
                
                SampledSpectrum extCoeff = curMedium->extinctionCoefficient(queryPoint, wls);
                float probRealCollision = extCoeff[wls.selectedLambda] / majorant;
                if (sampler.getSample() < probRealCollision) {
                    *mi = MediumInteraction(ray.time, sampledDistance, queryPoint);
                    *medThroughput = SampledSpectrum::Zero;
                    (*medThroughput)[wls.selectedLambda] = 1.0f / extCoeff[wls.selectedLambda];
                    *singleWavelength = true;
                    return true;
                }
                sampledDistance += -std::log(sampler.getSample()) / majorant;
            }
            
            sampledDistance = distToNextBoundary;
            curMedium = boundaryMedium;
            
            if (sampledDistance >= ray.distMax) {
                *medThroughput = SampledSpectrum::Zero;
                (*medThroughput)[wls.selectedLambda] = 1.0f;
                *singleWavelength = true;
                return false;
            }
        }
        
        SLRAssert(false, "This code path should never be executed.");
        return true;
    }
    
    SampledSpectrum MediumObjectAggregate::evaluateTransmittance(Ray &ray, const WavelengthSamples &wls, LightPathSampler &pathSampler, bool *singleWavelength) const {
        FreePathSampler sampler = pathSampler.getFreePathSampler();
        FloatSum sampledDistance = ray.distMin;
        
        // Current implementation does not allow overlapped media.
        // If allowing overlapping, it might require to get medium list containing queryPoint.
        Point3D queryPoint = ray.org + ray.distMin * ray.dir;
        const SingleMediumObject* curMedium = nullptr;
        for (int i = 0; i < m_objLists.size(); ++i) {
            m_objLists[i]->getMediumContaining(queryPoint, &curMedium);
            if (curMedium)
                break;
        }
        
        SampledSpectrum transmittance = SampledSpectrum::Zero;
        transmittance[wls.selectedLambda] = 1.0f;
        while (true) {
            // find the next boundary distance and medium
            ray.distMin = sampledDistance;
            bool enter;
            float distToNextBoundary = ray.distMax;
            const SingleMediumObject* boundaryMedium = nullptr;
            if (curMedium) {
                ray.distMin *= 1 + Ray::Epsilon;
                float distToBoundary = INFINITY;
                curMedium->intersectBoundary(ray, &distToBoundary, &enter, &boundaryMedium);
            }
            else {
                for (int i = 0; i < m_objLists.size(); ++i) {
                    float distToBoundary = INFINITY;
                    if (m_objLists[i]->intersectBoundary(ray, &distToBoundary, &enter, &boundaryMedium)) {
                        // Because overlapped media are not allowed, exitting in this code path is invalid.
                        if (enter == false)
                            continue;
                        distToNextBoundary = std::min(distToNextBoundary, distToBoundary);
                    }
                }
            }
            
            float majorant = 0.0f;
            if (curMedium)
                majorant = curMedium->majorantExtinctionCoefficient();
            
            sampledDistance += -std::log(sampler.getSample()) / majorant;
            while (sampledDistance < distToNextBoundary) {
                queryPoint = ray.org + sampledDistance * ray.dir;
                
                SampledSpectrum extCoeff = curMedium->extinctionCoefficient(queryPoint, wls);
                float probRealCollision = extCoeff[wls.selectedLambda] / majorant;
                transmittance[wls.selectedLambda] *= (1.0f - probRealCollision);
                sampledDistance += -std::log(sampler.getSample()) / majorant;
            }
            
            sampledDistance = distToNextBoundary;
            curMedium = boundaryMedium;
            
            if (sampledDistance >= ray.distMax) {
                *singleWavelength = true;
                return transmittance;
            }
        }
        
        SLRAssert(false, "This code path should never be executed.");
        return transmittance;
    }
}

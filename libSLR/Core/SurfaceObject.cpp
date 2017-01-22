//
//  SurfaceObject.cpp
//
//  Created by 渡部 心 on 2015/07/15.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "SurfaceObject.h"
#include "surface_material.h"
#include "Transform.h"
#include "distributions.h"
#include "textures.h"
#include "../Memory/ArenaAllocator.h"
#include "../Accelerator/StandardBVH.h"
#include "../Accelerator/SBVH.h"
#include "../Accelerator/QBVH.h"
#include "../Core/light_path_samplers.h"
#include "../Surface/InfiniteSphere.h"
#include "../SurfaceMaterials/IBLEmission.h"
#include "../Scene/Scene.h"
#include "../BSDFs/basic_BSDFs.h"

namespace SLR {
    SurfaceLight::SurfaceLight(const SurfaceInteraction &si) : m_hierarchy(si.getHierarchy()) {
        SLRAssert(m_hierarchy.top()->isEmitting(), "This is not a light.");
    }
    
    SampledSpectrum SurfaceLight::sample(const LightPosQuery &query, const SurfaceLightPosSample &smp, SurfaceLightPosQueryResult* result) const {
        return m_hierarchy.top()->sample(*this, query, smp, result);
    }
    
    Ray SurfaceLight::sampleRay(const LightPosQuery &lightPosQuery, const SurfaceLightPosSample &lightPosSample,
                                const EDFQuery &edfQuery, const EDFSample &edfSample,
                                ArenaAllocator &mem,
                                SurfaceLightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                                EDFQueryResult* edfResult, SampledSpectrum* Le1) const {
        return m_hierarchy.top()->sampleRay(*this, lightPosQuery, lightPosSample, edfQuery, edfSample, mem, lightPosResult, Le0, edf, edfResult, Le1);
    }
    
    SampledSpectrum SurfaceLight::sample(const LightPosQuery &query, LightPathSampler &pathSampler, ArenaAllocator &mem, LightPosQueryResult** lpResult) const {
        SurfaceLightPosQueryResult* result = mem.create<SurfaceLightPosQueryResult>();
        SampledSpectrum ret = sample(query, pathSampler.getSurfaceLightPosSample(), result);
        *lpResult = result;
        return ret;
    }
    
    
    
    bool SurfaceObject::intersect(Ray &ray, SurfacePoint *surfPt) const {
        SurfaceInteraction si;
        if (!intersect(ray, &si))
            return false;
        getSurfacePoint(si, surfPt);
        return true;
    }
    
    bool SurfaceObject::testVisibility(const SurfacePoint &shdP, const SurfacePoint &lightP, float time) const {
        SLRAssert(shdP.atInfinity() == false && lightP.atInfinity() == false, "Points must be in finite region.");
        float dist = distance(lightP.getPosition(), shdP.getPosition());
        Ray ray(shdP.getPosition(), (lightP.getPosition() - shdP.getPosition()) / dist, time, Ray::Epsilon, dist * (1 - Ray::Epsilon));
        SurfaceInteraction si;
        return !intersect(ray, &si);
    }
    
    
    
    bool SingleSurfaceObject::isEmitting() const {
        return m_material->isEmitting();
    }
    
    float SingleSurfaceObject::importance() const {
        return 1.0f;// TODO: consider a total power emitted from this object.
    }
    
    void SingleSurfaceObject::selectLight(float u, SurfaceLight* light, float* prob) const {
        light->push(this);
        *prob = 1.0f;
    }
    
    float SingleSurfaceObject::evaluateProb(const SurfaceLight &light) const {
        SLRAssert(light.top() == this, "Object stored in Light does not match this object.");
        return 1.0f;
    }
    
    SampledSpectrum SingleSurfaceObject::sample(const SurfaceLight &light,
                                                const LightPosQuery &query, const SurfaceLightPosSample &smp, SurfaceLightPosQueryResult* result) const {
        SLRAssert(light.top() == this, "Object stored in Intersection does not match this object.");
        
        m_surface->sample(smp.uPos[0], smp.uPos[1], &result->surfPt, &result->areaPDF);
        result->posType = DirectionType::LowFreq;// TODO: consider sampling delta function. and dedicated enum?
        result->surfPt.setObject(this);
        return m_material->emittance(result->surfPt, query.wls);
    }
    
    Ray SingleSurfaceObject::sampleRay(const SurfaceLight &light,
                                       const LightPosQuery &lightPosQuery, const SurfaceLightPosSample &lightPosSample,
                                       const EDFQuery &edfQuery, const EDFSample &edfSample,
                                       ArenaAllocator &mem,
                                       SurfaceLightPosQueryResult *lightPosResult, SampledSpectrum *Le0, EDF **edf,
                                       EDFQueryResult *edfResult, SampledSpectrum *Le1) const {
        SLRAssert(light.top() == this, "Object stored in Intersection does not match this object.");
        
        // sample a position with emittance on the selected light's surface.
        *Le0 = sample(light, lightPosQuery, lightPosSample, lightPosResult);
        *edf = lightPosResult->surfPt.createEDF(lightPosQuery.wls, mem);
        SLRAssert(!std::isnan(lightPosResult->areaPDF)/* && !std::isinf(lightResult)*/, "areaPDF: unexpected value detected: %f", lightPosResult->areaPDF);
        // sample a direction from EDF.
        *Le1 = (*edf)->sample(edfQuery, edfSample, edfResult);
        return Ray(lightPosResult->surfPt.getPosition(), lightPosResult->surfPt.fromLocal(edfResult->dir_sn), lightPosQuery.time, Ray::Epsilon);
    }
    
    bool SingleSurfaceObject::intersect(Ray &ray, SurfaceInteraction* si) const {
#ifdef DEBUG
        if (Accelerator::traceTraverse) {
            debugPrintf("%s%p, SingleSurfaceObject::intersect()\n",
                        Accelerator::traceTraversePrefix.c_str(), this);
            Accelerator::traceTraversePrefix += "  ";
        }
#endif
        if (!m_surface->intersect(ray, si)) {
#ifdef DEBUG
            if (Accelerator::traceTraverse) {
                debugPrintf("%snot found\n", Accelerator::traceTraversePrefix.c_str());
                size_t newLength = Accelerator::traceTraversePrefix.length() - 2;
                Accelerator::traceTraversePrefix.resize(newLength);
            }
#endif
            return false;
        }
        si->push(this);
#ifdef DEBUG
        if (Accelerator::traceTraverse) {
            debugPrintf("%sfound: %g, %g\n",
                        Accelerator::traceTraversePrefix.c_str(), ray.distMax, si->getDistance());
            size_t newLength = Accelerator::traceTraversePrefix.length() - 2;
            Accelerator::traceTraversePrefix.resize(newLength);
        }
#endif
        return true;
    }
    
    void SingleSurfaceObject::getSurfacePoint(const SurfaceInteraction &si, SurfacePoint* surfPt) const {
        m_surface->getSurfacePoint(si, surfPt);
        surfPt->setObject(this);
    }
    
    BSDF* SingleSurfaceObject::createBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const {
        return m_material->getBSDF(surfPt, wls, mem);
    }
    
    float SingleSurfaceObject::evaluateAreaPDF(const SurfacePoint& surfPt) const {
        return m_surface->evaluateAreaPDF(surfPt);
    }
    
    SampledSpectrum SingleSurfaceObject::emittance(const SurfacePoint& surfPt, const WavelengthSamples &wls) const {
        return m_material->emittance(surfPt, wls);
    }
    
    EDF* SingleSurfaceObject::createEDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const {
        return m_material->getEDF(surfPt, wls, mem);
    }
    
    
    
    void BumpSingleSurfaceObject::getSurfacePoint(const SurfaceInteraction &si, SurfacePoint *surfPt) const {
        SingleSurfaceObject::getSurfacePoint(si, surfPt);
        
        const ReferenceFrame &originalFrame = surfPt->getShadingFrame();
        
        Vector3D nLocal = m_normalMap->evaluate(*surfPt);
        Vector3D tLocal = Vector3D::Ex - dot(nLocal, Vector3D::Ex) * nLocal;
        Vector3D bLocal = Vector3D::Ey - dot(nLocal, Vector3D::Ey) * nLocal;
        Vector3D t = normalize(originalFrame.fromLocal(tLocal));
        Vector3D b = normalize(originalFrame.fromLocal(bLocal));
        Vector3D n = normalize(originalFrame.fromLocal(nLocal));
        
        ReferenceFrame bumpFrame{t, b, n};
        surfPt->setShadingFrame(bumpFrame);
    }
    
    
    
    InfiniteSphereSurfaceObject::InfiniteSphereSurfaceObject(const Scene* scene, const IBLEmission* emitter) :
    m_scene(scene) {
        m_surface = new InfiniteSphere();
        m_material = new EmitterSurfaceMaterial(nullptr, emitter);
        m_dist = emitter->createIBLImportanceMap();
    }
    
    InfiniteSphereSurfaceObject::~InfiniteSphereSurfaceObject() {
        delete m_dist;
        delete m_material;
        delete m_surface;
    }
    
    bool InfiniteSphereSurfaceObject::isEmitting() const {
        //    return m_material->isEmitting();
        return true;
    }
    
    float InfiniteSphereSurfaceObject::importance() const {
        return 1.0f;// TODO: consider a total power emitted from this object.
    }
    
    SampledSpectrum InfiniteSphereSurfaceObject::sample(const SurfaceLight &light,
                                                        const LightPosQuery &query, const SurfaceLightPosSample &smp, SurfaceLightPosQueryResult* result) const {
        SLRAssert(light.top() == this, "Object stored in Intersection does not match this object.");
        
        float uvPDF;
        float theta, phi;
        m_dist->sample(smp.uPos[0], smp.uPos[1], &phi, &theta, &uvPDF);
        phi *= 2 * M_PI;
        theta *= M_PI;
        
        Point3D pos = Point3D(-std::sin(phi) * std::sin(theta), std::cos(theta), std::cos(phi) * std::sin(theta));
        
        Vector3D texCoord0Dir = normalize(Vector3D(-std::cos(phi), 0.0f, -std::sin(phi)));
        Normal3D geometricNormal = -(Vector3D)pos;
        
        ReferenceFrame shadingFrame;
        shadingFrame.x = texCoord0Dir;
        shadingFrame.z = geometricNormal;
        shadingFrame.y = cross(shadingFrame.z, shadingFrame.x);
        SLRAssert(absDot(shadingFrame.z, shadingFrame.x) < 0.01f, "shading normal and tangent must be orthogonal.");
        
        result->surfPt = SurfacePoint(pos, // ---------------------------------------- position in world coodinate
                                      true, // --------------------------------------- atInfinity
                                      shadingFrame, // ------------------------------- shading frame
                                      geometricNormal, // ---------------------------- geometric normal in world coordinate
                                      phi, theta, // --------------------------------- surface parameter
                                      TexCoord2D(phi / (2 * M_PI), theta / M_PI), // - texture coordinate
                                      texCoord0Dir // -------------------------------- direction of texture coordinate 0
                                      );
        result->surfPt.setObject(this);
        result->posType = DirectionType::LowFreq;
        // The true value is: lim_{l to inf} uvPDF / (2 * M_PI * M_PI * std::sin(theta)) / l^2
        result->areaPDF = uvPDF / (2 * M_PI * M_PI * std::sin(theta));
        return m_material->emittance(result->surfPt, query.wls);
    }
    
    Ray InfiniteSphereSurfaceObject::sampleRay(const SurfaceLight &light,
                                               const LightPosQuery &lightPosQuery, const SurfaceLightPosSample &lightPosSample,
                                               const EDFQuery &edfQuery, const EDFSample &edfSample,
                                               ArenaAllocator &mem,
                                               SurfaceLightPosQueryResult *lightPosResult, SampledSpectrum *Le0, EDF **edf,
                                               EDFQueryResult *edfResult, SampledSpectrum *Le1) const {
        SLRAssert(light.top() == this, "Object stored in Intersection does not match this object.");
        
        // sample a position with emittance on the selected light's surface.
        *Le0 = sample(light, lightPosQuery, lightPosSample, lightPosResult);
        *edf = lightPosResult->surfPt.createEDF(lightPosQuery.wls, mem);
        SLRAssert(!std::isnan(lightPosResult->areaPDF)/* && !std::isinf(lightResult)*/, "areaPDF: unexpected value detected: %f", lightPosResult->areaPDF);
        
        // sample a direction from EDF.
        // Sampled directions for a certain sampled position (on the infinite sphere) must be parallel,
        // but be able to reach any position in the scene.
        // Therefore, it requires modification to ray's origin.
        *Le1 = (*edf)->sample(edfQuery, edfSample, edfResult);
        Vector3D vx, vy;
        Vector3D vz = lightPosResult->surfPt.fromLocal(edfResult->dir_sn);
        vz.makeCoordinateSystem(&vx, &vy);
        float dx, dy;
        concentricSampleDisk(edfSample.uDir[0], edfSample.uDir[1], &dx, &dy);
        
        float worldRadius = m_scene->getWorldRadius();
        Point3D org = m_scene->getWorldCenter() + 1.1f * worldRadius * lightPosResult->surfPt.getPosition() + worldRadius * (dx * vx + dy * vy);
        return Ray(org, vz, lightPosQuery.time, 0);
    }
    
    BSDF* InfiniteSphereSurfaceObject::createBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const {
        //SLRAssert(false, "InfiniteSphereSurfaceObject::createBSDF() should not be called.");
        //return nullptr;
        return mem.create<NullBSDF>();
    }
    
    float InfiniteSphereSurfaceObject::evaluateAreaPDF(const SurfacePoint& surfPt) const {
        float phi, theta;
        surfPt.getSurfaceParameter(&phi, &theta);
        float uvPDF = m_dist->evaluatePDF(phi / (2 * M_PI), theta / M_PI);
        return uvPDF / (2 * M_PI * M_PI * std::sin(theta));
    }
    
    
    
    BoundingBox3D TransformedSurfaceObject::bounds() const {
        return m_transform->motionBounds(m_surfObj->bounds());
    }
    
    bool TransformedSurfaceObject::isEmitting() const {
        return m_surfObj->isEmitting();
    }
    
    float TransformedSurfaceObject::importance() const {
        return m_surfObj->importance();
    }
    
    void TransformedSurfaceObject::selectLight(float u, SurfaceLight* light, float* prob) const {
        m_surfObj->selectLight(u, light, prob);
        light->push(this);
    }
    
    float TransformedSurfaceObject::evaluateProb(const SurfaceLight &light) const {
        SLRAssert(light.top() == this, "Object stored in Intersection does not match this object.");
        
        ScopedPop<const SurfaceObject*> sp = light.scopedPop();
        return m_surfObj->evaluateProb(light);
    }
    
    SampledSpectrum TransformedSurfaceObject::sample(const SurfaceLight &light,
                                                     const LightPosQuery &query, const SurfaceLightPosSample &smp, SurfaceLightPosQueryResult* result) const {
        SLRAssert(light.top() == this, "Object stored in Intersection does not match this object.");
        
        ScopedPop<const SurfaceObject*> sp = light.scopedPop();
        SampledSpectrum M = light.top()->sample(light, query, smp, result);
        StaticTransform sampledTF;
        m_transform->sample(query.time, &sampledTF);
        result->surfPt.applyTransform(sampledTF);
        return M;
    }
    
    Ray TransformedSurfaceObject::sampleRay(const SurfaceLight &light,
                                            const LightPosQuery &lightPosQuery, const SurfaceLightPosSample &lightPosSample,
                                            const EDFQuery &edfQuery, const EDFSample &edfSample,
                                            ArenaAllocator &mem,
                                            SurfaceLightPosQueryResult *lightPosResult, SampledSpectrum *Le0, EDF **edf,
                                            EDFQueryResult *edfResult, SampledSpectrum *Le1) const {
        SLRAssert(light.top() == this, "Object stored in Intersection does not match this object.");
        
        ScopedPop<const SurfaceObject*> sp = light.scopedPop();
        Ray ray = light.top()->sampleRay(light, lightPosQuery, lightPosSample, edfQuery, edfSample, mem, lightPosResult, Le0, edf, edfResult, Le1);
        StaticTransform sampledTF;
        m_transform->sample(lightPosQuery.time, &sampledTF);
        lightPosResult->surfPt.applyTransform(sampledTF);
        return sampledTF * ray;
    }
    
    bool TransformedSurfaceObject::intersect(Ray &ray, SurfaceInteraction* si) const {
#ifdef DEBUG
        if (Accelerator::traceTraverse) {
            debugPrintf("%s%p, TransformedSurfaceObject::intersect()\n",
                        Accelerator::traceTraversePrefix.c_str(), this);
            Accelerator::traceTraversePrefix += "  ";
        }
#endif
        Ray localRay;
        StaticTransform sampledTF;
        m_transform->sample(ray.time, &sampledTF);
        localRay = invert(sampledTF) * ray;
        if (!m_surfObj->intersect(localRay, si)) {
#ifdef DEBUG
            if (Accelerator::traceTraverse) {
                debugPrintf("%snot found\n", Accelerator::traceTraversePrefix.c_str());
                size_t newLength = Accelerator::traceTraversePrefix.length() - 2;
                Accelerator::traceTraversePrefix.resize(newLength);
            }
#endif
            return false;
        }
        ray.distMax = localRay.distMax;
        si->push(this);
#ifdef DEBUG
        if (Accelerator::traceTraverse) {
            debugPrintf("%sfound: %g\n",
                        Accelerator::traceTraversePrefix.c_str(), ray.distMax);
            size_t newLength = Accelerator::traceTraversePrefix.length() - 2;
            Accelerator::traceTraversePrefix.resize(newLength);
        }
#endif
        return true;
    }
    
    void TransformedSurfaceObject::getSurfacePoint(const SurfaceInteraction &si, SurfacePoint* surfPt) const {
        SLRAssert(si.top() == this, "Object stored in Intersection does not match this object.");
        
        ScopedPop<const SurfaceObject*> sp = si.scopedPop();
        si.getSurfacePoint(surfPt);
        StaticTransform sampledTF;
        m_transform->sample(si.getTime(), &sampledTF);
        surfPt->applyTransform(sampledTF);
    }
    
    
    
    SurfaceObjectAggregate::SurfaceObjectAggregate(std::vector<SurfaceObject*> &objs) {
//        SBVH sbvh(objs);
//        m_accelerator = new QBVH(sbvh);
        m_accelerator = new SBVH(objs);
//        m_accelerator = new StandardBVH(objs, StandardBVH::Partitioning::BinnedSAH);
        
        std::vector<const SurfaceObject*> lights;
        std::vector<float> lightImportances;
        for (int i = 0; i < objs.size(); ++i) {
            const SurfaceObject* obj = objs[i];
            if (obj->isEmitting()) {
                lights.push_back(obj);
                lightImportances.push_back(obj->importance());
            }
        }
        
        m_lightList = new const SurfaceObject*[lightImportances.size()];
        m_lightDist1D = new RegularConstantDiscrete1D(lightImportances);
        
        for (int i = 0; i < lights.size(); ++i) {
            const SurfaceObject* light = lights[i];
            m_lightList[i] = light;
            m_revMap[light] = i;
        }
    }
    
    SurfaceObjectAggregate::~SurfaceObjectAggregate() {
        delete m_accelerator;
        
        delete m_lightList;
        delete m_lightDist1D;
    };
    
    BoundingBox3D SurfaceObjectAggregate::bounds() const {
        return m_accelerator->bounds();
    }
    
    bool SurfaceObjectAggregate::isEmitting() const {
        return m_revMap.size() > 0;
    }
    
    float SurfaceObjectAggregate::importance() const {
        return m_lightDist1D->integral();
    }
    
    void SurfaceObjectAggregate::selectLight(float u, SurfaceLight* light, float* prob) const {
        uint32_t lIdx = m_lightDist1D->sample(u, prob, &u);
        const SurfaceObject* obj = m_lightList[lIdx];
        float cProb;
        obj->selectLight(u, light, &cProb);
        *prob *= cProb;
        light->push(this);
    }
    
    float SurfaceObjectAggregate::evaluateProb(const SurfaceLight &light) const {
        SLRAssert(light.top() == this, "Object stored in Intersection does not match this object.");
        
        ScopedPop<const SurfaceObject*> sp = light.scopedPop();
        SLRAssert(m_revMap.count(light.top()) > 0, "Specified light is not found.");
        float prob = m_lightDist1D->evaluatePMF(m_revMap.at(light.top()));
        return prob * light.top()->evaluateProb(light);
    }
    
    float SurfaceObjectAggregate::costForIntersect() const {
        return m_accelerator->costForIntersect();
    }
    
    bool SurfaceObjectAggregate::intersect(Ray &ray, SurfaceInteraction* si) const {
#ifdef DEBUG
        if (Accelerator::traceTraverse) {
            debugPrintf("%s%p, SurfaceObjectAggregate::intersect()\n",
                        Accelerator::traceTraversePrefix.c_str(), this);
            Accelerator::traceTraversePrefix += "  ";
        }
#endif
        if (!m_accelerator->intersect(ray, si)) {
#ifdef DEBUG
            if (Accelerator::traceTraverse) {
                debugPrintf("%snot found\n", Accelerator::traceTraversePrefix.c_str());
                size_t newLength = Accelerator::traceTraversePrefix.length() - 2;
                Accelerator::traceTraversePrefix.resize(newLength);
            }
#endif
            return false;
        }
        si->push(this);
#ifdef DEBUG
        if (Accelerator::traceTraverse) {
            debugPrintf("%sfound: %g\n",
                        Accelerator::traceTraversePrefix.c_str(), ray.distMax);
            size_t newLength = Accelerator::traceTraversePrefix.length() - 2;
            Accelerator::traceTraversePrefix.resize(newLength);
        }
#endif
        return true;
    }
}

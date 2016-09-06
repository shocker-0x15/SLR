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
#include "../Accelerator/StandardBVH.h"
#include "../Accelerator/SBVH.h"
#include "../Accelerator/QBVH.h"
#include "textures.h"
#include "../Surface/InfiniteSphere.h"
#include "../SurfaceMaterials/IBLEmission.h"
#include "../Memory/ArenaAllocator.h"
#include "../BSDFs/basic_BSDFs.h"

namespace SLR {
    SampledSpectrum Light::sample(const LightPosQuery &query, const LightPosSample &smp, LightPosQueryResult* result) const {
        return m_hierarchy.top()->sample(*this, query, smp, result);
    }
    
    Ray Light::sampleRay(const LightPosQuery &lightPosQuery, const LightPosSample &lightPosSample, LightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                         const EDFQuery &edfQuery, const EDFSample &edfSample, EDFQueryResult* edfResult, SampledSpectrum* Le1,
                         ArenaAllocator &mem) const {
        return m_hierarchy.top()->sampleRay(*this, lightPosQuery, lightPosSample, lightPosResult, Le0, edf, edfQuery, edfSample, edfResult, Le1, mem);
    }
    
    
    
    bool SurfaceObject::intersect(Ray &ray, SurfacePoint *surfPt) const {
        Intersection isect;
        if (!intersect(ray, &isect))
            return false;
        getSurfacePoint(isect, surfPt);
        return true;
    }
    
    bool SurfaceObject::testVisibility(const SurfacePoint &shdP, const SurfacePoint &lightP, float time) const {
        SLRAssert(shdP.atInfinity == false && lightP.atInfinity == false, "Points must be in finite region.");
        float dist = distance(lightP.p, shdP.p);
        Ray ray(shdP.p, (lightP.p - shdP.p) / dist, time, Ray::Epsilon, dist * (1 - Ray::Epsilon));
        Intersection isect;
        return !intersect(ray, &isect);
    }
    
    
    
    bool SingleSurfaceObject::intersect(Ray &ray, Intersection* isect) const {
        if (!m_surface->intersect(ray, isect))
            return false;
        isect->time = ray.time;
        isect->obj.push(this);
        return true;
    }
    
    void SingleSurfaceObject::getSurfacePoint(const Intersection &isect, SurfacePoint* surfPt) const {
        m_surface->getSurfacePoint(isect, surfPt);
        surfPt->obj = this;
    }
    
    bool SingleSurfaceObject::isEmitting() const {
        return m_material->isEmitting();
    }
    
    float SingleSurfaceObject::importance() const {
        return 1.0f;// TODO: consider a total power emitted from this object.
    }
    
    void SingleSurfaceObject::selectLight(float u, Light* light, float* prob) const {
        light->push(this);
        *prob = 1.0f;
    }
    
    float SingleSurfaceObject::evaluateProb(const Light &light) const {
        return light.top() == this ? 1.0f : 0.0f;
    }
    
    SampledSpectrum SingleSurfaceObject::sample(const Light &light, const LightPosQuery &query, const LightPosSample &smp, LightPosQueryResult* result) const {
        if (light.top()!= this) {
            result->areaPDF = 0.0f;
            return SampledSpectrum::Zero;
        }
        m_surface->sample(smp.uPos[0], smp.uPos[1], &result->surfPt, &result->areaPDF);
        result->posType = DirectionType::LowFreq;// TODO: consider sampling delta function. and dedicated enum?
        result->surfPt.obj = this;
        return m_material->emittance(result->surfPt, query.wls);
    }
    
    Ray SingleSurfaceObject::sampleRay(const Light &light,
                                       const LightPosQuery &lightPosQuery, const LightPosSample &lightPosSample, LightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                                       const EDFQuery &edfQuery, const EDFSample &edfSample, EDFQueryResult* edfResult, SampledSpectrum* Le1,
                                       ArenaAllocator &mem) const {
        // sample a position with emittance on the selected light's surface.
        *Le0 = sample(light, lightPosQuery, lightPosSample, lightPosResult);
        *edf = lightPosResult->surfPt.createEDF(lightPosQuery.wls, mem);
        SLRAssert(!std::isnan(lightPosResult->areaPDF)/* && !std::isinf(lightResult)*/, "areaPDF: unexpected value detected: %f", lightPosResult->areaPDF);
        // sample a direction from EDF.
        *Le1 = (*edf)->sample(edfQuery, edfSample, edfResult);
        return Ray(lightPosResult->surfPt.p, lightPosResult->surfPt.shadingFrame.fromLocal(edfResult->dir_sn), lightPosQuery.time, Ray::Epsilon);
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
    
    
    void BumpSingleSurfaceObject::getSurfacePoint(const Intersection &isect, SurfacePoint *surfPt) const {
        SingleSurfaceObject::getSurfacePoint(isect, surfPt);
        Vector3D nLocal = m_normalMap->evaluate(*surfPt);
        Vector3D tLocal = Vector3D::Ex - dot(nLocal, Vector3D::Ex) * nLocal;
        Vector3D bLocal = Vector3D::Ey - dot(nLocal, Vector3D::Ey) * nLocal;
        Vector3D t = normalize(surfPt->shadingFrame.fromLocal(tLocal));
        Vector3D b = normalize(surfPt->shadingFrame.fromLocal(bLocal));
        Vector3D n = normalize(surfPt->shadingFrame.fromLocal(nLocal));
        surfPt->shadingFrame.x = t;
        surfPt->shadingFrame.y = b;
        surfPt->shadingFrame.z = n;
    }
    
    
    InfiniteSphereSurfaceObject::InfiniteSphereSurfaceObject(const Scene* scene, const IBLEmission* emitter) : m_scene(scene) {
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
    
    SampledSpectrum InfiniteSphereSurfaceObject::sample(const Light &light, const LightPosQuery &query, const LightPosSample &smp, LightPosQueryResult* result) const {
        if (light.top()!= this) {
            result->areaPDF = 0.0f;
            return SampledSpectrum::Zero;
        }
        float uvPDF;
        float theta, phi;
        m_dist->sample(smp.uPos[0], smp.uPos[1], &phi, &theta, &uvPDF);
        phi *= 2 * M_PI;
        theta *= M_PI;
        SurfacePoint &surfPt = result->surfPt;
        surfPt.p = Point3D(-std::sin(phi) * std::sin(theta), std::cos(theta), std::cos(phi) * std::sin(theta));
        surfPt.atInfinity = true;
        surfPt.gNormal = -(Vector3D)surfPt.p;
        surfPt.u = phi;
        surfPt.v = theta;
        surfPt.texCoord = TexCoord2D(phi / (2 * M_PI), theta / M_PI);
        surfPt.texCoord0Dir = normalize(Vector3D(-std::cos(phi), 0.0f, -std::sin(phi)));
        surfPt.shadingFrame.x = surfPt.texCoord0Dir;
        surfPt.shadingFrame.z = surfPt.gNormal;
        surfPt.shadingFrame.y = cross(surfPt.shadingFrame.z, surfPt.shadingFrame.x);
        SLRAssert(absDot(surfPt.shadingFrame.z, surfPt.shadingFrame.x) < 0.01f, "shading normal and tangent must be orthogonal.");
        surfPt.obj = this;
        result->posType = DirectionType::LowFreq;
        // The true value is: lim_{l to inf} uvPDF / (2 * M_PI * M_PI * std::sin(theta)) / l^2
        result->areaPDF = uvPDF / (2 * M_PI * M_PI * std::sin(theta));
        return m_material->emittance(result->surfPt, query.wls);
    }
    
    Ray InfiniteSphereSurfaceObject::sampleRay(const Light &light,
                                               const LightPosQuery &lightPosQuery, const LightPosSample &lightPosSample, LightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                                               const EDFQuery &edfQuery, const EDFSample &edfSample, EDFQueryResult* edfResult, SampledSpectrum* Le1,
                                               ArenaAllocator &mem) const {
        // sample a position with emittance on the selected light's surface.
        *Le0 = light.sample(lightPosQuery, lightPosSample, lightPosResult);
        *edf = lightPosResult->surfPt.createEDF(lightPosQuery.wls, mem);
        SLRAssert(!std::isnan(lightPosResult->areaPDF)/* && !std::isinf(lightResult)*/, "areaPDF: unexpected value detected: %f", lightPosResult->areaPDF);
        
        // sample a direction from EDF.
        // Sampled directions for a certain sampled position (on the infinite sphere) must be parallel,
        // but be able to reach any position in the scene.
        // Therefore, it requires modification to ray's origin.
        *Le1 = (*edf)->sample(edfQuery, edfSample, edfResult);
        Vector3D vx, vy;
        Vector3D vz = lightPosResult->surfPt.shadingFrame.fromLocal(edfResult->dir_sn);
        vz.makeCoordinateSystem(&vx, &vy);
        float dx, dy;
        concentricSampleDisk(edfSample.uDir[0], edfSample.uDir[1], &dx, &dy);
        
        float worldRadius = m_scene->getWorldRadius();
        return Ray(m_scene->getWorldCenter() + 1.1f * worldRadius * lightPosResult->surfPt.p + worldRadius * (dx * vx + dy * vy), vz, lightPosQuery.time, 0);
    }
    
    BSDF* InfiniteSphereSurfaceObject::createBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const {
        //SLRAssert(false, "InfiniteSphereSurfaceObject::createBSDF() should not be called.");
        //return nullptr;
        return mem.create<NullBSDF>();
    }
    
    float InfiniteSphereSurfaceObject::evaluateAreaPDF(const SurfacePoint& surfPt) const {
        float phi = surfPt.u;
        float theta = surfPt.v;
        float uvPDF = m_dist->evaluatePDF(phi / (2 * M_PI), theta / M_PI);
        return uvPDF / (2 * M_PI * M_PI * std::sin(theta));
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
    
    float SurfaceObjectAggregate::costForIntersect() const {
        return m_accelerator->costForIntersect();
    }
    
    BoundingBox3D SurfaceObjectAggregate::bounds() const {
        return m_accelerator->bounds();
    }
    
    bool SurfaceObjectAggregate::intersect(Ray &ray, Intersection* isect) const {
        return m_accelerator->intersect(ray, isect);
    }
    
    bool SurfaceObjectAggregate::isEmitting() const {
        return m_revMap.size() > 0;
    }
    
    float SurfaceObjectAggregate::importance() const {
        return m_lightDist1D->integral();
    }
    
    void SurfaceObjectAggregate::selectLight(float u, Light* light, float* prob) const {
        uint32_t lIdx = m_lightDist1D->sample(u, prob, &u);
        const SurfaceObject* obj = m_lightList[lIdx];
        float cProb;
        obj->selectLight(u, light, &cProb);
        *prob *= cProb;
        //    light->push(this);
    }
    
    float SurfaceObjectAggregate::evaluateProb(const Light &light) const {
        // if (light.top() != this)
        //     return 0.0f;
        // light.pop();
        // float prob = calcProb(light.top()) * light.top()->evaluateProb(light);
        // light.push(this);
        // return prob;
        if (m_revMap.count(light.top()) == 0)
            return 0.0f;
        float prob = m_lightDist1D->evaluatePMF(m_revMap.at(light.top()));
        return prob * light.top()->evaluateProb(light);
    }
    
    
    
    BoundingBox3D TransformedSurfaceObject::bounds() const {
        return m_transform->motionBounds(m_surfObj->bounds());
    }
    
    bool TransformedSurfaceObject::intersect(Ray &ray, Intersection *isect) const {
        Ray localRay;
        StaticTransform sampledTF;
        m_transform->sample(ray.time, &sampledTF);
        localRay = invert(sampledTF) * ray;
        if (!m_surfObj->intersect(localRay, isect))
            return false;
        ray.distMax = localRay.distMax;
        isect->obj.push(this);
        return true;
    }
    
    Point3D TransformedSurfaceObject::getIntersectionPoint(const Intersection &isect) const {
        isect.obj.pop();
        Point3D ret = isect.obj.top()->getIntersectionPoint(isect);
        StaticTransform sampledTF;
        m_transform->sample(isect.time, &sampledTF);
        ret = sampledTF * ret;
        isect.obj.push(this);
        return ret;
    }
    
    void TransformedSurfaceObject::getSurfacePoint(const Intersection &isect, SurfacePoint *surfPt) const {
        isect.obj.pop();
        isect.obj.top()->getSurfacePoint(isect, surfPt);
        StaticTransform sampledTF;
        m_transform->sample(isect.time, &sampledTF);
        *surfPt = sampledTF * *surfPt;
        isect.obj.push(this);
    }
    
    bool TransformedSurfaceObject::isEmitting() const {
        return m_surfObj->isEmitting();
    }
    
    float TransformedSurfaceObject::importance() const {
        return m_surfObj->importance();
    }
    
    void TransformedSurfaceObject::selectLight(float u, Light* light, float* prob) const {
        m_surfObj->selectLight(u, light, prob);
        light->push(this);
    }
    
    float TransformedSurfaceObject::evaluateProb(const Light &light) const {
        if (light.top() != this)
            return 0.0f;
        light.pop();
        float prob = m_surfObj->evaluateProb(light);
        light.push(this);
        return prob;
    }
    
    SampledSpectrum TransformedSurfaceObject::sample(const Light &light, const LightPosQuery &query, const LightPosSample &smp, LightPosQueryResult* result) const {
        if (light.top() != this) {
            result->areaPDF = 0.0f;
            return SampledSpectrum::Zero;
        }
        light.pop();
        SampledSpectrum M = light.top()->sample(light, query, smp, result);
        StaticTransform sampledTF;
        m_transform->sample(query.time, &sampledTF);
        result->surfPt = sampledTF * result->surfPt;
        light.push(this);
        return M;
    }
    
    Ray TransformedSurfaceObject::sampleRay(const Light &light,
                                            const LightPosQuery &lightPosQuery, const LightPosSample &lightPosSample, LightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                                            const EDFQuery &edfQuery, const EDFSample &edfSample, EDFQueryResult* edfResult, SampledSpectrum* Le1,
                                            ArenaAllocator &mem) const {
        if (light.top() != this) {
            lightPosResult->areaPDF = 0.0f;
            *Le0 = SampledSpectrum::Zero;
            edfResult->dirPDF = 0.0f;
            *Le1 = SampledSpectrum::Zero;
            return Ray();
        }
        light.pop();
        Ray ray = light.top()->sampleRay(light, lightPosQuery, lightPosSample, lightPosResult, Le0, edf, edfQuery, edfSample, edfResult, Le1, mem);
        StaticTransform sampledTF;
        m_transform->sample(lightPosQuery.time, &sampledTF);
        lightPosResult->surfPt = sampledTF * lightPosResult->surfPt;
        light.push(this);
        return sampledTF * ray;
    }
    
    
    
    void Scene::build(const SurfaceObjectAggregate* aggregate, const InfiniteSphereSurfaceObject* envSphere, const Camera* camera) {
        m_aggregate = aggregate;
        m_envSphere = envSphere;
        m_camera = camera;
        
        // calculate world bounding sphere and store its radius and disc area.
        BoundingBox3D worldBounds = m_aggregate->bounds();
        m_worldCenter = worldBounds.centroid();
        m_worldRadius = (worldBounds.maxP - m_worldCenter).length();
        m_worldDiscArea = M_PI * m_worldRadius * m_worldRadius;
    };
    
    bool Scene::intersect(Ray &ray, Intersection *isect) const {
        if (m_aggregate->intersect(ray, isect))
            return true;
        if (m_envSphere) {
            if (m_envSphere->intersect(ray, isect))
                return true;
        }
        return false;
    }
    
    bool Scene::testVisibility(const SurfacePoint &shdP, const SurfacePoint &lightP, float time) const {
        SLRAssert(shdP.atInfinity == false, "Shading point must be in finite region.");
        Ray ray;
        if (lightP.atInfinity) {
            ray = Ray(shdP.p, normalize(lightP.p - Point3D::Zero), time, Ray::Epsilon, FLT_MAX);
        }
        else {
            float dist = distance(lightP.p, shdP.p);
            ray = Ray(shdP.p, (lightP.p - shdP.p) / dist, time, Ray::Epsilon, dist * (1 - Ray::Epsilon));
        }
        Intersection isect;
        return !intersect(ray, &isect);
    }
    
    void Scene::selectLight(float u, Light* light, float* prob) const {
        if (m_envSphere) {
            float sumImps = m_aggregate->importance() + m_envSphere->importance();
            float su = sumImps * u;
            if (su < m_aggregate->importance()) {
                u = u / (m_aggregate->importance() / sumImps);
                m_aggregate->selectLight(u, light, prob);
                *prob *= m_aggregate->importance() / sumImps;
            }
            else {
                u = (u - m_aggregate->importance()) / (m_envSphere->importance() / sumImps);
                m_envSphere->selectLight(u, light, prob);
                *prob *= m_envSphere->importance() / sumImps;
            }
        }
        else {
            m_aggregate->selectLight(u, light, prob);
        }
    }
    
    float Scene::evaluateProb(const Light &light) const {
        float ret = 0.0f;
        if (m_envSphere) {
            float sumImps = m_aggregate->importance() + m_envSphere->importance();
            if (light.top() == m_envSphere)
                ret = m_envSphere->importance() / sumImps;
            else
                ret = m_aggregate->importance()  / sumImps * m_aggregate->evaluateProb(light);
        }
        else {
            ret = m_aggregate->evaluateProb(light);
        }
        SLRAssert(!std::isnan(ret) && !std::isinf(ret), "%g", ret);
        return ret;
    }
}

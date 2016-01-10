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
#include "../Accelerator/BBVH.h"
#include "textures.h"
#include "../SurfaceMaterials/IBLEmission.h"

namespace SLR {
    SampledSpectrum Light::sample(const LightPosQuery &query, const LightPosSample &smp, LightPosQueryResult* result) const {
        return m_hierarchy.top()->sample(*this, query, smp, result);
    }
    
    
    bool SurfaceObject::intersect(Ray &ray, SurfacePoint *surfPt) const {
        Intersection isect;
        if (!intersect(ray, &isect))
            return false;
        getSurfacePoint(isect, surfPt);
        return true;
    }
    
    bool SurfaceObject::testVisiblility(const SurfacePoint &shdP, const SurfacePoint &lightP, float time) const {
        SLRAssert(shdP.atInfinity == false && lightP.atInfinity == false, "Points must be in finite region.");
        float dist = distance(lightP.p, shdP.p);
        Ray ray(shdP.p, (lightP.p - shdP.p) / dist, time, Ray::Epsilon, dist * (1 - Ray::Epsilon));
        Intersection isect;
        return !intersect(ray, &isect);
    }
    
    
    BoundingBox3D SingleSurfaceObject::bounds() const {
        return m_surface->bounds();
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
        result->isDeltaPos = false;// TODO: consider sampling delta function.
        result->surfPt.obj = this;
        return m_material->emittance(result->surfPt, query.wls);
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
        Vector3D nLocal = m_normalMap->evaluate(surfPt->texCoord);
        const Vector3D rawSN = surfPt->shadingFrame.z;
        const Vector3D tc0Dir = normalize(surfPt->texCoord0Dir - dot(rawSN, surfPt->texCoord0Dir) * rawSN);
        const Vector3D tc1Dir = cross(rawSN, tc0Dir);
        Vector3D n = Vector3D(dot(Vector3D(tc0Dir.x, tc1Dir.x, rawSN.x), nLocal),
                              dot(Vector3D(tc0Dir.y, tc1Dir.y, rawSN.y), nLocal),
                              dot(Vector3D(tc0Dir.z, tc1Dir.z, rawSN.z), nLocal));
        Vector3D t = normalize(surfPt->shadingFrame.x - dot(n, surfPt->shadingFrame.x) * n);
        Vector3D b = cross(n, t);
        surfPt->shadingFrame.x = t;
        surfPt->shadingFrame.y = b;
        surfPt->shadingFrame.z = n;
    }
    
    
    InfiniteSphereSurfaceObject::InfiniteSphereSurfaceObject(const Surface* surf, const IBLEmission* emitter) {
        m_surface = surf;
        m_material = new EmitterSurfaceMaterial(nullptr, emitter);
        m_dist = emitter->createIBLImportanceMap();
    }
    
    InfiniteSphereSurfaceObject::~InfiniteSphereSurfaceObject() {
        delete m_dist;
        delete m_material;
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
        surfPt.texCoord0Dir = normalize(Vector3D(-std::cos(phi), 0.0f, -std::sin(theta)));
        surfPt.shadingFrame.x = surfPt.texCoord0Dir;
        surfPt.shadingFrame.z = surfPt.gNormal;
        surfPt.shadingFrame.y = cross(surfPt.shadingFrame.z, surfPt.shadingFrame.x);
        surfPt.obj = this;
        result->isDeltaPos = false;
        result->areaPDF = uvPDF / (2 * M_PI * M_PI * std::sin(theta));
        return m_material->emittance(result->surfPt, query.wls);
    }
    
    BSDF* InfiniteSphereSurfaceObject::createBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const {
        SLRAssert(false, "InfiniteSphereSurfaceObject::createBSDF() should not be called.");
        return nullptr;
    }
    
    float InfiniteSphereSurfaceObject::evaluateAreaPDF(const SurfacePoint& surfPt) const {
        float phi = surfPt.u;
        float theta = surfPt.v;
        float uvPDF = m_dist->evaluatePDF(phi / (2 * M_PI), theta / M_PI);
        return uvPDF / (2 * M_PI * M_PI * std::sin(theta));
    }
    
    
    SurfaceObjectAggregate::SurfaceObjectAggregate(std::vector<SurfaceObject*> &objs) {
        m_accelerator = new BBVH(objs, BBVH::Partitioning::BinnedSAH);
        
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
    
    bool SurfaceObjectAggregate::intersect(Ray &ray, Intersection* isect) const {
        return m_accelerator->intersect(ray, isect);
    }
    
    void SurfaceObjectAggregate::getSurfacePoint(const Intersection &isect, SurfacePoint* surfPt) const {
        isect.obj.top()->getSurfacePoint(isect, surfPt);
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
        // 	return 0.0f;
        // light.pop();
        // float prob = calcProb(light.top()) * light.top()->evaluateProb(light);
        // light.push(this);
        // return prob;
        if (m_revMap.count(light.top()) == 0)
            return 0.0f;
        float prob = m_lightDist1D->evaluatePMF(m_revMap.at(light.top()));
        return prob * light.top()->evaluateProb(light);
    }
    
    SampledSpectrum SurfaceObjectAggregate::sample(const Light &light, const LightPosQuery &query, const LightPosSample &smp, LightPosQueryResult* result) const {
        //    if (light.top() != this) {
        //        result->areaPDF = 0.0f;
        //        return SampledSpectrum::Zero;
        //    }
        //    light.pop();
        //    SampledSpectrum M = light.top()->sample(light, query, result);
        //    light.push(this);
        //    return M;
        return light.top()->sample(light, query, smp, result);
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
    
    bool Scene::testVisiblility(const SurfacePoint &shdP, const SurfacePoint &lightP, float time) const {
        SLRAssert(shdP.atInfinity == false, "Shading point must be in finite region.");
        Ray ray;
        if (lightP.atInfinity) {
            ray = Ray(shdP.p, normalize(lightP.p - Point3D::Zero), time, Ray::Epsilon, m_worldRadius);
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
        if (m_envSphere) {
            float sumImps = m_aggregate->importance() + m_envSphere->importance();
            if (light.top() == m_envSphere) {
                return m_envSphere->importance() / sumImps;
            }
            else {
                return m_aggregate->importance()  / sumImps * m_aggregate->evaluateProb(light);
            }
        }
        else {
            return m_aggregate->evaluateProb(light);
        }
    }
}

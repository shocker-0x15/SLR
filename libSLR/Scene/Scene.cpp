//
//  Scene.cpp
//
//  Created by 渡部 心 on 2016/09/19.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "Scene.h"
#include "../Memory/ArenaAllocator.h"
#include "../Core/distributions.h"
#include "../Core/Transform.h"
#include "../Core/cameras.h"
#include "../Core/light_path_samplers.h"

namespace SLR {
    void Scene::build(Allocator* sceneMem) {
        m_sceneMem = sceneMem;
        
        RenderingData renderingData(this);
        m_rootNode->createRenderingData(sceneMem, nullptr, &renderingData);
        if (m_envNode)
            m_envNode->createRenderingData(sceneMem, nullptr, &renderingData);
        
        m_surfaceAggregate = sceneMem->create<SurfaceObjectAggregate>(renderingData.surfObjs);
        m_mediumAggregate = sceneMem->create<MediumObjectAggregate>(renderingData.medObjs);
        m_envSphere = m_envNode ? renderingData.envObj : nullptr;
        
        m_camera = renderingData.camera;
        if (m_camera)
            m_camera->setTransform(renderingData.camTransform);
        
        // calculate world bounding sphere and store its radius and disc area.
        BoundingBox3D worldBounds = m_surfaceAggregate->bounds();
        m_worldCenter = worldBounds.centroid();
        m_worldRadius = (worldBounds.maxP - m_worldCenter).length();
        m_worldDiscArea = M_PI * m_worldRadius * m_worldRadius;
    };
    
    void Scene::destory() {
        m_sceneMem->destroy(m_mediumAggregate);
        m_sceneMem->destroy(m_surfaceAggregate);
        if (m_envNode)
            m_envNode->destroyRenderingData(m_sceneMem);
        m_rootNode->destroyRenderingData(m_sceneMem);
    }
    
    bool Scene::intersect(Ray &ray, SurfaceInteraction *si) const {
        if (m_surfaceAggregate->intersect(ray, si))
            return true;
        if (m_envSphere) {
            if (m_envSphere->intersect(ray, si))
                return true;
        }
        return false;
    }
    
    bool Scene::interact(Ray &ray, const WavelengthSamples &wls, LightPathSampler &pathSampler, ArenaAllocator &mem,
                         Interaction** interact, SampledSpectrum* medThroughput, bool* singleWavelength) const {
        SurfaceInteraction si;
        bool hitSurface = m_surfaceAggregate->intersect(ray, &si);
        MediumInteraction mi;
        bool hitMedium = m_mediumAggregate->interact(ray, wls, pathSampler, &mi, medThroughput, singleWavelength);
        if (hitMedium) {
            if (interact)
                *interact = mem.create<MediumInteraction>(mi);
            return true;
        }
        if (hitSurface) {
            if (interact)
                *interact = mem.create<SurfaceInteraction>(si);
            return true;
        }
        if (m_envSphere) {
            if (m_envSphere->intersect(ray, &si)) {
                if (interact)
                    *interact = mem.create<SurfaceInteraction>(si);
                return true;
            }
        }
        return false;
    }
    
    bool Scene::testVisibility(const SurfacePoint &shdP, const SurfacePoint &lightP, float time) const {
        SLRAssert(shdP.atInfinity() == false, "Shading point must be in finite region.");
        Ray ray;
        if (lightP.atInfinity()) {
            ray = Ray(shdP.getPosition(), normalize(lightP.getPosition() - Point3D::Zero), time, Ray::Epsilon, FLT_MAX);
        }
        else {
            float dist = distance(lightP.getPosition(), shdP.getPosition());
            ray = Ray(shdP.getPosition(), (lightP.getPosition() - shdP.getPosition()) / dist, time, Ray::Epsilon, dist * (1 - Ray::Epsilon));
        }
        SurfaceInteraction si;
        return !intersect(ray, &si);
    }
    
    bool Scene::testVisibility(const InteractionPoint* shdP, const InteractionPoint* lightP, float time,
                               const WavelengthSamples &wls, LightPathSampler &pathSampler, SampledSpectrum* fractionalVisibility, bool* singleWavelength) const {
        SLRAssert(shdP->atInfinity() == false, "Shading point must be in finite region.");
        Ray ray;
        if (lightP->atInfinity()) {
            ray = Ray(shdP->getPosition(), normalize(lightP->getPosition() - Point3D::Zero), time, Ray::Epsilon, FLT_MAX);
        }
        else {
            float dist = distance(lightP->getPosition(), shdP->getPosition());
            ray = Ray(shdP->getPosition(), (lightP->getPosition() - shdP->getPosition()) / dist, time, Ray::Epsilon, dist * (1 - Ray::Epsilon));
        }
        *fractionalVisibility = SampledSpectrum::Zero;
        SurfaceInteraction si;
        if (m_surfaceAggregate->intersect(ray, &si))
            return false;
        *fractionalVisibility = m_mediumAggregate->evaluateTransmittance(ray, wls, pathSampler, singleWavelength);
        return true;
    }
    
    void Scene::selectSurfaceLight(float u, SurfaceLight* light, float* prob) const {
        if (m_envSphere) {
            float sumImps = m_surfaceAggregate->importance() + m_envSphere->importance();
            float su = sumImps * u;
            if (su < m_surfaceAggregate->importance()) {
                u = u / (m_surfaceAggregate->importance() / sumImps);
                m_surfaceAggregate->selectLight(u, light, prob);
                *prob *= m_surfaceAggregate->importance() / sumImps;
            }
            else {
                u = (u - m_surfaceAggregate->importance()) / (m_envSphere->importance() / sumImps);
                m_envSphere->selectLight(u, light, prob);
                *prob *= m_envSphere->importance() / sumImps;
            }
        }
        else {
            m_surfaceAggregate->selectLight(u, light, prob);
        }
    }
    
    float Scene::evaluateSurfaceLightProb(const SurfaceLight &light) const {
        float ret = 0.0f;
        if (m_envSphere) {
            float sumImps = m_surfaceAggregate->importance() + m_envSphere->importance();
            if (light.top() == m_envSphere)
                ret = m_envSphere->importance() / sumImps;
            else
                ret = m_surfaceAggregate->importance()  / sumImps * m_surfaceAggregate->evaluateProb(light);
        }
        else {
            ret = m_surfaceAggregate->evaluateProb(light);
        }
        SLRAssert(std::isfinite(ret), "%g", ret);
        return ret;
    }
    
    void Scene::selectLight(float u, ArenaAllocator &mem, Light** light, float *prob) const {
        float importances[3] = {m_surfaceAggregate->importance(), m_mediumAggregate->importance(), 0.0f};
        if (m_envSphere)
            importances[2] = m_envSphere->importance();
        float sumImportances, base;
        uint32_t idx = sampleDiscrete(importances, &sumImportances, &base, 3, u);
        
        float prob1st = importances[idx] / sumImportances;
        u = (u * sumImportances - base) / importances[idx];
        switch (idx) {
            case 0: {
                SurfaceLight* surfLight = mem.create<SurfaceLight>();
                m_surfaceAggregate->selectLight(u, surfLight, prob);
                *light = surfLight;
                break;
            }
            case 1: {
                VolumetricLight* volLight = mem.create<VolumetricLight>();
                m_mediumAggregate->selectLight(u, volLight, prob);
                *light = volLight;
                break;
            }
            case 2: {
                SurfaceLight* surfLight = mem.create<SurfaceLight>();
                m_envSphere->selectLight(u, surfLight, prob);
                *light = surfLight;
                break;
            }
            default:
                break;
        }
        *prob *= prob1st;
    }
    
    float Scene::evaluateLightProb(const Light* light) const {
        float importances[3] = {m_surfaceAggregate->importance(), m_mediumAggregate->importance(), 0.0f};
        if (m_envSphere)
            importances[2] = m_envSphere->importance();
        if (m_surfaceAggregate->contains(light))
            return m_surfaceAggregate->evaluateProbability(light) * evaluateProbability(importances, 3, 0);
        else if (m_mediumAggregate->contains(light))
            return m_mediumAggregate->evaluateProbability(light) * evaluateProbability(importances, 3, 1);
        else
            return /*m_envSphere->evaluateProb(light)*/1.0f * evaluateProbability(importances, 3, 2);
    }
}

//
//  Scene.cpp
//
//  Created by 渡部 心 on 2016/09/19.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "Scene.h"

#include "../MemoryAllocators/ArenaAllocator.h"
#include "../Core/distributions.h"
#include "../Core/transform.h"
#include "../Core/camera.h"
#include "../Core/light_path_sampler.h"

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
        BoundingBox3D worldBounds = calcUnion(m_surfaceAggregate->bounds(), m_mediumAggregate->bounds());
        m_worldCenter = worldBounds.isValid() ? worldBounds.centroid() : Point3D::Zero;
        m_worldRadius = worldBounds.isValid() ? (worldBounds.maxP - m_worldCenter).length() : 0.0f;
        m_worldDiscArea = M_PI * m_worldRadius * m_worldRadius;
    }
    
    void Scene::destory() {
        m_sceneMem->destroy(m_mediumAggregate);
        m_sceneMem->destroy(m_surfaceAggregate);
        if (m_envNode)
            m_envNode->destroyRenderingData(m_sceneMem);
        m_rootNode->destroyRenderingData(m_sceneMem);
    }
    
    bool Scene::intersect(const Ray &ray, const RaySegment &segment, SurfaceInteraction *si) const {
        float importances[2] = {m_surfaceAggregate->importance(), 0.0f};
        if (m_envSphere)
            importances[1] = m_envSphere->importance();
        
        if (m_surfaceAggregate->intersect(ray, segment, si)) {
            si->setLightProb(evaluateProbability(importances, 2, 0) * si->getLightProb());
            return true;
        }
        if (m_envSphere) {
            if (m_envSphere->intersect(ray, segment, si)) {
                si->setLightProb(evaluateProbability(importances, 2, 1) * si->getLightProb());
                return true;
            }
        }
        return false;
    }
    
    bool Scene::interact(const Ray &ray, const RaySegment &segment, const WavelengthSamples &wls, LightPathSampler &pathSampler, ArenaAllocator &mem,
                         Interaction** interact, SampledSpectrum* medThroughput, bool* singleWavelength) const {
        float importances[3] = {m_surfaceAggregate->importance(), m_mediumAggregate->importance(), 0.0f};
        if (m_envSphere)
            importances[2] = m_envSphere->importance();
        
        SurfaceInteraction si;
        bool hitSurface = m_surfaceAggregate->intersect(ray, segment, &si);
        MediumInteraction mi;
        bool hitMedium = m_mediumAggregate->interact(ray, RaySegment(segment.distMin, si.getDistance()), wls, pathSampler, &mi, medThroughput, singleWavelength);
        if (hitMedium) {
            mi.setLightProb(evaluateProbability(importances, 3, 1) * mi.getLightProb());
            *interact = mem.create<MediumInteraction>(mi);
            return true;
        }
        if (hitSurface) {
            mi.setLightProb(evaluateProbability(importances, 3, 0) * mi.getLightProb());
            *interact = mem.create<SurfaceInteraction>(si);
            return true;
        }
        if (m_envSphere) {
            if (m_envSphere->intersect(ray, segment, &si)) {
                mi.setLightProb(evaluateProbability(importances, 3, 2) * mi.getLightProb());
                *interact = mem.create<SurfaceInteraction>(si);
                return true;
            }
        }
        return false;
    }
    
    bool Scene::testVisibility(const SurfacePoint &shdP, const SurfacePoint &lightP, float time) const {
        SLRAssert(shdP.atInfinity() == false, "Shading point must be in finite region.");
        Ray ray;
        RaySegment segment;
        if (lightP.atInfinity()) {
            ray = Ray(shdP.getPosition(), normalize(lightP.getPosition() - Point3D::Zero), time);
            segment = RaySegment(Ray::Epsilon, FLT_MAX);
        }
        else {
            float dist = distance(lightP.getPosition(), shdP.getPosition());
            ray = Ray(shdP.getPosition(), (lightP.getPosition() - shdP.getPosition()) / dist, time);
            segment = RaySegment(Ray::Epsilon, dist * (1 - Ray::Epsilon));
        }
        SurfaceInteraction si;
        return !intersect(ray, segment, &si);
    }
    
    bool Scene::testVisibility(const InteractionPoint* shdP, const InteractionPoint* lightP, float time,
                               const WavelengthSamples &wls, LightPathSampler &pathSampler, SampledSpectrum* fractionalVisibility, bool* singleWavelength) const {
        SLRAssert(shdP->atInfinity() == false, "Shading point must be in finite region.");
        Ray ray;
        RaySegment segment;
        if (lightP->atInfinity()) {
            ray = Ray(shdP->getPosition(), normalize(lightP->getPosition() - Point3D::Zero), time);
            segment = RaySegment(Ray::Epsilon, FLT_MAX);
        }
        else {
            float dist = distance(lightP->getPosition(), shdP->getPosition());
            ray = Ray(shdP->getPosition(), (lightP->getPosition() - shdP->getPosition()) / dist, time);
            segment = RaySegment(Ray::Epsilon, dist * (1 - Ray::Epsilon));
        }
        *fractionalVisibility = SampledSpectrum::Zero;
        SurfaceInteraction si;
        if (m_surfaceAggregate->intersect(ray, segment, &si))
            return false;
        *fractionalVisibility = m_mediumAggregate->evaluateTransmittance(ray, segment, wls, pathSampler, singleWavelength);
        return true;
    }
    
    void Scene::selectSurfaceLight(float u, float time, SurfaceLight* light, float* prob) const {
        if (m_envSphere) {
            float sumImps = m_surfaceAggregate->importance() + m_envSphere->importance();
            float su = sumImps * u;
            if (su < m_surfaceAggregate->importance()) {
                u = u / (m_surfaceAggregate->importance() / sumImps);
                m_surfaceAggregate->selectLight(u, time, light, prob);
                *prob *= m_surfaceAggregate->importance() / sumImps;
            }
            else {
                u = (u - m_surfaceAggregate->importance()) / (m_envSphere->importance() / sumImps);
                m_envSphere->selectLight(u, time, light, prob);
                *prob *= m_envSphere->importance() / sumImps;
            }
        }
        else {
            m_surfaceAggregate->selectLight(u, time, light, prob);
        }
    }
    
    void Scene::selectLight(float u, float time, ArenaAllocator &mem, Light** light, float *prob) const {
        float importances[3] = {m_surfaceAggregate->importance(), m_mediumAggregate->importance(), 0.0f};
        if (m_envSphere)
            importances[2] = m_envSphere->importance();
        float prob1st;
        float sumImportances;
        uint32_t idx = sampleDiscrete(importances, 3, u, &prob1st, &sumImportances, &u);
        
        switch (idx) {
            case 0: {
                SurfaceLight* surfLight = mem.create<SurfaceLight>();
                m_surfaceAggregate->selectLight(u, time, surfLight, prob);
                *light = surfLight;
                break;
            }
            case 1: {
                VolumetricLight* volLight = mem.create<VolumetricLight>();
                m_mediumAggregate->selectLight(u, time, volLight, prob);
                *light = volLight;
                break;
            }
            case 2: {
                SurfaceLight* surfLight = mem.create<SurfaceLight>();
                m_envSphere->selectLight(u, time, surfLight, prob);
                *light = surfLight;
                break;
            }
            default:
                break;
        }
        *prob *= prob1st;
    }
}

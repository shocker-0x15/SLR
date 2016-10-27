//
//  Scene.cpp
//
//  Created by 渡部 心 on 2016/09/19.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "Scene.h"

namespace SLR {
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

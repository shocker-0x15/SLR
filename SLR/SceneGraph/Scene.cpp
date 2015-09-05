//
//  Scene.cpp
//
//  Created by 渡部 心 on 2015/07/11.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "Scene.h"
#include "nodes.h"
#include "../Memory/ArenaAllocator.h"
#include "../Core/cameras.h"
#include "../Core/SurfaceObject.h"
#include "../SceneGraph/InfiniteSphereNode.h"

Scene::Scene() : m_mainCamera(nullptr), m_envNode(nullptr), m_envSphere(nullptr) {
    m_rootNode = createShared<InternalNode>();
    m_rootNode->setName("root");
    m_rootNode->setTransform(createShared<StaticTransform>(Matrix4x4::Identity));
    
    m_mem = new ArenaAllocator();
};

Scene::~Scene() {
    if (m_sceneObj)
        delete m_sceneObj;
    if (m_mem)
        delete m_mem;
}

void Scene::build() {    
    RenderingData renderingData;
    m_rootNode->getRenderingData(*m_mem, nullptr, &renderingData);
    SLRAssert(renderingData.camera != nullptr, "Camera is not set.");
    
    // calculate camera transformation.
    std::vector<const Transform*> camTFs;
    InternalNode* parentNode = renderingData.camParent;
    while (parentNode) {
        camTFs.push_back(parentNode->getTransform().get());
        parentNode = parentNode->getParent();
    }
    ChainedTransform* chain = m_mem->create<ChainedTransform>(nullptr, camTFs[0]);
    ChainedTransform* prevChain = chain;
    ChainedTransform* head = chain;
    for (int i = 1; i < camTFs.size() - 1; ++i) {
        chain = m_mem->create<ChainedTransform>(nullptr, camTFs[i]);
        prevChain->setParent(chain);
        prevChain = chain;
    }
    prevChain->setParent(camTFs.back());
    
    m_mainCamera = renderingData.camera;
    m_mainCamera->setTransform(head->reduce(*m_mem));
    
    // calculate scene BVH.
    m_sceneObj = new SurfaceObjectAggregate(renderingData.surfObjs);
    if (m_envNode)
        m_envSphere = m_envNode->getSurfaceObject();
    
    // calculate world bounding sphere and store its radius and disc area.
    BoundingBox3D worldBounds = m_sceneObj->bounds();
    m_worldCenter = worldBounds.centroid();
    m_worldRadius = (worldBounds.maxP - m_worldCenter).length();
    m_worldDiscArea = M_PI * m_worldRadius * m_worldRadius;
}

bool Scene::intersect(Ray &ray, Intersection *isect) const {
    if (m_sceneObj->intersect(ray, isect))
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
    if (m_envNode) {
        float sumImps = m_sceneObj->importance() + m_envSphere->importance();
        float su = sumImps * u;
        if (su < m_sceneObj->importance()) {
            u = u / (m_sceneObj->importance() / sumImps);
            m_sceneObj->selectLight(u, light, prob);
            *prob *= m_sceneObj->importance() / sumImps;
        }
        else {
            u = (u - m_sceneObj->importance()) / (m_envSphere->importance() / sumImps);
            m_envSphere->selectLight(u, light, prob);
            *prob *= m_envSphere->importance() / sumImps;
        }
    }
    else {
        m_sceneObj->selectLight(u, light, prob);
    }
}

float Scene::evaluateProb(const Light &light) const {
    if (m_envNode) {
        float sumImps = m_sceneObj->importance() + m_envSphere->importance();
        if (light.top() == m_envSphere) {
            return m_envSphere->importance() / sumImps;
        }
        else {
            return m_sceneObj->importance()  / sumImps * m_sceneObj->evaluateProb(light);
        }
    }
    else {
        return m_sceneObj->evaluateProb(light);
    }
}

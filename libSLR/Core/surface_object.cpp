//
//  surface_object.cpp
//
//  Created by 渡部 心 on 2015/07/15.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "surface_object.h"

#include "../MemoryAllocators/ArenaAllocator.h"
#include "distributions.h"
#include "transform.h"
#include "light_path_sampler.h"
#include "textures.h"
#include "surface_material.h"
#include "../Accelerator/StandardBVH.h"
#include "../Accelerator/SBVH.h"
#include "../Accelerator/QBVH.h"
#include "../SurfaceShape/InfiniteSphereSurfaceShape.h"
#include "../BSDF/basic_bsdfs.h"
#include "../SurfaceMaterial/IBLEmitterSurfaceProperty.h"
#include "../Scene/Scene.h"

namespace SLR {
    SampledSpectrum SurfaceLight::sample(const LightPosQuery &query, const SurfaceLightPosSample &smp, SurfaceLightPosQueryResult* result) const {
        return m_obj->sample(m_appliedTransform, query, smp, result);
    }
    
    void SurfaceLight::sampleRay(const LightPosQuery &lightPosQuery, const SurfaceLightPosSample &lightPosSample,
                                 const EDFQuery &edfQuery, const EDFSample &edfSample,
                                 ArenaAllocator &mem,
                                 SurfaceLightPosQueryResult* lightPosResult, SampledSpectrum* Le0, EDF** edf,
                                 EDFQueryResult* edfResult, SampledSpectrum* Le1, Ray* ray, float* epsilon) const {
        m_obj->sampleRay(m_appliedTransform, lightPosQuery, lightPosSample, edfQuery, edfSample, mem, lightPosResult, Le0, edf, edfResult, Le1, ray, epsilon);
    }
    
    SampledSpectrum SurfaceLight::sample(const LightPosQuery &query, LightPathSampler &pathSampler, ArenaAllocator &mem, LightPosQueryResult** lpResult) const {
        *lpResult = mem.create<SurfaceLightPosQueryResult>();
        SampledSpectrum ret = sample(query, pathSampler.getSurfaceLightPosSample(), (SurfaceLightPosQueryResult*)*lpResult);
        return ret;
    }
    
    void SurfaceLight::sampleRay(const LightPosQuery &lightPosQuery, LightPathSampler &pathSampler, const EDFQuery &edfQuery, ArenaAllocator &mem,
                                 LightPosQueryResult** lightPosResult, SampledSpectrum* Le0, EDF** edf,
                                 EDFQueryResult* edfResult, SampledSpectrum* Le1, Ray* ray, float* epsilon) const {
        *lightPosResult = mem.create<SurfaceLightPosQueryResult>();
        sampleRay(lightPosQuery, pathSampler.getSurfaceLightPosSample(), edfQuery, pathSampler.getEDFSample(), mem,
                  (SurfaceLightPosQueryResult*)*lightPosResult, Le0, edf, edfResult, Le1, ray, epsilon);
    }
    
    
    
    bool SingleSurfaceObject::isEmitting() const {
        return m_material->isEmitting();
    }
    
    float SingleSurfaceObject::importance() const {
        return 1.0f;// TODO: consider a total power emitted from this object.
    }
    
    void SingleSurfaceObject::selectLight(float u, float time, SurfaceLight* light, float* prob) const {
        light->setObject(this);
        *prob = 1.0f;
    }
    
    SampledSpectrum SingleSurfaceObject::sample(const StaticTransform &transform,
                                                const LightPosQuery &query, const SurfaceLightPosSample &smp, SurfaceLightPosQueryResult* result) const {
        m_surface->sample(smp.uPos[0], smp.uPos[1], &result->surfPt, &result->areaPDF, &result->posType);
        result->surfPt.setObject(this);
        result->surfPt.applyTransform(transform);
        return m_material->emittance(result->surfPt, query.wls);
    }
    
    void SingleSurfaceObject::sampleRay(const StaticTransform &transform,
                                        const LightPosQuery &lightPosQuery, const SurfaceLightPosSample &lightPosSample,
                                        const EDFQuery &edfQuery, const EDFSample &edfSample,
                                        ArenaAllocator &mem,
                                        SurfaceLightPosQueryResult *lightPosResult, SampledSpectrum *Le0, EDF **edf,
                                        EDFQueryResult *edfResult, SampledSpectrum *Le1, Ray* ray, float* epsilon) const {
        // JP: 光源面上の位置を放射発散度とともにサンプルする。
        // EN: sample a position with emittance on the light's surface.
        *Le0 = sample(transform, lightPosQuery, lightPosSample, lightPosResult);
        *edf = lightPosResult->surfPt.createEDF(lightPosQuery.wls, mem);
        SLRAssert(!std::isnan(lightPosResult->areaPDF)/* && !std::isinf(lightResult)*/, "areaPDF: unexpected value detected: %f", lightPosResult->areaPDF);
        
        // JP: EDFから方向とその値をサンプルする。
        // EN: sample a direction from EDF with its value.
        *Le1 = (*edf)->sample(edfQuery, edfSample, edfResult);
        *ray = Ray(lightPosResult->surfPt.getPosition(), lightPosResult->surfPt.fromLocal(edfResult->dir_sn), lightPosQuery.time);
        *epsilon = Ray::Epsilon;
    }
    
    bool SingleSurfaceObject::intersectWithoutAlpha(const Ray &ray, const RaySegment &segment, SurfaceInteraction* si) const {
        if (!m_surface->intersectWithoutAlpha(ray, segment, si))
            return false;
        
        si->setObject(this);
        return true;
    }
    
    bool SingleSurfaceObject::intersect(const Ray &ray, const RaySegment &segment, LightPathSampler &pathSampler, SurfaceInteraction* si) const {
#ifdef DEBUG
        if (Accelerator::traceTraverse) {
            debugPrintf("%s%p, SingleSurfaceObject::intersect()\n",
                        Accelerator::traceTraversePrefix.c_str(), this);
            Accelerator::traceTraversePrefix += "  ";
        }
#endif
        if (!m_surface->intersect(ray, segment, pathSampler, si)) {
#ifdef DEBUG
            if (Accelerator::traceTraverse) {
                debugPrintf("%snot found\n", Accelerator::traceTraversePrefix.c_str());
                size_t newLength = Accelerator::traceTraversePrefix.length() - 2;
                Accelerator::traceTraversePrefix.resize(newLength);
            }
#endif
            return false;
        }
        si->setObject(this);
        si->setLightProb(isEmitting() ? 1.0f : 0.0f);
#ifdef DEBUG
        if (Accelerator::traceTraverse) {
            debugPrintf("%sfound: %g, %g\n",
                        Accelerator::traceTraversePrefix.c_str(), segment.distMax, si->getDistance());
            size_t newLength = Accelerator::traceTraversePrefix.length() - 2;
            Accelerator::traceTraversePrefix.resize(newLength);
        }
#endif
        return true;
    }
    
    float SingleSurfaceObject::testVisibility(const Ray &ray, const RaySegment &segment) const {
        return m_surface->testVisibility(ray, segment);
    }
    
    void SingleSurfaceObject::calculateSurfacePoint(const SurfaceInteraction &si, SurfacePoint* surfPt) const {
        m_surface->calculateSurfacePoint(si, surfPt);
        surfPt->setObject(this);
        surfPt->applyTransform(si.getAppliedTransform());
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
    
    
    
    void BumpSingleSurfaceObject::calculateSurfacePoint(const SurfaceInteraction &si, SurfacePoint *surfPt) const {
        SingleSurfaceObject::calculateSurfacePoint(si, surfPt);
        
        const ReferenceFrame &originalFrame = surfPt->getShadingFrame();
        
        // JP: 法線マップに従ってシェーディングフレームを変更する。
        // EN: perturb the shading frame according to the normal map.
        Vector3D nLocal = m_normalMap->evaluate(*surfPt);
        Vector3D tLocal = Vector3D::Ex - dot(nLocal, Vector3D::Ex) * nLocal;
        Vector3D bLocal = Vector3D::Ey - dot(nLocal, Vector3D::Ey) * nLocal;
        Vector3D t = normalize(originalFrame.fromLocal(tLocal));
        Vector3D b = normalize(originalFrame.fromLocal(bLocal));
        Vector3D n = normalize(originalFrame.fromLocal(nLocal));
        
        ReferenceFrame bumpFrame(t, b, n);
        surfPt->setShadingFrame(bumpFrame);
    }
    
    
    
    InfiniteSphereSurfaceObject::InfiniteSphereSurfaceObject(const Scene* scene, const IBLEmitterSurfaceProperty* emitter) :
    m_scene(scene) {
        m_surface = new InfiniteSphereSurfaceShape();
        m_material = new EmitterSurfaceMaterial(nullptr, emitter);
        m_dist = emitter->createIBLImportanceMap();
    }
    
    InfiniteSphereSurfaceObject::~InfiniteSphereSurfaceObject() {
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
    
    SampledSpectrum InfiniteSphereSurfaceObject::sample(const StaticTransform &transform,
                                                        const LightPosQuery &query, const SurfaceLightPosSample &smp, SurfaceLightPosQueryResult* result) const {
        // JP: テクスチャー空間からサンプルする。
        // EN: sample from texture space.
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
        result->surfPt.applyTransform(transform);
        // JP: テクスチャー空間中のPDFを面積に関するものに変換する。
        // EN: convert the PDF in texture space to one with respect to area.
        // The true value is: lim_{l to inf} uvPDF / (2 * M_PI * M_PI * std::sin(theta)) / l^2
        result->areaPDF = uvPDF / (2 * M_PI * M_PI * std::sin(theta));
        result->posType = DirectionType::LowFreq;
        return m_material->emittance(result->surfPt, query.wls);
    }
    
    void InfiniteSphereSurfaceObject::sampleRay(const StaticTransform &transform,
                                                const LightPosQuery &lightPosQuery, const SurfaceLightPosSample &lightPosSample,
                                                const EDFQuery &edfQuery, const EDFSample &edfSample,
                                                ArenaAllocator &mem,
                                                SurfaceLightPosQueryResult *lightPosResult, SampledSpectrum *Le0, EDF **edf,
                                                EDFQueryResult *edfResult, SampledSpectrum *Le1, Ray* ray, float* epsilon) const {
        // JP: 光源面上の位置を放射発散度とともにサンプルする。
        // EN: sample a position with emittance on the light's surface.
        *Le0 = sample(transform, lightPosQuery, lightPosSample, lightPosResult);
        *edf = lightPosResult->surfPt.createEDF(lightPosQuery.wls, mem);
        SLRAssert(!std::isnan(lightPosResult->areaPDF)/* && !std::isinf(lightResult)*/, "areaPDF: unexpected value detected: %f", lightPosResult->areaPDF);
        
        // JP: EDFから方向とその値をサンプルする。
        //     (無限に大きいの球上の)ある点からサンプルされるあらゆる方向は平行にならなければならないが、
        //     一方でシーン中のあらゆる点に到達できる必要がある。
        //     そのためレイの原点を修正する必要がある。
        // EN: sample a direction from EDF with its value.
        //     All the sampled directions for a certain position (on the infinite sphere) must be parallel,
        //     but be able to reach any position in the scene.
        //     Therefore, it requires modification to ray's origin.
        *Le1 = (*edf)->sample(edfQuery, edfSample, edfResult);
        Vector3D vx, vy;
        Vector3D vz = lightPosResult->surfPt.fromLocal(edfResult->dir_sn);
        vz.makeCoordinateSystem(&vx, &vy);
        float dx, dy;
        concentricSampleDisk(edfSample.uDir[0], edfSample.uDir[1], &dx, &dy);
        
        // JP: シーンのバウンディングスフィアの大円中からサンプルする。
        // EN: sample within a great circle of the bounding sphere of the scene.
        float worldRadius = m_scene->getWorldRadius();
        Point3D org = m_scene->getWorldCenter() + 1.1f * worldRadius * lightPosResult->surfPt.getPosition() + worldRadius * (dx * vx + dy * vy);
        *ray = Ray(org, vz, lightPosQuery.time);
        *epsilon = 0.0f;
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
    
    void TransformedSurfaceObject::selectLight(float u, float time, SurfaceLight* light, float* prob) const {
        m_surfObj->selectLight(u, time, light, prob);
        StaticTransform tf;
        m_transform->sample(time, &tf);
        light->applyTransformFromLeft(tf);
    }
    
    bool TransformedSurfaceObject::contains(const Point3D &p, float time) const {
        StaticTransform sampledTF;
        m_transform->sample(time, &sampledTF);
        Point3D localP = invert(sampledTF) * p;
        return m_surfObj->contains(localP, time);
    }
    
    bool TransformedSurfaceObject::intersectWithoutAlpha(const Ray &ray, const RaySegment &segment, SLR::SurfaceInteraction *si) const {
        Ray localRay;
        StaticTransform sampledTF;
        m_transform->sample(ray.time, &sampledTF);
        localRay = invert(sampledTF) * ray;
        
        if (!m_surfObj->intersectWithoutAlpha(localRay, segment, si))
            return false;
            
        si->applyTransformFromLeft(sampledTF);
        return true;
    }
    
    bool TransformedSurfaceObject::intersect(const Ray &ray, const RaySegment &segment, LightPathSampler &pathSampler, SurfaceInteraction* si) const {
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
        if (!m_surfObj->intersect(localRay, segment, pathSampler, si)) {
#ifdef DEBUG
            if (Accelerator::traceTraverse) {
                debugPrintf("%snot found\n", Accelerator::traceTraversePrefix.c_str());
                size_t newLength = Accelerator::traceTraversePrefix.length() - 2;
                Accelerator::traceTraversePrefix.resize(newLength);
            }
#endif
            return false;
        }
        si->applyTransformFromLeft(sampledTF);
#ifdef DEBUG
        if (Accelerator::traceTraverse) {
            debugPrintf("%sfound: %g\n",
                        Accelerator::traceTraversePrefix.c_str(), segment.distMax);
            size_t newLength = Accelerator::traceTraversePrefix.length() - 2;
            Accelerator::traceTraversePrefix.resize(newLength);
        }
#endif
        return true;
    }
    
    float TransformedSurfaceObject::testVisibility(const Ray &ray, const RaySegment &segment) const {
        Ray localRay;
        StaticTransform sampledTF;
        m_transform->sample(ray.time, &sampledTF);
        localRay = invert(sampledTF) * ray;
        return m_surfObj->testVisibility(ray, segment);
    }
    
    
    
    SurfaceObjectAggregate::SurfaceObjectAggregate(std::vector<SurfaceObject*> &objs) {
//        SBVH sbvh(objs);
//        m_accelerator = new QBVH(sbvh);
        m_accelerator = new SBVH(objs);
//        m_accelerator = new StandardBVH(objs, StandardBVH::Partitioning::BinnedSAH);
        
        std::vector<uint32_t> lightIndices;
        std::vector<float> lightImportances;
        for (int i = 0; i < objs.size(); ++i) {
            const SurfaceObject* obj = objs[i];
            if (obj->isEmitting()) {
                lightIndices.push_back(i);
                lightImportances.push_back(obj->importance());
            }
        }
        
        m_numLights = (uint32_t)lightImportances.size();
        m_lightList = new const SurfaceObject*[m_numLights];
        m_lightDist1D = new DiscreteDistribution1D(lightImportances);
        
        for (int i = 0; i < m_numLights; ++i) {
            uint32_t objIdx = lightIndices[i];
            const SurfaceObject* light = objs[objIdx];
            m_lightList[i] = light;
            m_objToLightMap[light] = i;
        }
    }
    
    SurfaceObjectAggregate::~SurfaceObjectAggregate() {
        delete m_accelerator;
        
        delete m_lightDist1D;
        delete[] m_lightList;
    };
    
    BoundingBox3D SurfaceObjectAggregate::bounds() const {
        return m_accelerator->bounds();
    }
    
    bool SurfaceObjectAggregate::isEmitting() const {
        return m_numLights > 0;
    }
    
    float SurfaceObjectAggregate::importance() const {
        return m_lightDist1D->integral();
    }
    
    void SurfaceObjectAggregate::selectLight(float u, float time, SurfaceLight* light, float* prob) const {
        uint32_t lIdx = m_lightDist1D->sample(u, prob, &u);
        const SurfaceObject* obj = m_lightList[lIdx];
        float cProb;
        obj->selectLight(u, time, light, &cProb);
        *prob *= cProb;
    }
    
    float SurfaceObjectAggregate::costForIntersect() const {
        return m_accelerator->costForIntersect();
    }
    
    bool SurfaceObjectAggregate::contains(const Point3D &p, float time) const {
        Ray probeRay(p, Vector3D::Ex, time);
        SurfaceInteraction si;
        bool hit = m_accelerator->intersectWithoutAlpha(probeRay, RaySegment(), &si);
        if (!hit)
            return false;
        return dot(si.getGeometricNormal(), probeRay.dir) >= 0.0f;
    }
    
    bool SurfaceObjectAggregate::intersectWithoutAlpha(const Ray &ray, const RaySegment &segment, SurfaceInteraction* si) const {
        return m_accelerator->intersectWithoutAlpha(ray, segment, si);
    }
    
    bool SurfaceObjectAggregate::intersect(const Ray &ray, const RaySegment &segment, LightPathSampler &pathSampler, SurfaceInteraction* si) const {
#ifdef DEBUG
        if (Accelerator::traceTraverse) {
            debugPrintf("%s%p, SurfaceObjectAggregate::intersect()\n",
                        Accelerator::traceTraversePrefix.c_str(), this);
            Accelerator::traceTraversePrefix += "  ";
        }
#endif
        const SurfaceObject* hitObj;
        if (!m_accelerator->intersect(ray, segment, pathSampler, si, &hitObj)) {
#ifdef DEBUG
            if (Accelerator::traceTraverse) {
                debugPrintf("%snot found\n", Accelerator::traceTraversePrefix.c_str());
                size_t newLength = Accelerator::traceTraversePrefix.length() - 2;
                Accelerator::traceTraversePrefix.resize(newLength);
            }
#endif
            return false;
        }
        if (m_objToLightMap.count(hitObj) > 0) {
            uint32_t lightIdx = m_objToLightMap.at(hitObj);
            si->setLightProb(m_lightDist1D->evaluatePMF(lightIdx) * si->getLightProb());
        }
#ifdef DEBUG
        if (Accelerator::traceTraverse) {
            debugPrintf("%sfound: %g\n",
                        Accelerator::traceTraversePrefix.c_str(), segment.distMax);
            size_t newLength = Accelerator::traceTraversePrefix.length() - 2;
            Accelerator::traceTraversePrefix.resize(newLength);
        }
#endif
        return true;
    }
    
    float SurfaceObjectAggregate::testVisibility(const Ray &ray, const RaySegment &segment) const {
        return m_accelerator->testVisibility(ray, segment);
    }
}

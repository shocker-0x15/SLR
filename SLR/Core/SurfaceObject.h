//
//  SurfaceObject.h
//
//  Created by 渡部 心 on 2015/07/15.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__SurfaceObject__
#define __SLR__SurfaceObject__

#include "../defines.h"
#include "../references.h"
#include "geometry.h"

struct LightPosQuery {
    float time;
    WavelengthSamples wls;
    LightPosQuery(float t, const WavelengthSamples &lambdas) : time(t), wls(lambdas) { };
};

struct LightPosSample {
    float uPos[2];
    LightPosSample(float up0, float up1) : uPos{up0, up1} { }
};

struct LightPosQueryResult {
    SurfacePoint surfPt;
    float areaPDF;
    bool isDeltaPos;
};

class Light {
    mutable std::stack<const SurfaceObject*> m_hierarchy;
public:
    Light() { };
    Light(const std::stack<const SurfaceObject*> &hierarchy) : m_hierarchy(hierarchy) { };
    
    void push(const SurfaceObject* obj) const { m_hierarchy.push(obj); };
    void pop() const { m_hierarchy.pop(); };
    const SurfaceObject* top() const { return m_hierarchy.top(); };
    
    SampledSpectrum sample(const LightPosQuery &query, const LightPosSample &smp, LightPosQueryResult* result) const;
};

class SurfaceObject {
public:
    SurfaceObject() { };
    virtual ~SurfaceObject() { };
    
    virtual BoundingBox3D bounds() const = 0;
    virtual bool intersect(Ray &ray, Intersection* isect) const = 0;
    virtual void getSurfacePoint(const Intersection &isect, SurfacePoint* surfPt) const = 0;

    virtual bool isEmitting() const = 0;
    virtual float importance() const = 0;
    virtual void selectLight(float u, Light* light, float* prob) const = 0;
    virtual float evaluateProb(const Light &light) const = 0;
    
    virtual SampledSpectrum sample(const Light &light, const LightPosQuery &query, const LightPosSample &smp, LightPosQueryResult* result) const = 0;
    
    bool intersect(Ray &ray, SurfacePoint* surfPt) const;
    bool testVisiblility(const SurfacePoint &shdP, const SurfacePoint &lightP, float time) const;
};

class SingleSurfaceObject : public SurfaceObject {
protected:
    const Surface* m_surface;
    const SurfaceMaterial* m_material;
public:
    SingleSurfaceObject() { };
    SingleSurfaceObject(const Surface* surf, const SurfaceMaterial* mat) : m_surface(surf), m_material(mat) { };
    virtual ~SingleSurfaceObject() { };
    
    BoundingBox3D bounds() const override;
    bool intersect(Ray &ray, Intersection* isect) const override;
    void getSurfacePoint(const Intersection &isect, SurfacePoint* surfPt) const override;
    
    bool isEmitting() const override;
    float importance() const override;
    void selectLight(float u, Light* light, float* prob) const override;
    float evaluateProb(const Light &light) const override;
    
    SampledSpectrum sample(const Light &light, const LightPosQuery &query, const LightPosSample &smp, LightPosQueryResult* result) const override;
    
    virtual BSDF* createBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const;
    
    virtual float evaluateAreaPDF(const SurfacePoint& surfPt) const;
    virtual SampledSpectrum emittance(const SurfacePoint& surfPt, const WavelengthSamples &wls) const;
    virtual EDF* createEDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const;
};

class BumpSingleSurfaceObject : public SingleSurfaceObject {
    const Normal3DTexture* m_normalMap;
public:
    BumpSingleSurfaceObject(const Surface* surf, const SurfaceMaterial* mat, const Normal3DTexture* normalMap) :
    SingleSurfaceObject(surf, mat), m_normalMap(normalMap) { };
    
    void getSurfacePoint(const Intersection &isect, SurfacePoint* surfPt) const override;
};

class InfiniteSphereSurfaceObject : public SingleSurfaceObject {
    const RegularConstantContinuous2D* m_dist;
public:
    InfiniteSphereSurfaceObject(const Surface* surf, const SurfaceMaterial* mat, const RegularConstantContinuous2D* dist) :
    SingleSurfaceObject(surf, mat), m_dist(dist) { };
    ~InfiniteSphereSurfaceObject() { };
    
    bool isEmitting() const override;
    float importance() const override;
    
    SampledSpectrum sample(const Light &light, const LightPosQuery &query, const LightPosSample &smp, LightPosQueryResult* result) const override;
    
    BSDF* createBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const override;
    
    float evaluateAreaPDF(const SurfacePoint& surfPt) const override;
};

class SurfaceObjectAggregate : public SurfaceObject {
    BBVH* m_accelerator;
    const SurfaceObject** m_lightList;
    RegularConstantDiscrete1D* m_lightDist1D;
    std::map<const SurfaceObject*, uint32_t> m_revMap;
public:
    SurfaceObjectAggregate(std::vector<SurfaceObject*> &objs);
    ~SurfaceObjectAggregate();
    
    BoundingBox3D bounds() const override;
    bool intersect(Ray &ray, Intersection* isect) const override;
    void getSurfacePoint(const Intersection &isect, SurfacePoint* surfPt) const override;
    
    bool isEmitting() const override;
    float importance() const override;
    void selectLight(float u, Light* light, float* prob) const override;
    float evaluateProb(const Light &light) const override;
    
    SampledSpectrum sample(const Light &light, const LightPosQuery &query, const LightPosSample &smp, LightPosQueryResult* result) const override;
};

class TransformedSurfaceObject : public SurfaceObject {
    const SurfaceObject* m_surfObj;
    const Transform* m_transform;
    friend class Light;
public:
    TransformedSurfaceObject(const SurfaceObject* surfObj, const Transform* transform) : m_surfObj(surfObj), m_transform(transform) { };
    
    BoundingBox3D bounds() const override;
    bool intersect(Ray &ray, Intersection* isect) const override;
    void getSurfacePoint(const Intersection &isect, SurfacePoint* surfPt) const override;
    
    bool isEmitting() const override;
    float importance() const override;
    void selectLight(float u, Light* light, float* prob) const override;
    float evaluateProb(const Light &light) const override;
    
    SampledSpectrum sample(const Light &light, const LightPosQuery &query, const LightPosSample &smp, LightPosQueryResult* result) const override;
    
    void setTransform(const Transform* t) { m_transform = t; };
};

#endif

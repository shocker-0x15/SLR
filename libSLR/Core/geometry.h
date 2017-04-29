//
//  geometry.h
//
//  Created by 渡部 心 on 2015/05/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_geometry__
#define __SLR_geometry__

#include "../defines.h"
#include "../declarations.h"
#include "../BasicTypes/rgb_types.h"
#include "../BasicTypes/spectrum_types.h"
#include "../BasicTypes/TexCoord2D.h"
#include "transform.h"

namespace SLR {
    class SLR_API Interaction {
        friend class SurfacePoint;
        friend class MediumPoint;
        
        float m_time;
        float m_dist;
        Point3D m_p;
        StaticTransform m_appliedTransform;
        float m_lightProb;
    public:
        Interaction(float time, float dist, const Point3D &p) :
        m_time(time), m_dist(dist), m_p(p)
        {}
        
        float getTime() const { return m_time; }
        float getDistance() const { return m_dist; }
        void applyTransformFromLeft(const StaticTransform &transform) {
            m_appliedTransform = transform * m_appliedTransform;
        }
        const StaticTransform &getAppliedTransform() const { return m_appliedTransform; }
        void setLightProb(float prob) { m_lightProb = prob; }
        float getLightProb() const { return m_lightProb; }
        
        virtual InteractionPoint* createInteractionPoint(ArenaAllocator &mem) const = 0;
    };
    
    
    
    class SLR_API SurfaceInteraction : public Interaction {
        friend class SurfacePoint;
        
        const SingleSurfaceObject* m_obj;
        Normal3D m_gNormal;
        float m_u, m_v;
        TexCoord2D m_texCoord;
    public:
        SurfaceInteraction() : Interaction(0.0f, INFINITY, Point3D::Zero)
        {}
        SurfaceInteraction(float time, float dist, const Point3D &p,
                           const Normal3D &gNormal, float u, float v, const TexCoord2D &texCoord) :
        Interaction(time, dist, p), m_gNormal(gNormal), m_u(u), m_v(v), m_texCoord(texCoord)
        {}
        
        void setObject(const SingleSurfaceObject* obj) { m_obj = obj; }
        
        const Normal3D &getGeometricNormal() const { return m_gNormal; }
        void getSurfaceParameter(float* u, float* v) const {
            *u = m_u;
            *v = m_v;
        }
        
        void calculateSurfacePoint(SurfacePoint* surfPt) const;
        
        InteractionPoint* createInteractionPoint(ArenaAllocator &mem) const override;
    };
    
    
    
    class SLR_API MediumInteraction : public Interaction {
        friend class MediumPoint;
        
        const SingleMediumObject* m_obj;
        Vector3D m_dirIn;
        float m_u, m_v, m_w;
    public:
        MediumInteraction() : Interaction(0.0f, INFINITY, Point3D::Zero)
        {}
        MediumInteraction(float time, float dist, const Point3D &p, const Vector3D &dirIn, float u, float v, float w) :
        Interaction(time, dist, p), m_dirIn(dirIn), m_u(u), m_v(v), m_w(w)
        {}
        
        void setObject(const SingleMediumObject* obj) { m_obj = obj; }
        
        Vector3D getIncomingDirection() const { return m_dirIn; }
        void getMediumParameter(float* u, float* v, float* w) const {
            *u = m_u;
            *v = m_v;
            *w = m_w;
        }
        
        void calculateMediumPoint(MediumPoint* medPt) const;
        
        InteractionPoint* createInteractionPoint(ArenaAllocator &mem) const override;
    };
    
    
    
    struct SLR_API ReferenceFrame {
        Vector3D x, y, z;
        
        ReferenceFrame() { }
        ReferenceFrame(const Vector3D &t, const Vector3D &b, const Vector3D &n) : x(t), y(b), z(n) { }
        ReferenceFrame(const Vector3D &zz) : z(zz) {
            z.makeCoordinateSystem(&x, &y);
        }
        
        Vector3D toLocal(const Vector3D &v) const { return Vector3D(dot(x, v), dot(y, v), dot(z, v)); }
        Vector3D fromLocal(const Vector3D &v) const {
            // assume orthonormal basis
            return Vector3D(dot(Vector3D(x.x, y.x, z.x), v),
                            dot(Vector3D(x.y, y.y, z.y), v),
                            dot(Vector3D(x.z, y.z, z.z), v));
        }
    };
    
    
    
    class SLR_API InteractionPoint {
    protected:
        Point3D m_p;
        bool m_atInfinity;
        ReferenceFrame m_shadingFrame;
    public:
        InteractionPoint() {}
        InteractionPoint(const Point3D &p, bool atInfinity, const ReferenceFrame &shadingFrame) :
        m_p(p), m_atInfinity(atInfinity), m_shadingFrame(shadingFrame) {}
        
        const Point3D &getPosition() const { return m_p; }
        bool atInfinity() const { return m_atInfinity; }
        const ReferenceFrame &getShadingFrame() const { return m_shadingFrame; }
        void setShadingFrame(const ReferenceFrame &shadingFrame) {
            m_shadingFrame = shadingFrame;
        }
        
        float getSquaredDistance(const Point3D &shadingPoint) const {
            return m_atInfinity ? 1.0f : sqDistance(m_p, shadingPoint);
        }
        Vector3D getDirectionFrom(const Point3D &shadingPoint, float* dist2) const {
            if (m_atInfinity) {
                *dist2 = 1.0f;
                return normalize(m_p - Point3D::Zero);
            }
            else {
                Vector3D ret(m_p - shadingPoint);
                *dist2 = ret.sqLength();
                return ret / std::sqrt(*dist2);
            }
        }
        Vector3D toLocal(const Vector3D &vecWorld) const { return m_shadingFrame.toLocal(vecWorld); }
        Vector3D fromLocal(const Vector3D &vecLocal) const { return m_shadingFrame.fromLocal(vecLocal); }
        
        virtual bool isEmitting() const = 0;
        virtual SampledSpectrum emittance(const WavelengthSamples &wls) const = 0;
        virtual float evaluateSpatialPDF() const = 0;
        virtual EDF* createEDF(const WavelengthSamples &wls, ArenaAllocator &mem) const = 0;
        
        virtual ABDFQuery* createABDFQuery(const Vector3D &dirLocal, int16_t selectedWL, DirectionType dirType, bool reqRev, bool adjoint, ArenaAllocator &mem) const = 0;
        virtual AbstractBDF* createAbstractBDF(const WavelengthSamples &wls, ArenaAllocator &mem) const = 0;
        virtual SampledSpectrum evaluateExtinctionCoefficient(const WavelengthSamples &wls) const = 0;
        virtual float calcCosTerm(const Vector3D &vecWorld) const = 0;
        
        virtual void applyTransform(const StaticTransform &transform);
    };
    
    inline float squaredDistance(const InteractionPoint* p0, const InteractionPoint* p1) {
        return (p0->atInfinity() || p1->atInfinity()) ? 1.0f : sqDistance(p0->getPosition(), p1->getPosition());
    }
    
    
    
    class SLR_API SurfacePoint : public InteractionPoint {
        Normal3D m_gNormal;
        float m_u, m_v;
        TexCoord2D m_texCoord;
        Vector3D m_texCoord0Dir;
        const SingleSurfaceObject* m_obj;
    public:
        SurfacePoint() { }
        SurfacePoint(const Point3D &p, bool atInfinity, const ReferenceFrame &shadingFrame,
                     const Normal3D &gNormal, float u, float v, const TexCoord2D &texCoord, const Vector3D &texCoord0Dir) :
        InteractionPoint(p, atInfinity, shadingFrame),
        m_gNormal(gNormal), m_u(u), m_v(v), m_texCoord(texCoord), m_texCoord0Dir(texCoord0Dir) { }
        SurfacePoint(const SurfaceInteraction &si,
                     bool atInfinity, const ReferenceFrame &shadingFrame,
                     const Vector3D &texCoord0Dir) :
        InteractionPoint(si.m_p, atInfinity, shadingFrame),
        m_gNormal(si.m_gNormal), m_u(si.m_u), m_v(si.m_v), m_texCoord(si.m_texCoord), m_texCoord0Dir(texCoord0Dir) { }
        
        void setObject(const SingleSurfaceObject* obj) { m_obj = obj; }
        
        const Normal3D &getGeometricNormal() const { return m_gNormal; }
        void getSurfaceParameter(float* u, float* v) const {
            *u = m_u;
            *v = m_v;
        }
        const TexCoord2D &getTextureCoordinate() const { return m_texCoord; }
        const void setTextureCoordinate(const TexCoord2D &texCoord) { m_texCoord = texCoord; }
        
        Normal3D getLocalGeometricNormal() const {
            return m_shadingFrame.toLocal(m_gNormal);
        }
        
        float evaluateAreaPDF() const;
        BSDF* createBSDF(const WavelengthSamples &wls, ArenaAllocator &mem) const;
        
        //----------------------------------------------------------------
        // InteractionPoint's methods
        
        bool isEmitting() const override;
        SampledSpectrum emittance(const WavelengthSamples &wls) const override;
        float evaluateSpatialPDF() const override {
            return evaluateAreaPDF();
        }
        EDF* createEDF(const WavelengthSamples &wls, ArenaAllocator &mem) const override;
        
        ABDFQuery* createABDFQuery(const Vector3D &dirLocal, int16_t selectedWL, DirectionType filter, bool reqRev, bool adjoint, ArenaAllocator &mem) const override;
        AbstractBDF* createAbstractBDF(const WavelengthSamples &wls, ArenaAllocator &mem) const override {
            return (AbstractBDF*)createBSDF(wls, mem);
        }
        SampledSpectrum evaluateExtinctionCoefficient(const WavelengthSamples &wls) const override {
            return SampledSpectrum::One;
        }
        float calcCosTerm(const Vector3D &vecWorld) const override {
            return absDot(vecWorld, m_gNormal);
        }
        
        void applyTransform(const StaticTransform &transform) override;
        
        // END: InteractionPoint's methods
        //----------------------------------------------------------------
    };
    
    inline float squaredDistance(const SurfacePoint &p0, const SurfacePoint &p1) {
        return (p0.atInfinity() || p1.atInfinity()) ? 1.0f : sqDistance(p0.getPosition(), p1.getPosition());
    }
    
    
    
    class SLR_API MediumPoint : public InteractionPoint {
        float m_u, m_v, m_w;
        const SingleMediumObject* m_obj;
    public:
        MediumPoint() { }
        MediumPoint(const MediumInteraction &mi,
                    bool atInfinity, const ReferenceFrame &shadingFrame) :
        InteractionPoint(mi.m_p, atInfinity, shadingFrame), m_u(mi.m_u), m_v(mi.m_v), m_w(mi.m_w) { }
        
        void setObject(const SingleMediumObject* obj) { m_obj = obj; }
        
        void getMediumParameter(float* u, float* v, float* w) const {
            *u = m_u;
            *v = m_v;
            *w = m_w;
        }
        
        float evaluateVolumePDF() const;
        BSDF* createPhaseFunction(const WavelengthSamples &wls, ArenaAllocator &mem) const;
        
        //----------------------------------------------------------------
        // InteractionPoint's methods
        
        bool isEmitting() const override;
        SampledSpectrum emittance(const WavelengthSamples &wls) const override;
        float evaluateSpatialPDF() const override {
            return evaluateVolumePDF();
        }
        EDF* createEDF(const WavelengthSamples &wls, ArenaAllocator &mem) const override;
        
        ABDFQuery* createABDFQuery(const Vector3D &dirLocal, int16_t selectedWL, DirectionType filter, bool reqRev, bool adjoint, ArenaAllocator &mem) const override;
        AbstractBDF* createAbstractBDF(const WavelengthSamples &wls, ArenaAllocator &mem) const override;
        SampledSpectrum evaluateExtinctionCoefficient(const WavelengthSamples &wls) const override;
        float calcCosTerm(const Vector3D &vecWorld) const override {
            return 1.0f;
        }
        
        void applyTransform(const StaticTransform &transform) override;
        
        // END: InteractionPoint's methods
        //----------------------------------------------------------------
    };
    
    
    
    class SLR_API SurfaceShape {
    public:
        virtual ~SurfaceShape() { }
        
        virtual float costForIntersect() const = 0;
        virtual BoundingBox3D bounds() const = 0;
        virtual BoundingBox3D choppedBounds(BoundingBox3D::Axis chopAxis, float minChopPos, float maxChopPos) const {
            BoundingBox3D baseBBox = bounds();
            if (maxChopPos < baseBBox.minP[chopAxis])
                return BoundingBox3D();
            if (minChopPos > baseBBox.maxP[chopAxis])
                return BoundingBox3D();
            if (minChopPos < baseBBox.minP[chopAxis] && maxChopPos > baseBBox.maxP[chopAxis])
                return baseBBox;
            BoundingBox3D ret = baseBBox;
            ret.minP[chopAxis] = std::max(minChopPos, ret.minP[chopAxis]);
            ret.maxP[chopAxis] = std::min(maxChopPos, ret.maxP[chopAxis]);
            return ret;
        }
        virtual void splitBounds(BoundingBox3D::Axis splitAxis, float splitPos, BoundingBox3D* bbox0, BoundingBox3D* bbox1) const {
            BoundingBox3D baseBBox = bounds();
            if (splitPos < baseBBox.minP[splitAxis]) {
                *bbox0 = BoundingBox3D();
                *bbox1 = baseBBox;
                return;
            }
            if (splitPos > baseBBox.maxP[splitAxis]) {
                *bbox0 = baseBBox;
                *bbox1 = BoundingBox3D();
                return;
            }
            *bbox0 = baseBBox;
            bbox0->maxP[splitAxis] = std::min(bbox0->maxP[splitAxis], splitPos);
            *bbox1 = baseBBox;
            bbox1->minP[splitAxis] = std::max(bbox1->minP[splitAxis], splitPos);
        }
        virtual bool preTransformed() const = 0;
        virtual bool intersect(const Ray &ray, const RaySegment &segment, SurfaceInteraction* si) const = 0;
        virtual void calculateSurfacePoint(const SurfaceInteraction &si, SurfacePoint* surfPt) const = 0;
        virtual float area() const = 0;
        virtual void sample(float u0, float u1, SurfacePoint* surfPt, float* areaPDF, DirectionType* posType) const = 0;
        virtual float evaluateAreaPDF(const SurfacePoint& surfPt) const = 0;
    };
    
    
    
    class SLR_API MediumDistribution {
    protected:
    public:
        MediumDistribution() { }
        virtual ~MediumDistribution() { }

        virtual float majorantExtinctionCoefficientAtWavelength(float wl) const = 0;
        
        virtual bool subdivide(Allocator* mem, MediumDistribution** fragments, uint32_t* numFragments) const = 0;
        
        virtual BoundingBox3D bounds() const = 0;
        virtual bool contains(const Point3D &p) const = 0;
        virtual bool intersectBoundary(const Ray &ray, const RaySegment &segment, float* distToBoundary, bool* enter) const = 0;
        virtual bool interact(const Ray &ray, const RaySegment &segment, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                              MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const = 0;
        virtual SampledSpectrum evaluateTransmittance(const Ray &ray, const RaySegment &segment, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                                                      bool* singleWavelength) const = 0;
        virtual void calculateMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const = 0;
        virtual SampledSpectrum evaluateExtinctionCoefficient(const Point3D &param, const WavelengthSamples &wls) const = 0;
        virtual SampledSpectrum evaluateAlbedo(const Point3D &param, const WavelengthSamples &wls) const = 0;
        virtual float volume() const = 0;
        virtual void sample(float u0, float u1, float u2, MediumPoint* medPt, float* volumePDF) const = 0;
        virtual float evaluateVolumePDF(const MediumPoint& medPt) const = 0;
    };
}

#endif /* __SLR_geometry__ */

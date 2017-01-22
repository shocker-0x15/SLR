//
//  geometry.h
//
//  Created by 渡部 心 on 2015/05/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef SLR_geometry_h
#define SLR_geometry_h

#include "../defines.h"
#include "../references.h"
#include "../BasicTypes/Vector3.h"
#include "../BasicTypes/Vector4.h"
#include "../BasicTypes/Normal3.h"
#include "../BasicTypes/Point3.h"
#include "../BasicTypes/Matrix3x3.h"
#include "../BasicTypes/Matrix4x4.h"
#include "../BasicTypes/Quaternion.h"
#include "../BasicTypes/TexCoord2.h"
#include "../BasicTypes/RGBTypes.h"
#include "../BasicTypes/SpectrumTypes.h"

namespace SLR {
    struct SLR_API Ray {
        static const float Epsilon;
        
        Point3D org;
        Vector3D dir;
        float distMin, distMax;
        float time;
        
        Ray() { }
        Ray(const Point3D &o, const Vector3D &d, float t, float dMin = 0.0f, float dMax = INFINITY) :
        org(o), dir(d), distMin(dMin), distMax(dMax), time(t) { }
    };
    
    
    
    struct SLR_API BoundingBox3D {
        enum Axis : uint8_t {
            Axis_X = 0,
            Axis_Y,
            Axis_Z,
        };
        
        Point3D minP, maxP;
        
        BoundingBox3D() {
            minP.x = minP.y = minP.z = INFINITY;
            maxP.x = maxP.y = maxP.z = -INFINITY;
        }
        BoundingBox3D(const Point3D &p) : minP(p), maxP(p) { }
        BoundingBox3D(const Point3D &pmin, const Point3D &pmax) : minP(pmin), maxP(pmax) { }
        
        Point3D centroid() const {
            return (minP + maxP) * 0.5f;
        }
        
        float surfaceArea() const {
            Vector3D d = maxP - minP;
            return 2 * (d.x * d.y + d.y * d.z + d.z * d.x);
        }
        
        float volume() const {
            Vector3D d = maxP - minP;
            return d.x * d.y * d.z;
        }
        
        Point3D corner(uint32_t c) const {
            SLRAssert(c >= 0 && c < 8, "\"c\" is out of range [0, 8]");
            const size_t offset = sizeof(Point3D);
            return Point3D(*(&minP.x + offset * (c & 0x01)),
                           *(&minP.y + offset * (c & 0x02)),
                           *(&minP.z + offset * (c & 0x04)));
        }
        
        float centerOfAxis(Axis axis) const {
            return (minP[axis] + maxP[axis]) * 0.5f;
        }
        
        float width(Axis axis) const {
            return maxP[axis] - minP[axis];
        }
        
        Axis widestAxis() const {
            Vector3D d = maxP - minP;
            if (d.x > d.y && d.x > d.z)
                return Axis_X;
            else if (d.y > d.z)
                return Axis_Y;
            else
                return Axis_Z;
        }
        
        bool isValid() const {
            Vector3D d = maxP - minP;
            return d.x >= 0 && d.y >= 0 && d.z >= 0;
        }
        
        BoundingBox3D &unify(const Point3D &p) {
            minP = min(minP, p);
            maxP = max(maxP, p);
            return *this;
        }
        
        BoundingBox3D unify(const BoundingBox3D &b) {
            minP = min(minP, b.minP);
            maxP = max(maxP, b.maxP);
            return *this;
        }
        
        bool contains(const Point3D &p) const {
            return ((p.x >= minP.x && p.x < maxP.x) &&
                    (p.y >= minP.y && p.y < maxP.y) &&
                    (p.z >= minP.z && p.z < maxP.z));
        }
        
        // check intersection between ray segment and bounding volume.
        bool intersect(const Ray &r) const {
            float dist0 = r.distMin, dist1 = r.distMax;
            Vector3D invRayDir = r.dir.reciprocal();
            Vector3D tNear = (minP - r.org) * invRayDir;
            Vector3D tFar = (maxP - r.org) * invRayDir;
            for (int i = 0; i < 3; ++i) {
                if (tNear[i] > tFar[i])
                    std::swap(tNear[i], tFar[i]);
                dist0 = tNear[i] > dist0 ? tNear[i] : dist0;
                dist1 = tFar[i] < dist1 ? tFar[i] : dist1;
                if (dist0 > dist1)
                    return false;
            }
            return true;
        }
        
        // check intersection between ray segment and boundary.
        bool intersectBoundary(const Ray &r, float* distToBoundary, bool* enter) const {
            float dist0 = -INFINITY, dist1 = INFINITY;
            Vector3D invRayDir = r.dir.reciprocal();
            Vector3D tNear = (minP - r.org) * invRayDir;
            Vector3D tFar = (maxP - r.org) * invRayDir;
            for (int i = 0; i < 3; ++i) {
                if (tNear[i] > tFar[i])
                    std::swap(tNear[i], tFar[i]);
                dist0 = tNear[i] > dist0 ? tNear[i] : dist0;
                dist1 = tFar[i] < dist1 ? tFar[i] : dist1;
                if (dist0 > dist1)
                    return false;
            }
            if (contains(r.org)) {
                *distToBoundary = dist1;
                *enter = false;
                return dist1 >= r.distMin && dist1 < r.distMax;
            }
            else {
                *distToBoundary = dist0;
                *enter = true;
                return dist0 >= r.distMin && dist0 < r.distMax;
            }
        }
        
        bool hasNaN() const {
            return minP.hasNaN() || maxP.hasNaN();
        }
        
        bool hasInf() const {
            return minP.hasInf() || maxP.hasInf();
        }
    };
    
    inline BoundingBox3D calcUnion(const BoundingBox3D &b0, const BoundingBox3D &b1) {
        return BoundingBox3D(min(b0.minP, b1.minP), max(b0.maxP, b1.maxP));
    }
    
    inline BoundingBox3D intersection(const BoundingBox3D &b0, const BoundingBox3D &b1) {
        return BoundingBox3D(max(b0.minP, b1.minP), min(b0.maxP, b1.maxP));
    }
    
    
    
    struct SLR_API Vertex {
        Point3D position;
        Normal3D normal;
        Tangent3D tangent;
        TexCoord2D texCoord;
        
        Vertex() { }
        Vertex(const Point3D &pos, const Normal3D &norm, const Tangent3D &tang, const TexCoord2D &tc) : position(pos), normal(norm), tangent(tang), texCoord(tc) { }
    };
    
    
    
    template <typename T, uint32_t MaxDepth = 10>
    struct FixedStack {
        T data[MaxDepth];
        uint32_t curSize;
        
        FixedStack() : curSize(0) { }
        FixedStack(const FixedStack &fs) {
            for (int i = 0; i < fs.curSize; ++i)
                data[i] = fs.data[i];
            curSize = fs.curSize;
        }
        FixedStack& operator=(const FixedStack &fs) {
            for (int i = 0; i < fs.curSize; ++i)
                data[i] = fs.data[i];
            curSize = fs.curSize;
            return *this;
        }
        
        void push(const T &obj) { data[curSize++] = obj; }
        const T &top() const { return data[curSize - 1]; }
        void pop() { --curSize; }
        uint32_t size() { return curSize; }
    };
    
    template <typename T>
    struct ScopedPop {
        FixedStack<T> &stack;
        T item;

        ScopedPop(FixedStack<T> &_stack) : stack(_stack) {
            item = stack.top();
            stack.pop();
        }
        ~ScopedPop() {
            stack.push(item);
        }
    };
    
    class SLR_API Interaction {
        friend class SurfacePoint;
        friend class MediumPoint;
        
        float m_time;
        float m_dist;
        Point3D m_p;
    public:
        Interaction(float time, float dist, const Point3D &p) :
        m_time(time), m_dist(dist), m_p(p)
        {}
        
        float getTime() const { return m_time; }
        float getDistance() const { return m_dist; }
        
        virtual InteractionPoint* createInteractionPoint(ArenaAllocator &mem) const = 0;
        virtual Light* createLight(ArenaAllocator &mem) const = 0;
    };
    
    class SLR_API SurfaceInteraction : public Interaction {
        friend class SurfacePoint;
        
        mutable FixedStack<const SurfaceObject*> m_hierarchy;
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
        
        void push(const SurfaceObject* obj) { m_hierarchy.push(obj); }
        ScopedPop<const SurfaceObject*> scopedPop() const { return ScopedPop<const SurfaceObject*>(m_hierarchy); }
        const SurfaceObject* top() const { return m_hierarchy.top(); }
        FixedStack<const SurfaceObject*> getHierarchy() const { return m_hierarchy; }
        uint32_t getObjectDepth() const { return m_hierarchy.size(); }
        
        const Normal3D &getGeometricNormal() const { return m_gNormal; }
        void getSurfaceParameter(float* u, float* v) const {
            *u = m_u;
            *v = m_v;
        }
        
        void getSurfacePoint(SurfacePoint* surfPt) const;
        InteractionPoint* createInteractionPoint(ArenaAllocator &mem) const override;
        Light* createLight(ArenaAllocator &mem) const override;
    };
    
    class SLR_API MediumInteraction : public Interaction {
        friend class MediumPoint;
        
        mutable FixedStack<const MediumObject*> m_hierarchy;
        float m_u, m_v, m_t;
    public:
        MediumInteraction() : Interaction(0.0f, INFINITY, Point3D::Zero)
        {}
        MediumInteraction(float time, float dist, const Point3D &p) :
        Interaction(time, dist, p)
        {}
        
        ScopedPop<const MediumObject*> scopedPop() const { return ScopedPop<const MediumObject*>(m_hierarchy); }
        const MediumObject* top() const { return m_hierarchy.top(); }
        FixedStack<const MediumObject*> getHierarchy() const { return m_hierarchy; }
        uint32_t getObjectDepth() const { return m_hierarchy.size(); }
        
        void getMediumParameter(float* u, float* v, float* t) const {
            *u = m_u;
            *v = m_v;
            *t = m_t;
        }
        
        void getMediumPoint(MediumPoint* medPt) const;
        InteractionPoint* createInteractionPoint(ArenaAllocator &mem) const override;
        Light* createLight(ArenaAllocator &mem) const override;
    };
    
    
    
    struct SLR_API ReferenceFrame {
        Vector3D x, y, z;
        
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
        virtual SampledSpectrum fluxDensity(const WavelengthSamples &wls) const = 0;
        virtual float evaluateSpatialPDF() const = 0;
        virtual EDF* createEDF(const WavelengthSamples &wls, ArenaAllocator &mem) const = 0;
        
        virtual ABDFQuery* createABDFQuery(const Vector3D &dirLocal, int16_t selectedWL, DirectionType dirType, bool adjoint, ArenaAllocator &mem) const = 0;
        virtual AbstractBDF* createAbstractBDF(const WavelengthSamples &wls, ArenaAllocator &mem) const = 0;
        virtual SampledSpectrum evaluateInteractance(const WavelengthSamples &wls) const = 0;
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
        
        SampledSpectrum emittance(const WavelengthSamples &wls) const;
        float evaluateAreaPDF() const;
        BSDF* createBSDF(const WavelengthSamples &wls, ArenaAllocator &mem) const;
        
        //----------------------------------------------------------------
        // InteractionPoint's methods
        bool isEmitting() const override;
        SampledSpectrum fluxDensity(const WavelengthSamples &wls) const override {
            return emittance(wls);
        }
        float evaluateSpatialPDF() const override {
            return evaluateAreaPDF();
        }
        EDF* createEDF(const WavelengthSamples &wls, ArenaAllocator &mem) const override;
        
        ABDFQuery* createABDFQuery(const Vector3D &dirLocal, int16_t selectedWL, DirectionType filter, bool adjoint, ArenaAllocator &mem) const override;
        AbstractBDF* createAbstractBDF(const WavelengthSamples &wls, ArenaAllocator &mem) const override {
            return (AbstractBDF*)createBSDF(wls, mem);
        }
        SampledSpectrum evaluateInteractance(const WavelengthSamples &wls) const override {
            return SampledSpectrum::One;
        }
        float calcCosTerm(const Vector3D &vecWorld) const override {
            return absDot(vecWorld, m_gNormal);
        }
        
        void applyTransform(const StaticTransform &transform) override;
        //----------------------------------------------------------------
    };
    
    inline float squaredDistance(const SurfacePoint &p0, const SurfacePoint &p1) {
        return (p0.atInfinity() || p1.atInfinity()) ? 1.0f : sqDistance(p0.getPosition(), p1.getPosition());
    }
    
    class SLR_API MediumPoint : public InteractionPoint {
        const SingleMediumObject* m_obj;
    public:
        void setObject(const SingleMediumObject* obj) { m_obj = obj; }
        
        SampledSpectrum emittance(const WavelengthSamples &wls) const;
        float evaluateVolumePDF() const;
        BSDF* createPhaseFunction(const WavelengthSamples &wls, ArenaAllocator &mem) const;
        
        //----------------------------------------------------------------
        // InteractionPoint's methods
        bool isEmitting() const override;
        SampledSpectrum fluxDensity(const WavelengthSamples &wls) const override {
            return emittance(wls);
        }
        float evaluateSpatialPDF() const override {
            return evaluateVolumePDF();
        }
        EDF* createEDF(const WavelengthSamples &wls, ArenaAllocator &mem) const override;
        
        ABDFQuery* createABDFQuery(const Vector3D &dirLocal, int16_t selectedWL, DirectionType filter, bool adjoint, ArenaAllocator &mem) const override;
        AbstractBDF* createAbstractBDF(const WavelengthSamples &wls, ArenaAllocator &mem) const override;
        SampledSpectrum evaluateInteractance(const WavelengthSamples &wls) const override;
        float calcCosTerm(const Vector3D &vecWorld) const override {
            return 1.0f;
        }
        
        void applyTransform(const StaticTransform &transform) override;
        //----------------------------------------------------------------
    };
    
    
    
    // SurfaceShape?
    class SLR_API Surface {
    public:
        virtual ~Surface() { }
        
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
        virtual bool intersect(const Ray &ray, SurfaceInteraction* si) const = 0;
        virtual void getSurfacePoint(const SurfaceInteraction &si, SurfacePoint* surfPt) const = 0;
        virtual float area() const = 0;
        virtual void sample(float u0, float u1, SurfacePoint* surfPt, float* areaPDF) const = 0;
        virtual float evaluateAreaPDF(const SurfacePoint& surfPt) const = 0;
    };
    
    // MediumDistribution?
    class SLR_API Medium {
    protected:
        float m_maxExtinctionCoefficient;
    public:
        Medium(float maxExtinctionCoefficient) : m_maxExtinctionCoefficient(maxExtinctionCoefficient) { }
        virtual ~Medium() { }
        
        float majorantExtinctionCoefficient() const { return m_maxExtinctionCoefficient; };
        
        virtual bool subdivide(Allocator* mem, Medium** fragments, uint32_t* numFragments) const = 0;
        
        virtual BoundingBox3D bounds() const = 0;
        virtual bool contains(const Point3D &p) const = 0;
        virtual bool intersectBoundary(const Ray &ray, float* distToBoundary, bool* enter) const = 0;
        // The term "majorant" comes from the paper of Residual Ratio Tracking.
        virtual SampledSpectrum extinctionCoefficient(const Point3D &p, const WavelengthSamples &wls) const = 0;
        virtual bool interact(const Ray &ray, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                              MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const = 0;
        virtual void getMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const = 0;
        virtual void queryCoefficients(const Point3D &p, const WavelengthSamples &wls, SampledSpectrum* sigma_s, SampledSpectrum* sigma_e) const = 0;
        virtual float volume() const = 0;
        virtual void sample(float u0, float u1, float u2, MediumPoint* medPt, float* volumePDF) const = 0;
        virtual float evaluateVolumePDF(const MediumPoint& medPt) const = 0;
    };
}

#endif

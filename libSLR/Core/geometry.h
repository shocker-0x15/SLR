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
    
    
    
    template <typename T>
    class ScopedPop {
        std::vector<T> &m_stack;
        T m_item;
    public:
        ScopedPop(std::vector<T> &stack) : m_stack(stack) {
            m_item = m_stack.back();
            m_stack.pop_back();
        }
        ~ScopedPop() {
            m_stack.push_back(m_item);
        }
    };
    
    // It is not safe to directly use the point and normal because they need to be applied several transforms.
    struct SLR_API Intersection {
    private:
        mutable std::vector<const SurfaceObject*> m_hierarchy;
    public:
        float time;
        float dist;
        Point3D p;
        Normal3D gNormal;
        float u, v;
        TexCoord2D texCoord;
        
        Intersection() : dist(INFINITY) { }
        
        void push(const SurfaceObject* obj) { m_hierarchy.push_back(obj); }
        ScopedPop<const SurfaceObject*> scopedPop() const { return ScopedPop<const SurfaceObject*>(m_hierarchy); }
        const SurfaceObject* top() const { return m_hierarchy.back(); }
        const std::vector<const SurfaceObject*> getHierarchy() const { return m_hierarchy; }
        size_t getObjectDepth() const { return m_hierarchy.size(); }
        
        const SurfaceMaterial* getSurfaceMaterial() const;
        void getSurfacePoint(SurfacePoint* surfPt) const;
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
    
    
    
    struct SLR_API SurfacePoint {
        Point3D p;
        bool atInfinity;
        Normal3D gNormal;
        float u, v;
        TexCoord2D texCoord;
        Vector3D texCoord0Dir;
        ReferenceFrame shadingFrame;
        const SingleSurfaceObject* obj;
        
        float getSquaredDistance(const Point3D &shadingPoint) const { return atInfinity ? 1.0f : sqDistance(p, shadingPoint); }
        Vector3D getDirectionFrom(const Point3D &shadingPoint, float* dist2) const;
        bool isEmitting() const;
        SampledSpectrum emittance(const WavelengthSamples &wls) const;
        float evaluateAreaPDF() const;
        BSDF* createBSDF(const WavelengthSamples &wls, ArenaAllocator &mem) const;
        EDF* createEDF(const WavelengthSamples &wls, ArenaAllocator &mem) const;
        
        friend SurfacePoint operator*(const StaticTransform &transform, const SurfacePoint &surfPt);
    };
    
    inline float squaredDistance(const SurfacePoint &p0, const SurfacePoint &p1) {
        return (p0.atInfinity || p1.atInfinity) ? 1.0f : sqDistance(p0.p, p1.p);
    }
    
    struct SLR_API MediumPoint {
        Point3D p;
        const MediumObject* obj;
        
        float getSquaredDistance(const Point3D &shadingPoint) const { return sqDistance(p, shadingPoint); }
        Vector3D getDirectionFrom(const Point3D &shadingPoint, float* dist2) const;
        bool isEmitting() const;
        SampledSpectrum fluxDensity(const WavelengthSamples &wls) const;
        float evaluateVolumePDF() const;
        BSDF* createPhaseFunction(const WavelengthSamples &wls, ArenaAllocator &mem) const;
        EDF* createEDF(const WavelengthSamples &wls, ArenaAllocator &mem) const;
    };
    
    
    
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
        virtual bool intersect(const Ray &ray, Intersection* isect) const = 0;
        virtual void getSurfacePoint(const Intersection &isect, SurfacePoint* surfPt) const = 0;
        virtual float area() const = 0;
        virtual void sample(float u0, float u1, SurfacePoint* surfPt, float* areaPDF) const = 0;
        virtual float evaluateAreaPDF(const SurfacePoint& surfPt) const = 0;
    };
}

#endif

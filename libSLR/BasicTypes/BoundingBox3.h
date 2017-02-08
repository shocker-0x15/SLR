//
//  BoundingBox3.h
//
//  Created by 渡部 心 on 2017/01/23.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef __SLR_BoundingBox3__
#define __SLR_BoundingBox3__

#include "Point3.h"
#include "Ray.h"

namespace SLR {
    template <typename RealType>
    struct SLR_API BoundingBox3Template {
        enum Axis : uint8_t {
            Axis_X = 0,
            Axis_Y,
            Axis_Z,
        };
        
        Point3Template<RealType> minP, maxP;
        
        BoundingBox3Template() {
            minP.x = minP.y = minP.z = INFINITY;
            maxP.x = maxP.y = maxP.z = -INFINITY;
        }
        BoundingBox3Template(const Point3Template<RealType> &p) : minP(p), maxP(p) { }
        BoundingBox3Template(const Point3Template<RealType> &pmin, const Point3Template<RealType> &pmax) : minP(pmin), maxP(pmax) { }
        
        Point3Template<RealType> centroid() const {
            return (minP + maxP) * 0.5f;
        }
        
        RealType surfaceArea() const {
            Vector3Template<RealType> d = maxP - minP;
            return 2 * (d.x * d.y + d.y * d.z + d.z * d.x);
        }
        
        RealType volume() const {
            Vector3Template<RealType> d = maxP - minP;
            return d.x * d.y * d.z;
        }
        
        Point3Template<RealType> corner(uint32_t c) const {
            SLRAssert(c >= 0 && c < 8, "\"c\" is out of range [0, 8]");
            const size_t offset = sizeof(Point3Template<RealType>);
            return Point3Template<RealType>(*(&minP.x + offset * (c & 0x01)),
                           *(&minP.y + offset * (c & 0x02)),
                           *(&minP.z + offset * (c & 0x04)));
        }
        
        RealType centerOfAxis(Axis axis) const {
            return (minP[axis] + maxP[axis]) * 0.5f;
        }
        
        RealType width(Axis axis) const {
            return maxP[axis] - minP[axis];
        }
        
        Axis widestAxis() const {
            Vector3Template<RealType> d = maxP - minP;
            if (d.x > d.y && d.x > d.z)
                return Axis_X;
            else if (d.y > d.z)
                return Axis_Y;
            else
                return Axis_Z;
        }
        
        bool isValid() const {
            Vector3Template<RealType> d = maxP - minP;
            return d.x >= 0 && d.y >= 0 && d.z >= 0;
        }
        
        BoundingBox3Template &unify(const Point3Template<RealType> &p) {
            minP = min(minP, p);
            maxP = max(maxP, p);
            return *this;
        }
        
        BoundingBox3Template unify(const BoundingBox3Template &b) {
            minP = min(minP, b.minP);
            maxP = max(maxP, b.maxP);
            return *this;
        }
        
        bool contains(const Point3Template<RealType> &p) const {
            return ((p.x >= minP.x && p.x < maxP.x) &&
                    (p.y >= minP.y && p.y < maxP.y) &&
                    (p.z >= minP.z && p.z < maxP.z));
        }
        
        void localCoordinates(const Point3Template<RealType> &p, Point3Template<RealType>* param) const {
            *param = (p - minP) / (maxP - minP);
        }
        
        // check intersection between ray segment and bounding volume.
        bool intersect(const RayTemplate<RealType> &r) const {
            RealType dist0 = r.distMin, dist1 = r.distMax;
            Vector3Template<RealType> invRayDir = r.dir.reciprocal();
            Vector3Template<RealType> tNear = (minP - r.org) * invRayDir;
            Vector3Template<RealType> tFar = (maxP - r.org) * invRayDir;
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
        bool intersectBoundary(const RayTemplate<RealType> &r, RealType* distToBoundary, bool* enter) const {
            RealType dist0 = -INFINITY, dist1 = INFINITY;
            Vector3Template<RealType> invRayDir = r.dir.reciprocal();
            Vector3Template<RealType> tNear = (minP - r.org) * invRayDir;
            Vector3Template<RealType> tFar = (maxP - r.org) * invRayDir;
            for (int i = 0; i < 3; ++i) {
                if (tNear[i] > tFar[i])
                    std::swap(tNear[i], tFar[i]);
                dist0 = tNear[i] > dist0 ? tNear[i] : dist0;
                dist1 = tFar[i] < dist1 ? tFar[i] : dist1;
                if (dist0 > dist1)
                    return false;
            }
            if (contains(r.org + r.distMin * r.dir)) {
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
    
    template <typename RealType>
    inline BoundingBox3Template<RealType> calcUnion(const BoundingBox3Template<RealType> &b0, const BoundingBox3Template<RealType> &b1) {
        return BoundingBox3Template<RealType>(min(b0.minP, b1.minP), max(b0.maxP, b1.maxP));
    }
    
    template <typename RealType>
    inline BoundingBox3Template<RealType> intersection(const BoundingBox3Template<RealType> &b0, const BoundingBox3Template<RealType> &b1) {
        return BoundingBox3Template<RealType>(max(b0.minP, b1.minP), min(b0.maxP, b1.maxP));
    }
}

#endif /* __SLR_BoundingBox3__ */

//
//  Transform.h
//
//  Created by 渡部 心 on 2015/05/04.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__Transform__
#define __SLR__Transform__

#include "../defines.h"
#include "../references.h"
#include "geometry.h"

namespace SLR {
    class Transform {
    public:
        virtual ~Transform() { };
        
        virtual bool isStatic() const = 0;
        virtual void sample(float time, StaticTransform* tf) const { SLRAssert(false, "Not implemented."); };
        virtual bool isChained() const = 0;
        virtual Transform* copy(ArenaAllocator &mem) const = 0;
        virtual BoundingBox3D motionBounds(const BoundingBox3D &bb) const { SLRAssert(false, "Not implemented."); return BoundingBox3D(); };
        
        Ray mul(const Ray &r) const;
        Point3D mul(const Point3D &p, float t) const;
        Normal3D mul(const Normal3D &n, float t) const;
        
        Ray mulInv(const Ray &r) const;
        Point3D mulInv(const Point3D &p, float t) const;
        Normal3D mulInv(const Normal3D &n, float t) const;
    };
    
    class StaticTransform : public Transform {
        Matrix4x4 mat, matInv;
    public:
        StaticTransform(const Matrix4x4 &m = Matrix4x4::Identity) : mat(m), matInv(invert(m)) { };
        StaticTransform(const Matrix4x4 &m, const Matrix4x4 &mInv) : mat(m), matInv(mInv) { };
        
        Vector3D operator*(const Vector3D &v) const { return mat * v; };
        Vector4D operator*(const Vector4D &v) const { return mat * v; };
        Point3D operator*(const Point3D &p) const { return mat * p; };
        Normal3D operator*(const Normal3D &n) const {
            // The length of the normal is changed if the transform has scaling, so it requires normalization.
            return Normal3D(matInv.m00 * n.x + matInv.m10 * n.y + matInv.m20 * n.z,
                            matInv.m01 * n.x + matInv.m11 * n.y + matInv.m21 * n.z,
                            matInv.m02 * n.x + matInv.m12 * n.y + matInv.m22 * n.z);
        };
        Ray operator*(const Ray &r) const { return Ray(mat * r.org, mat * r.dir, r.time, r.distMin, r.distMax); };
        BoundingBox3D operator*(const BoundingBox3D &bb) const {
            BoundingBox3D ret;
            ret.unify(mat * Point3D(bb.minP.x, bb.minP.y, bb.minP.z));
            ret.unify(mat * Point3D(bb.minP.x, bb.minP.y, bb.maxP.z));
            ret.unify(mat * Point3D(bb.minP.x, bb.maxP.y, bb.minP.z));
            ret.unify(mat * Point3D(bb.minP.x, bb.maxP.y, bb.maxP.z));
            ret.unify(mat * Point3D(bb.maxP.x, bb.minP.y, bb.minP.z));
            ret.unify(mat * Point3D(bb.maxP.x, bb.minP.y, bb.maxP.z));
            ret.unify(mat * Point3D(bb.maxP.x, bb.maxP.y, bb.minP.z));
            ret.unify(mat * Point3D(bb.maxP.x, bb.maxP.y, bb.maxP.z));
            return ret;
        };
        StaticTransform operator*(const Matrix4x4 &m) const { return StaticTransform(mat * m); };
        StaticTransform operator*(const StaticTransform &t) const { return StaticTransform(mat * t.mat); };
        bool operator==(const StaticTransform &t) const { return mat == t.mat; };
        bool operator!=(const StaticTransform &t) const { return mat != t.mat; };
        
        Matrix4x4 getMatrix4x4() const { return mat; };
        bool isIdentity() const { return mat.isIdentity(); };
        
        
        bool isStatic() const override { return true; };
        void sample(float time, StaticTransform* tf) const override { *tf = *this; };
        bool isChained() const override { return false; };
        Transform* copy(ArenaAllocator &mem) const override;
        BoundingBox3D motionBounds(const BoundingBox3D &bb) const override { return *this * bb; };
        
        
        friend StaticTransform invert(const StaticTransform &t) { return StaticTransform(t.matInv, t.mat); };
        friend StaticTransform transpose(const StaticTransform &t) { return StaticTransform(transpose(t.mat)); };
    };
    
    class AnimatedTransform : public Transform {
        StaticTransform* m_tfBegin;
        StaticTransform* m_tfEnd;
        float m_tBegin, m_tEnd;
        Vector3D m_T[2];
        Quaternion m_R[2];
        Matrix4x4 m_S[2];
    public:
        AnimatedTransform(StaticTransform* tfBegin, StaticTransform* tfEnd, float tBegin, float tEnd) :
        m_tfBegin(tfBegin), m_tfEnd(tfEnd), m_tBegin(tBegin), m_tEnd(tEnd) {
            decompose(tfBegin->getMatrix4x4(), &m_T[0], &m_R[0], &m_S[0]);
            decompose(tfEnd->getMatrix4x4(), &m_T[1], &m_R[1], &m_S[1]);
        };
        
        bool isStatic() const override { return *m_tfBegin == *m_tfEnd; };
        
        void sample(float time, StaticTransform* tf) const override {
            if (time <= m_tBegin) {
                *tf = *m_tfBegin;
                return;
            }
            if (time >= m_tEnd) {
                *tf = *m_tfEnd;
                return;
            }
            float t = (time - m_tBegin) / (m_tEnd - m_tBegin);
            Vector3D trans = (1 - t) * m_T[0] + t * m_T[1];
            
            Quaternion rotate = Slerp(t, m_R[0], m_R[1]);
            
            Matrix4x4 scale = (1 - t) * m_S[0] + t * m_S[1];
            
            *tf = translate(trans) * rotate.toMatrix() * scale;
        };
        
        bool isChained() const override { return false; };
        
        Transform* copy(ArenaAllocator &mem) const override;
        
        // FIXME: This sampling-based way does not guarantee to generate the actual bounds.
        BoundingBox3D motionBounds(const BoundingBox3D &bb) const override {
            BoundingBox3D ret;
            const uint32_t numIte = 128;
            for (uint32_t i = 0; i < numIte; ++i) {
                float t = (float)i / (numIte - 1);
                float sTime = (1 - t) * m_tBegin + t * m_tEnd;
                StaticTransform staticTF;
                sample(sTime, &staticTF);
                ret.unify(staticTF * bb);
            }
            return ret;
        };
    };
    
    class ChainedTransform : public SLR::Transform {
        const SLR::Transform* m_parent;
        const SLR::Transform* m_transform;
    public:
        ChainedTransform(const SLR::Transform* parent, const SLR::Transform* transform) : m_parent(parent), m_transform(transform) {
            SLRAssert(transform->isChained() == false, "\"transform\" cannot be a ChainedTransform.");
        };
        
        void setParent(const SLR::Transform* parent) { m_parent = parent; };
        void setTransform(const SLR::Transform* transform) { m_transform = transform; };
        
        void reduce(SLR::ArenaAllocator &mem, std::vector<const SLR::Transform*> &reducedTFs) const;
        const SLR::Transform* reduce(SLR::ArenaAllocator &mem) const;
        
        bool isStatic() const override {
            bool ret = true;
            const SLR::Transform* current = this;
            while (ret && current) {
                const SLR::Transform* parent = nullptr;
                if (current->isChained()) {
                    parent = ((ChainedTransform*)current)->m_parent;
                    current = ((ChainedTransform*)current)->m_transform;
                }
                ret &= current->isStatic();
                
                current = parent;
            }
            return ret;
        };
        
        bool isChained() const override { return true; };

        SLR::Transform* copy(SLR::ArenaAllocator &mem) const override;
    };
}

#endif

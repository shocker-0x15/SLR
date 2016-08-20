//
//  Transform.cpp
//
//  Created by 渡部 心 on 2015/05/04.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "Transform.h"
#include "../Memory/ArenaAllocator.h"

namespace SLR {
    Ray Transform::mul(const Ray &r) const {
        StaticTransform tf;
        sample(r.time, &tf);
        return tf * r;
    }
    
    Point3D Transform::mul(const Point3D &p, float t) const {
        StaticTransform tf;
        sample(t, &tf);
        return tf * p;
    }
    
    Normal3D Transform::mul(const Normal3D &n, float t) const {
        StaticTransform tf;
        sample(t, &tf);
        return tf * n;
    }
    
    Ray Transform::mulInv(const Ray &r) const {
        StaticTransform tf;
        sample(r.time, &tf);
        return invert(tf) * r;
    }
    
    Point3D Transform::mulInv(const Point3D &p, float t) const {
        StaticTransform tf;
        sample(t, &tf);
        return invert(tf) * p;
    }
    
    Normal3D Transform::mulInv(const Normal3D &n, float t) const {
        StaticTransform tf;
        sample(t, &tf);
        return invert(tf) * n;
    }
    
    Transform* StaticTransform::copy() const {
        return new StaticTransform(*this);
    }
    Transform* StaticTransform::copy(ArenaAllocator &mem) const {
        return mem.create<StaticTransform>(*this);
    }
    Transform* StaticTransform::createByMulLeft(const StaticTransform &staticTF, ArenaAllocator &mem) const {
        return mem.create<StaticTransform>(staticTF * *this);
    }
    Transform* StaticTransform::createByMulRight(const StaticTransform &staticTF, ArenaAllocator &mem) const {
        return mem.create<StaticTransform>(*this * staticTF);
    }
    
    Transform* AnimatedTransform::copy() const {
        return new AnimatedTransform(*this);
    }
    Transform* AnimatedTransform::copy(ArenaAllocator &mem) const {
        return mem.create<AnimatedTransform>(*this);
    }
    Transform* AnimatedTransform::createByMulLeft(const StaticTransform &staticTF, ArenaAllocator &mem) const {
        return mem.create<AnimatedTransform>(staticTF * m_tfBegin, staticTF * m_tfEnd, m_tBegin, m_tEnd);
    }
    Transform* AnimatedTransform::createByMulRight(const StaticTransform &staticTF, ArenaAllocator &mem) const {
        return mem.create<AnimatedTransform>(m_tfBegin * staticTF, m_tfEnd * staticTF, m_tBegin, m_tEnd);
    }
    
    void ChainedTransform::reduce(SLR::ArenaAllocator &mem, std::vector<const SLR::Transform*> &reducedTFs) const {
        using namespace SLR;
        const Transform* linearTFs[10];
        uint32_t idx = 0;
        const Transform* current = this;
        while (current) {
            const Transform* parent = nullptr;
            if (current->isChained()) {
                parent = ((ChainedTransform*)current)->m_parent;
                current = ((ChainedTransform*)current)->m_transform;
            }
            linearTFs[idx++] = current;
            
            current = parent;
        }
        
        uint32_t numLinearTFs = idx;
        for (int i = 0; i < numLinearTFs; ++i) {
            current = linearTFs[i];
            if (reducedTFs.size() == 0) {
                reducedTFs.push_back(current->copy(mem));
                continue;
            }
            
            const SLR::Transform* top = reducedTFs.back();
            bool mergeable = current->isStatic() || top->isStatic();
            if (mergeable) {
                if (current->isStatic() && top->isStatic()) {
                    StaticTransform curSample;
                    current->sample(0.0f, &curSample);
                    StaticTransform topSample;
                    top->sample(0.0f, &topSample);
                    reducedTFs.pop_back();
                    reducedTFs.push_back((curSample * topSample).copy(mem));
                }
                else if (current->isStatic()) {
                    StaticTransform curSample;
                    current->sample(0.0f, &curSample);
                    Transform* animated = top->createByMulLeft(curSample, mem);
                    reducedTFs.pop_back();
                    reducedTFs.push_back(animated);
                }
                else {
                    StaticTransform topSample;
                    top->sample(0.0f, &topSample);
                    Transform* animated = current->createByMulRight(topSample, mem);
                    reducedTFs.pop_back();
                    reducedTFs.push_back(animated);
                }
            }
            else {
                reducedTFs.push_back(current->copy(mem));
            }
        }
    }
    
    const SLR::Transform* ChainedTransform::reduce(SLR::ArenaAllocator &mem) const {
        using namespace SLR;
        std::vector<const Transform*> reducedTFs;
        reduce(mem, reducedTFs);
        
        const Transform* ret = nullptr;
        if (reducedTFs.size() > 1) {
            ChainedTransform* chain = mem.create<ChainedTransform>(nullptr, reducedTFs[0]);
            ChainedTransform* prevChain = chain;
            ret = chain;
            for (int i = 1; i < reducedTFs.size() - 1; ++i) {
                chain = mem.create<ChainedTransform>(nullptr, reducedTFs[i]);
                prevChain->setParent(chain);
                prevChain = chain;
            }
            prevChain->setParent(reducedTFs.back());
        }
        else if (reducedTFs.size() == 1) {
            ret = reducedTFs[0];
        }
        
        return ret;
    }
    
    SLR::Transform* ChainedTransform::copy() const {
        return new ChainedTransform(*this);
    }
    SLR::Transform* ChainedTransform::copy(SLR::ArenaAllocator &mem) const {
        return mem.create<ChainedTransform>(*this);
    }
}

//
//  Transform.cpp
//
//  Created by 渡部 心 on 2015/05/04.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "Transform.h"
#include "../Memory/ArenaAllocator.h"

Ray Transform::mul(const Ray &r) const {
    StaticTransform tf;
    sample(r.time, &tf);
    return tf * r;
};

Point3D Transform::mul(const Point3D &p, float t) const {
    StaticTransform tf;
    sample(t, &tf);
    return tf * p;
};

Normal3D Transform::mul(const Normal3D &n, float t) const {
    StaticTransform tf;
    sample(t, &tf);
    return tf * n;
};

Ray Transform::mulInv(const Ray &r) const {
    StaticTransform tf;
    sample(r.time, &tf);
    return invert(tf) * r;
};

Point3D Transform::mulInv(const Point3D &p, float t) const {
    StaticTransform tf;
    sample(t, &tf);
    return invert(tf) * p;
};

Normal3D Transform::mulInv(const Normal3D &n, float t) const {
    StaticTransform tf;
    sample(t, &tf);
    return invert(tf) * n;
};

Transform* StaticTransform::copy(ArenaAllocator &mem) const {
    return mem.create<StaticTransform>(*this);
};

Transform* AnimatedTransform::copy(ArenaAllocator &mem) const {
    return mem.create<AnimatedTransform>(*this);
};

Transform* ChainedTransform::copy(ArenaAllocator &mem) const {
    return mem.create<ChainedTransform>(*this);
};

void ChainedTransform::reduce(ArenaAllocator &mem, std::vector<const Transform*> &reducedTFs) const {
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
    current = linearTFs[idx];
    uint32_t numChaindS = 0;
    StaticTransform cumS;
    for (int i = 0; i < numLinearTFs; ++i) {
        current = linearTFs[i];
        if (current->isStatic()) {
            StaticTransform staticTF;
            current->sample(0.0f, &staticTF);
            cumS = staticTF * cumS;
            if (cumS.isIdentity()) {
                numChaindS = 0;
            }
            else {
                ++numChaindS;
                if (i + 1 >= numLinearTFs || linearTFs[i + 1]->isStatic() == false) {
                    if (numChaindS == 1)
                        reducedTFs.push_back(current);
                    else
                        reducedTFs.push_back(cumS.copy(mem));
                    cumS = StaticTransform();
                    numChaindS = 0;
                }
            }
        }
        else {
            reducedTFs.push_back(linearTFs[i]->copy(mem));
        }
    }
}

const Transform* ChainedTransform::reduce(ArenaAllocator &mem) const {
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

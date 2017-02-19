//
//  transform.cpp
//
//  Created by 渡部 心 on 2015/05/04.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "transform.h"

#include "../MemoryAllocators/Allocator.h"

namespace SLR {
    Transform* StaticTransform::copy(Allocator* mem) const {
        return mem->create<StaticTransform>(*this);
    }
    
    
    
    Transform* AnimatedTransform::copy(Allocator* mem) const {
        return mem->create<AnimatedTransform>(*this);
    }
    
    
    
    void ChainedTransform::reduce(Allocator* mem, std::vector<Transform*> &reducedTFs) const {
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
        
        current = linearTFs[0];
        reducedTFs.push_back(current->copy(mem));
        
        uint32_t numLinearTFs = idx;
        for (int i = 1; i < numLinearTFs; ++i) {
            current = linearTFs[i];
            const Transform* child = reducedTFs.back();
            bool mergeable = current->isStatic() || child->isStatic();
            if (mergeable) {
                mem->destroy(reducedTFs.back());
                reducedTFs.pop_back();
                if (current->isStatic() && child->isStatic()) {
                    StaticTransform curSample;
                    current->sample(0.0f, &curSample);
                    StaticTransform childSample;
                    child->sample(0.0f, &childSample);
                    reducedTFs.push_back((curSample * childSample).copy(mem));
                }
                else if (current->isStatic()) {
                    StaticTransform curSample;
                    current->sample(0.0f, &curSample);
                    AnimatedTransform &atf = *(AnimatedTransform*)child;
                    reducedTFs.push_back((curSample * atf).copy(mem));
                }
                else {
                    StaticTransform childSample;
                    child->sample(0.0f, &childSample);
                    AnimatedTransform &atf = *(AnimatedTransform*)current;
                    reducedTFs.push_back((atf * childSample).copy(mem));
                }
            }
            else {
                reducedTFs.push_back(current->copy(mem));
            }
        }
    }
    
    Transform* ChainedTransform::reduce(Allocator* mem) const {
        using namespace SLR;
        std::vector<Transform*> reducedTFs;
        reduce(mem, reducedTFs);
        
        Transform* ret = nullptr;
        if (reducedTFs.size() > 1) {
            ChainedTransform* chain = mem->create<ChainedTransform>(nullptr, reducedTFs[0]);
            ChainedTransform* prevChain = chain;
            ret = chain;
            for (int i = 1; i < reducedTFs.size() - 1; ++i) {
                chain = mem->create<ChainedTransform>(nullptr, reducedTFs[i]);
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
    
    Transform* ChainedTransform::copy(Allocator* mem) const {
        return mem->create<ChainedTransform>(this->m_parent, this->m_transform);
    }
}

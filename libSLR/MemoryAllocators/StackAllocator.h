//
//  StackAllocator.h
//
//  Created by 渡部 心 on 2015/04/28.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__StackAllocator__
#define __SLR__StackAllocator__

#include "Allocator.h"
#include <unordered_map>

class StackAllocator : public Allocator {
    bool m_initialized;
    void* m_poolHead;
    uintptr_t m_poolSize;
#define maxNumAllocs 8192
    uint8_t* m_nexts[maxNumAllocs];
    std::unordered_map<void*, uint32_t> m_idxMap;
    bool m_freeable[maxNumAllocs + 1];
    uint32_t m_idx;
public:
    StackAllocator();
    StackAllocator(uintptr_t poolSize);
    ~StackAllocator();
    
    void init(uintptr_t poolSize);
    
    void* alloc(uintptr_t size, uintptr_t align) override;
    void free(void* ptr) override;
};

#endif

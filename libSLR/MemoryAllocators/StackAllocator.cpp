//
//  StackAllocator.cpp
//
//  Created by 渡部 心 on 2015/04/28.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "StackAllocator.h"

StackAllocator::StackAllocator() : m_initialized(false), m_poolHead(nullptr), m_idx(0), m_idxMap(maxNumAllocs) { }

StackAllocator::StackAllocator(uintptr_t poolSize) : StackAllocator() {
    init(poolSize);
}

StackAllocator::~StackAllocator() {
    if (m_poolHead)
        delete[] (uint8_t*)m_poolHead;
}

void StackAllocator::init(uintptr_t poolSize) {
    if (m_poolHead)
        free(m_poolHead);
    m_poolSize = poolSize;
    m_poolHead = new uint8_t[m_poolSize];
    m_nexts[m_idx] = (uint8_t*)m_poolHead;
    m_freeable[0] = false;
    m_idxMap.clear();
    m_initialized = true;
}

void* StackAllocator::alloc(uintptr_t size, uintptr_t align) {
    SLRAssert(m_initialized, "Initialization is undone.");
    SLRAssert(size > 0, "Size \"size\" is zero.");
    SLRAssert((align & (align - 1)) == 0, "Alignment \"align\" must be a power of 2.");
    uintptr_t mask = align - 1;
    uint8_t* head = (uint8_t*)(((uintptr_t)m_nexts[m_idx] + mask) & ~mask);
    m_idxMap[head] = m_idx;
    m_freeable[m_idx] = false;
    m_nexts[++m_idx] = head + size;
    m_freeable[m_idx] = true;
    SLRAssert((uintptr_t)m_nexts[m_idx] - (uintptr_t)m_poolHead <= m_poolSize, "Allocation exceeds the pre-allocated pool size.");
    return head;
}

void StackAllocator::free(void *ptr) {
    SLRAssert(m_idxMap.count(ptr) > 0, "Pointer \"ptr\" does not exist or it has already released.");
    uint32_t idx = m_idxMap[ptr];
    m_idxMap.erase(ptr);
    if (idx + 1 == m_idx) {
        while (m_freeable[m_idx])
            --m_idx;
    }
    else {
        m_freeable[idx + 1] = true;
    }
}

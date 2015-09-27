//
//  MSpaceAllocator.cpp
//
//  Created by 渡部 心 on 2015/04/29.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "MSpaceAllocator.h"

MSpaceAllocator::MSpaceAllocator() : m_initialized(false), m_size(0), m_mspace(nullptr) {
    
}

MSpaceAllocator::MSpaceAllocator(uintptr_t poolSize) : MSpaceAllocator() {
    init(poolSize);
}

MSpaceAllocator::~MSpaceAllocator() {
    if (m_mspace)
        destroy_mspace(m_mspace);
}

void MSpaceAllocator::init(uintptr_t poolSize) {
    SLRAssert(poolSize > 0, "The pool size \"poolSize\" is zero.");
    if (m_mspace)
        destroy_mspace(m_mspace);
    m_size = poolSize;
    m_mspace = create_mspace(m_size, 0);
    m_initialized = true;
}

void* MSpaceAllocator::alloc(uintptr_t size, uint32_t align) {
    SLRAssert(m_initialized, "Initialization is undone.");
    SLRAssert(size > 0, "Size \"size\" is zero.");
    SLRAssert((align & (align - 1)) == 0, "Alignment \"align\" must be a power of 2.");
    return mspace_memalign(m_mspace, align, size);
}

void MSpaceAllocator::free(void *ptr) {
    SLRAssert(ptr, "\"ptr\" is nullptr.");
    mspace_free(m_mspace, ptr);
}

void MSpaceAllocator::printStats() const {
    mallinfo info = mspace_mallinfo(m_mspace);
    
    std::cout << "  - non-mmapped space allocated from system : " << info.arena << std::endl;
    std::cout << "  - number of free chunks                   : " << info.ordblks << std::endl;
    std::cout << "  - space in mmapped regions                : " << info.hblkhd << std::endl;
    std::cout << "  - maximum total allocated space           : " << info.usmblks << std::endl;
    std::cout << "  - total allocated space                   : " << info.uordblks << std::endl;
    std::cout << "  - total free space                        : " << info.fordblks << std::endl;
    std::cout << "  - releasable (via malloc_trim) space      : " << info.keepcost << std::endl;
}

uintptr_t MSpaceAllocator::offset(void* ptr) const {
    return (uintptr_t)ptr - (uintptr_t)m_mspace;
}

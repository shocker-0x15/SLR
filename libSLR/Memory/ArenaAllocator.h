//
//  ArenaAllocator.h
//
//  Created by 渡部 心 on 2015/07/08.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__ArenaAllocator__
#define __SLR__ArenaAllocator__

#include "../defines.h"
#include "../references.h"
#include "Allocator.h"
#include <list>

namespace SLR {
    // Highly based on MemoryArena used in PBRT-v3
    class SLR_API ArenaAllocator : public Allocator {
        const size_t m_blockSize;
        
        uint8_t *m_currentBlock = nullptr;
        size_t m_currentAllocSize = 0;
        size_t m_currentBlockPos = 0;
        
        std::list<std::pair<size_t, uint8_t*>> m_usedBlocks, m_availableBlocks;
        
        ArenaAllocator(const ArenaAllocator &) = delete;
        ArenaAllocator &operator=(const ArenaAllocator &) = delete;
        
        template <typename T>
        T* allocateNewBlock(size_t count) const {
            return (T*)SLR_memalign(count * sizeof(T), SLR_L1_Cacheline_Size);
        }
        template <typename T>
        T* allocateNewBlock(size_t count, size_t align) const {
            return (T*)SLR_memalign(count * sizeof(T), std::max<size_t>(align, SLR_L1_Cacheline_Size));
        }
        void deallocateBlock(void* ptr) const {
            if (!ptr)
                return;
            SLR_freealign(ptr);
        }
    public:
        ArenaAllocator(size_t blockSize = 262144) : m_blockSize(blockSize) { }
        ~ArenaAllocator();
        
        void* alloc(uintptr_t size, uintptr_t align) override;
        void free(void* ptr) override {
            // ArenaAllocator assumes a user should not use free().
        }

        void reset();
        size_t totalAllocated() const;
    };    
}

#endif

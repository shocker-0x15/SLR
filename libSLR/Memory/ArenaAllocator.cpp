//
//  ArenaAllocator.cpp
//
//  Created by 渡部 心 on 2015/07/08.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "ArenaAllocator.h"

namespace SLR {
    ArenaAllocator::~ArenaAllocator() {
        deallocateBlock(m_currentBlock);
        for (auto &block : m_usedBlocks)
            deallocateBlock(block.second);
        for (auto &block : m_availableBlocks)
            deallocateBlock(block.second);
    };
    
    void* ArenaAllocator::alloc(uintptr_t size, uintptr_t align) {
        SLRAssert(size > 0, "Size \"size\" is zero.");
        SLRAssert((align & (align - 1)) == 0, "Alignment \"align\" must be a power of 2.");
        uintptr_t mask = align - 1;
        uint8_t* curPos = m_currentBlock + m_currentBlockPos;
        uint8_t* head = (uint8_t*)(((uintptr_t)curPos + mask) & ~mask);
        size_t offset = head - curPos;
        
        if (offset + size > m_currentAllocSize) {
            if (m_currentBlock) {
                m_usedBlocks.push_back(std::make_pair(m_currentAllocSize, m_currentBlock));
                m_currentBlock = nullptr;
                m_currentAllocSize = 0;
            }
            
            for (auto iter = m_availableBlocks.begin(); iter != m_availableBlocks.end(); ++iter) {
                if (((uintptr_t)iter->second & mask) == 0 && iter->first >= size) {
                    m_currentAllocSize = iter->first;
                    m_currentBlock = iter->second;
                    m_availableBlocks.erase(iter);
                    break;
                }
            }
            if (!m_currentBlock) {
                m_currentAllocSize = std::max(size, m_blockSize);
                m_currentBlock = allocateNewBlock<uint8_t>(m_currentAllocSize, align);
            }
            m_currentBlockPos = 0;
        }
        m_currentBlockPos += size;
        return head;
    };
    
    void ArenaAllocator::free(void* ptr) {
        // ArenaAllocator assumes a user should not use free().
    };
    
    void* ArenaAllocator::alloc(size_t numBytes) {
        numBytes = ((numBytes + SLR_Minimum_Machine_Alignment) & (~(SLR_Minimum_Machine_Alignment)));
        if (m_currentBlockPos + numBytes > m_currentAllocSize) {
            if (m_currentBlock) {
                m_usedBlocks.push_back(std::make_pair(m_currentAllocSize, m_currentBlock));
                m_currentBlock = nullptr;
                m_currentAllocSize = 0;
            }
            
            for (auto iter = m_availableBlocks.begin(); iter != m_availableBlocks.end(); ++iter) {
                if (iter->first >= numBytes) {
                    m_currentAllocSize = iter->first;
                    m_currentBlock = iter->second;
                    m_availableBlocks.erase(iter);
                    break;
                }
            }
            if (!m_currentBlock) {
                m_currentAllocSize = std::max(numBytes, m_blockSize);
                m_currentBlock = allocateNewBlock<uint8_t>(m_currentAllocSize);
            }
            m_currentBlockPos = 0;
        }
        void* ret = m_currentBlock + m_currentBlockPos;
        m_currentBlockPos += numBytes;
        return ret;
    };
    
    void ArenaAllocator::reset() {
        m_currentBlockPos = 0;
        m_availableBlocks.splice(m_availableBlocks.begin(), m_usedBlocks);
    };
    
    size_t ArenaAllocator::totalAllocated() const {
        size_t total = m_currentAllocSize;
        for (const auto &alloc : m_usedBlocks)
            total += alloc.first;
        for (const auto &alloc : m_availableBlocks)
            total += alloc.first;
        return total;
    };    
}

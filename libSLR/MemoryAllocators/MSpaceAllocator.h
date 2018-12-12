//
//  MSpaceAllocator.h
//
//  Created by 渡部 心 on 2015/04/29.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__MSpaceAllocator__
#define __SLR__MSpaceAllocator__

#include "../defines.h"
#include "../references.h"
#include "dlmalloc.h"

class MSpaceAllocator {
    bool m_initialized;
    uint64_t m_size;
    mspace m_mspace;
    
    MSpaceAllocator(const MSpaceAllocator &mem) { };
public:
    MSpaceAllocator();
    MSpaceAllocator(uintptr_t poolSize);
    ~MSpaceAllocator();
    
    void init(uintptr_t poolSize);
    
    void* alloc(uintptr_t size, uint32_t align);
    template <typename T>
    T* alloc(uint64_t numElems = 1, uint32_t align = alignof(T)) {
        return (T*)alloc(numElems * sizeof(T), align);
    };
    template <typename T, typename ...ArgTypes>
    T* create(ArgTypes... args) {
        T* rawPtr = alloc<T>();
        new (rawPtr) T(args...);
        return rawPtr;
    };
    void free(void* ptr);
    
    void printStats() const;
    
    uintptr_t offset(void* ptr) const;
    
    template <typename T, typename ...ArgTypes>
    std::shared_ptr<T> createShared(ArgTypes&&... args) {
        T* rawPtr = alloc<T>();
        new (rawPtr) T(std::forward<ArgTypes>(args)...);
        std::function<void(void*)> deleter = std::bind(&MSpaceAllocator::free, this, std::placeholders::_1);
        return std::shared_ptr<T>(rawPtr, deleter);
    };
};


template <class T> class MSpaceSTLAllocator;

template <>
class MSpaceSTLAllocator<void> {
public:
    typedef void*             pointer;
    typedef const void*       const_pointer;
    typedef void              value_type;
    
    template <class U>
    struct rebind {
        typedef MSpaceSTLAllocator<U> other;
    };
};

template <class T>
class MSpaceSTLAllocator {
    MSpaceAllocator &m_allocator;
public:
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;
    typedef T value_type;
    
    template <class U>
    struct rebind {
        typedef MSpaceSTLAllocator<U> other;
    };
    
    template <class U>
    friend class MSpaceSTLAllocator;
    
    MSpaceSTLAllocator(MSpaceAllocator &allocator) throw() : m_allocator(allocator) {}
    ~MSpaceSTLAllocator() throw() {}
//    MSpaceSTLAllocator(const allocater&) throw() {}
    template <class U>
    MSpaceSTLAllocator(const MSpaceSTLAllocator<U>& stlAllocator) throw() : m_allocator(stlAllocator.m_allocator) {}
    
    pointer allocate(size_type num, MSpaceSTLAllocator<void>::const_pointer hint = 0) {
        return m_allocator.alloc<T>((uint32_t)num);
    }
    void construct(pointer p, const T& value) {
        new (p) T(value);
    }
    
    void deallocate(pointer p, size_type num) {
        m_allocator.free(p);
    }
    // 初期化済みの領域を削除する
    void destroy(pointer p) {
        p->~T();
    }
    
    // アドレスを返す
    pointer address(reference value) const { return &value; }
    const_pointer address(const_reference value) const { return &value; }
    
    // 割当てることができる最大の要素数を返す
    size_type max_size() const throw() {
        return std::numeric_limits<size_t>::max() / sizeof(T);
    }
};

template <class T> using MSpaceVector = std::vector<T, MSpaceSTLAllocator<T>>;
template <class T> using MSpaceSet = std::set<T, std::less<T>, MSpaceSTLAllocator<T>>;

#endif

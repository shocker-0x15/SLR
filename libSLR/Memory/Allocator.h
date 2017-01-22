//
//  Allocator.h
//
//  Created by 渡部 心 on 2015/06/29.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef SLR_Allocator_h
#define SLR_Allocator_h

#include "../defines.h"
#include "../references.h"

namespace SLR {
    class SLR_API Allocator {
    public:
        template <typename T> using DeleterType = std::function<void(T*)>;
        
        virtual void* alloc(uintptr_t size, uintptr_t align) = 0;
        virtual void free(void* ptr) = 0;
        
        template <typename T>
        T* alloc(size_t numElems = 1, uintptr_t align = alignof(T)) {
            return (T*)alloc(numElems * sizeof(T), align);
        }
        template <typename T, typename ...ArgTypes>
        T* create(ArgTypes&&... args) {
            T* ptr = alloc<T>();
            new (ptr) T(std::forward<ArgTypes>(args)...);
            return ptr;
        }
        template <typename T>
        void destroy(T* ptr) {
            ptr->~T();
            free(ptr);
        }
        
        template <typename T, typename ...ArgTypes>
        std::shared_ptr<T> createShared(ArgTypes&&... args) {
            T* rawPtr = alloc<T>();
            new (rawPtr) T(std::forward<ArgTypes>(args)...);
            std::function<void(void*)> deleter = std::bind(&Allocator::destroy<T>, this, std::placeholders::_1);
            return std::shared_ptr<T>(rawPtr, deleter);
        }
        
        template <typename T, typename ...ArgTypes>
        std::unique_ptr<T, DeleterType<T>> createUnique(ArgTypes&&... args) {
            T* rawPtr = alloc<T>();
            new (rawPtr) T(std::forward<ArgTypes>(args)...);
            DeleterType<T> deleter = std::bind(&Allocator::destroy<T>, this, std::placeholders::_1);
            return std::unique_ptr<T, DeleterType<T>>(rawPtr, deleter);
        }
    };
    
    class DefaultAllocator : public Allocator {
        DefaultAllocator() { };
    public:
        void* alloc(uintptr_t size, uintptr_t align) override { return SLR_memalign(size, align); }
        void free(void* ptr) override { SLR_freealign(ptr); }
        
        static DefaultAllocator &instance() {
            static DefaultAllocator s_instance;
            return s_instance;
        }
    };    
}

#endif

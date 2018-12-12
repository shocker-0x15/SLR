//
//  random_number_generator.h
//
//  Created by 渡部 心 on 2015/06/28.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_random_number_generator__
#define __SLR_random_number_generator__

#include "../defines.h"
#include "../declarations.h"

namespace SLR {
    struct SLR_API Types32bit {
        typedef int32_t Int;
        typedef uint32_t UInt;
        typedef float Float;
    };
    
    struct SLR_API Types64bit {
        typedef int64_t Int;
        typedef uint64_t UInt;
        typedef double Float;
    };
    
    template <typename TypeSet>
    class SLR_API RandomNumberGeneratorTemplate {
    public:
        virtual ~RandomNumberGeneratorTemplate() { }
        
        virtual typename TypeSet::UInt getUInt() = 0;
        virtual typename TypeSet::Float getFloat0cTo1o();
    };
    
    template <> Types32bit::Float RandomNumberGeneratorTemplate<Types32bit>::getFloat0cTo1o();
    template <> Types64bit::Float RandomNumberGeneratorTemplate<Types64bit>::getFloat0cTo1o();    
}

#endif /* __SLR_random_number_generator__ */

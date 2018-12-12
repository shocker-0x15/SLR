//
//  spectum_library.h
//
//  Created by 渡部 心 on 2015/09/20.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_spectum_library__
#define __SLR_spectum_library__

#include "../defines.h"

namespace SLR {
    namespace SpectrumLibrary {
        enum class DistributionType : uint32_t {
            Regular = 0,
            Irregular,
        };
        
        struct Data {
            DistributionType dType;
            uint32_t numSamples;
            float minLambdas;
            float maxLambdas;
            const float* lambdas;
            const float* values;
        };
        
        SLR_API bool queryIlluminantSpectrum(const std::string &name, uint32_t index, Data* data);
        SLR_API bool queryReflectanceSpectrum(const std::string &name, uint32_t index, Data* data);
        SLR_API bool queryIoRSpectrum(const std::string &name, uint32_t index, Data* data);
    }    
}

#endif /* __SLR_spectum_library__ */

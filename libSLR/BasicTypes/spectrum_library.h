//
//  spectum_library.h
//
//  Created by 渡部 心 on 2015/09/20.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef SLR_spectum_library_h
#define SLR_spectum_library_h

#include "defines.h"

namespace SLR {
    namespace SpectrumLibrary {
        enum class DistributionType : uint32_t {
            Regular = 0,
            Irregular,
        };
        
        // M. N. Polyanskiy. Refractive index database. Available at http://refractiveindex.info (accessed Feb. 29 2015).
        struct IndexOfRefraction {
            DistributionType dType;
            uint32_t numSamples;
            float minLambdas;
            float maxLambdas;
            const float* lambdas;
            const float* etas;
            const float* ks;
        };
        
        
        
        // Ciddor 1996
        extern const float Air_etas[];
        extern const IndexOfRefraction Air;
        
        extern const float Water_lambdas[];
        extern const float Water_etas[];
        extern const float Water_ks[];
        extern const IndexOfRefraction Water;
        
        
        
        // S-BSL7 (OHARA)
        extern const float Glass_BK7_etas[];
        extern const IndexOfRefraction Glass_BK7;
        
        // Peter 1923 - Diamond
        extern const float Diamond_etas[];
        extern const IndexOfRefraction Diamond;
        
        
        
        extern const float Aluminium_lambdas[];
        extern const float Aluminium_etas[];
        extern const float Aluminium_ks[];
        extern const IndexOfRefraction Aluminium;
        
        extern const float Copper_lambdas[];
        extern const float Copper_etas[];
        extern const float Copper_ks[];
        extern const IndexOfRefraction Copper;
        
        extern const float Gold_lambdas[];
        extern const float Gold_etas[];
        extern const float Gold_ks[];
        extern const IndexOfRefraction Gold;
        
        extern const float Iron_lambdas[];
        extern const float Iron_etas[];
        extern const float Iron_ks[];
        extern const IndexOfRefraction Iron;
        
        extern const float Lead_lambdas[];
        extern const float Lead_etas[];
        extern const float Lead_ks[];
        extern const IndexOfRefraction Lead;
        
        extern const float Mercury_lambdas[];
        extern const float Mercury_etas[];
        extern const float Mercury_ks[];
        extern const IndexOfRefraction Mercury;
        
        extern const float Platinum_lambdas[];
        extern const float Platinum_etas[];
        extern const float Platinum_ks[];
        extern const IndexOfRefraction Platinum;
        
        extern const float Silver_lambdas[];
        extern const float Silver_etas[];
        extern const float Silver_ks[];
        extern const IndexOfRefraction Silver;
        
        extern const float Titanium_lambdas[];
        extern const float Titanium_etas[];
        extern const float Titanium_ks[];
        extern const IndexOfRefraction Titanium;
    }    
}

#endif

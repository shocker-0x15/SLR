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
        
        namespace Illuminant {
            struct Data {
                DistributionType dType;
                uint32_t numSamples;
                float minLambdas;
                float maxLambdas;
                const float* lambdas;
                const float* values;
            };
            
            extern SLR_API const std::map<std::string, const Data&> database;
            
            extern const float D65_values[];
            extern const Data D65;
        }
        
        namespace Reflectance {
            struct Data {
                DistributionType dType;
                uint32_t numSamples;
                float minLambdas;
                float maxLambdas;
                const float* lambdas;
                const float* values;
            };
            
            extern SLR_API const std::map<std::string, const Data&> database;

            extern const float DarkSkin_values[];
            extern const Data DarkSkin;
            
            extern const float LightSkin_values[];
            extern const Data LightSkin;
            
            extern const float BlueSky_values[];
            extern const Data BlueSky;
            
            extern const float Foliage_values[];
            extern const Data Foliage;
            
            extern const float BlueFlower_values[];
            extern const Data BlueFlower;
            
            extern const float BluishGreen_values[];
            extern const Data BluishGreen;
            
            extern const float Orange_values[];
            extern const Data Orange;
            
            extern const float PurplishBlue_values[];
            extern const Data PurplishBlue;
            
            extern const float ModerateRed_values[];
            extern const Data ModerateRed;
            
            extern const float Purple_values[];
            extern const Data Purple;
            
            extern const float YellowGreen_values[];
            extern const Data YellowGreen;
            
            extern const float OrangeYellow_values[];
            extern const Data OrangeYellow;
            
            extern const float Blue_values[];
            extern const Data Blue;
            
            extern const float Green_values[];
            extern const Data Green;
            
            extern const float Red_values[];
            extern const Data Red;
            
            extern const float Yellow_values[];
            extern const Data Yellow;
            
            extern const float Magent_values[];
            extern const Data Magent;
            
            extern const float Cyan_values[];
            extern const Data Cyan;
            
            extern const float White95_values[];
            extern const Data White95;
            
            extern const float Neutral80_values[];
            extern const Data Neutral80;
            
            extern const float Neutral65_values[];
            extern const Data Neutral65;
            
            extern const float Neutral50_values[];
            extern const Data Neutral50;
            
            extern const float Neutral35_values[];
            extern const Data Neutral35;
            
            extern const float Black2_values[];
            extern const Data Black2;
        }
        
        namespace IoR {
            // M. N. Polyanskiy. Refractive index database. Available at http://refractiveindex.info (accessed Feb. 29 2015).
            struct Data {
                DistributionType dType;
                uint32_t numSamples;
                float minLambdas;
                float maxLambdas;
                const float* lambdas;
                const float* etas;
                const float* ks;
            };
            
            extern SLR_API const std::map<std::string, const Data&> database;
            
            
            
            // Ciddor 1996
            extern const float Air_etas[];
            extern const Data Air;
            
            extern const float Water_lambdas[];
            extern const float Water_etas[];
            extern const float Water_ks[];
            extern const Data Water;
            
            
            
            // S-BSL7 (OHARA)
            extern const float Glass_BK7_etas[];
            extern const Data Glass_BK7;
            
            // Peter 1923 - Diamond
            extern const float Diamond_etas[];
            extern const Data Diamond;
            
            
            
            extern const float Aluminium_lambdas[];
            extern const float Aluminium_etas[];
            extern const float Aluminium_ks[];
            extern const Data Aluminium;
            
            extern const float Copper_lambdas[];
            extern const float Copper_etas[];
            extern const float Copper_ks[];
            extern const Data Copper;
            
            extern const float Gold_lambdas[];
            extern const float Gold_etas[];
            extern const float Gold_ks[];
            extern const Data Gold;
            
            extern const float Iron_lambdas[];
            extern const float Iron_etas[];
            extern const float Iron_ks[];
            extern const Data Iron;
            
            extern const float Lead_lambdas[];
            extern const float Lead_etas[];
            extern const float Lead_ks[];
            extern const Data Lead;
            
            extern const float Mercury_lambdas[];
            extern const float Mercury_etas[];
            extern const float Mercury_ks[];
            extern const Data Mercury;
            
            extern const float Platinum_lambdas[];
            extern const float Platinum_etas[];
            extern const float Platinum_ks[];
            extern const Data Platinum;
            
            extern const float Silver_lambdas[];
            extern const float Silver_etas[];
            extern const float Silver_ks[];
            extern const Data Silver;
            
            extern const float Titanium_lambdas[];
            extern const float Titanium_etas[];
            extern const float Titanium_ks[];
            extern const Data Titanium;   
        }
    }    
}

#endif /* __SLR_spectum_library__ */

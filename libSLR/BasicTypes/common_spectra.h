//
//  common_spectra.h
//
//  Created by 渡部 心 on 2015/09/17.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef SLR_common_spectra_h
#define SLR_common_spectra_h

#include "../defines.h"

namespace SLR {
    namespace ColorChecker {
        static const float MinWavelength = 380.0f;
        static const float MaxWavelength = 730.0f;
        static const uint32_t NumSamples = 36;
        
        enum Name {
            DarkSkin = 0,
            LightSkin,
            BlueSky,
            Foliage,
            BlueFlower,
            BluishGreen,
            Orange,
            PurplishBlue,
            ModerateRed,
            Purple,
            YellowGreen,
            OrangeYellow,
            Blue,
            Green,
            Red,
            Yellow,
            Magent,
            Cyan,
            White95,
            Neutral80,
            Neutral65,
            Neutral50,
            Neutral35,
            Black2
        };
        
        extern SLR_API const float Spectra[24][NumSamples];
    }
    
    namespace StandardIlluminant {
        static const float MinWavelength = 300.0f;
        static const float MaxWavelength = 830.0f;
        static const uint32_t NumSamples = 531;
        
        extern SLR_API const float D65[NumSamples];
    }
}

#endif

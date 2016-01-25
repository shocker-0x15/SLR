//
//  RGBTypes.cpp
//
//  Created by 渡部 心 on 2016/01/26.
//  Copyright © 2016年 渡部 心. All rights reserved.
//

#include "RGBTypes.h"

namespace SLR {
    template struct RGBSamplesTemplate<float>;
    template struct RGBSamplesTemplate<double>;
    
    template struct RGBTemplate<float>;
    template struct RGBTemplate<double>;
    
    template <typename RealType>
    RGBTemplate<RealType> pow(const RGBTemplate<RealType> &s, RealType p) {
        return RGBTemplate<RealType>(std::pow(s.r, p), std::pow(s.g, p), std::pow(s.b, p));
    }
    template RGBTemplate<float> pow(const RGBTemplate<float> &s, float p);
    template RGBTemplate<double> pow(const RGBTemplate<double> &s, double p);
    
    template <typename RealType>
    RGBTemplate<RealType> exp(const RGBTemplate<RealType> &s) {
        return RGBTemplate<RealType>(std::exp(s.r), std::exp(s.g), std::exp(s.b));
    }
    template RGBTemplate<float> exp(const RGBTemplate<float> &s);
    template RGBTemplate<double> exp(const RGBTemplate<double> &s);
    
    template <typename RealType>
    RGBTemplate<RealType> inverseGammaCorrection(const RGBTemplate<RealType> &s, RealType gamma) {
        return RGBTemplate<RealType>(std::pow(s.r, gamma), std::pow(s.g, gamma), std::pow(s.b, gamma));
    }
    template RGBTemplate<float> inverseGammaCorrection(const RGBTemplate<float> &s, float gamma);
    template RGBTemplate<double> inverseGammaCorrection(const RGBTemplate<double> &s, double gamma);
    
    template struct RGBStorageTemplate<float>;
    template struct RGBStorageTemplate<double>;
}

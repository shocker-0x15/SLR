//
//  RGBTypes.cpp
//
//  Created by 渡部 心 on 2016/01/26.
//  Copyright c 2016年 渡部 心. All rights reserved.
//

#include "RGBTypes.h"

namespace SLR {
    template struct SLR_API RGBSamplesTemplate<float>;
//    template SLR_API const uint32_t RGBSamplesTemplate<float>::NumComponents;

    template struct SLR_API RGBSamplesTemplate<double>;
//    template SLR_API const uint32_t RGBSamplesTemplate<double>::NumComponents;


    template struct SLR_API RGBTemplate<float>;
//    template SLR_API const uint32_t RGBTemplate<float>::NumComponents;
//    template SLR_API const RGBTemplate<float> RGBTemplate<float>::Zero;
//    template SLR_API const RGBTemplate<float> RGBTemplate<float>::One;
//    template SLR_API const RGBTemplate<float> RGBTemplate<float>::Inf;
//    template SLR_API const RGBTemplate<float> RGBTemplate<float>::NaN;

    template struct SLR_API RGBTemplate<double>;
//    template SLR_API const uint32_t RGBTemplate<double>::NumComponents;
//    template SLR_API const RGBTemplate<double> RGBTemplate<double>::Zero;
//    template SLR_API const RGBTemplate<double> RGBTemplate<double>::One;
//    template SLR_API const RGBTemplate<double> RGBTemplate<double>::Inf;
//    template SLR_API const RGBTemplate<double> RGBTemplate<double>::NaN;

    
    template <typename RealType>
    RGBTemplate<RealType> pow(const RGBTemplate<RealType> &s, RealType p) {
        return RGBTemplate<RealType>(std::pow(s.r, p), std::pow(s.g, p), std::pow(s.b, p));
    }
    template SLR_API RGBTemplate<float> pow(const RGBTemplate<float> &s, float p);
    template SLR_API RGBTemplate<double> pow(const RGBTemplate<double> &s, double p);
    
    template <typename RealType>
    RGBTemplate<RealType> exp(const RGBTemplate<RealType> &s) {
        return RGBTemplate<RealType>(std::exp(s.r), std::exp(s.g), std::exp(s.b));
    }
    template SLR_API RGBTemplate<float> exp(const RGBTemplate<float> &s);
    template SLR_API RGBTemplate<double> exp(const RGBTemplate<double> &s);
    
    template <typename RealType>
    RGBTemplate<RealType> inverseGammaCorrection(const RGBTemplate<RealType> &s, RealType gamma) {
        return RGBTemplate<RealType>(std::pow(s.r, gamma), std::pow(s.g, gamma), std::pow(s.b, gamma));
    }
    template SLR_API RGBTemplate<float> inverseGammaCorrection(const RGBTemplate<float> &s, float gamma);
    template SLR_API RGBTemplate<double> inverseGammaCorrection(const RGBTemplate<double> &s, double gamma);
    
    template struct SLR_API RGBStorageTemplate<float>;
    template struct SLR_API RGBStorageTemplate<double>;
}

//
//  SpectrumType.cpp
//
//  Created by 渡部 心 on 2016/01/26.
//  Copyright c 2016年 渡部 心. All rights reserved.
//

#include "SpectrumTypes.h"

namespace SLR {
    template struct SLR_API WavelengthSamplesTemplate<float, NumSpectralSamples>;
//    template SLR_API const uint32_t WavelengthSamplesTemplate<float, NumSpectralSamples>::NumComponents;

    template struct SLR_API WavelengthSamplesTemplate<double, NumSpectralSamples>;
//    template SLR_API const uint32_t WavelengthSamplesTemplate<double, NumSpectralSamples>::NumComponents;

    
    template struct SLR_API ContinuousSpectrumTemplate<float, NumSpectralSamples>;
    template struct SLR_API ContinuousSpectrumTemplate<double, NumSpectralSamples>;

    
    template struct SLR_API RegularContinuousSpectrumTemplate<float, NumSpectralSamples>;
    template struct SLR_API RegularContinuousSpectrumTemplate<double, NumSpectralSamples>;

    
    template struct SLR_API IrregularContinuousSpectrumTemplate<float, NumSpectralSamples>;
    template struct SLR_API IrregularContinuousSpectrumTemplate<double, NumSpectralSamples>;

    
    template struct SLR_API UpsampledContinuousSpectrumTemplate<float, NumSpectralSamples>;
    template struct SLR_API UpsampledContinuousSpectrumTemplate<double, NumSpectralSamples>;
    

    template struct SLR_API SampledSpectrumTemplate<float, NumSpectralSamples>;
//    template SLR_API const uint32_t SampledSpectrumTemplate<float, NumSpectralSamples>::NumComponents;
//    template SLR_API const SampledSpectrumTemplate<float, NumSpectralSamples> SampledSpectrumTemplate<float, NumSpectralSamples>::Zero;
//    template SLR_API const SampledSpectrumTemplate<float, NumSpectralSamples> SampledSpectrumTemplate<float, NumSpectralSamples>::One;
//    template SLR_API const SampledSpectrumTemplate<float, NumSpectralSamples> SampledSpectrumTemplate<float, NumSpectralSamples>::Inf;
//    template SLR_API const SampledSpectrumTemplate<float, NumSpectralSamples> SampledSpectrumTemplate<float, NumSpectralSamples>::NaN;

    template struct SLR_API SampledSpectrumTemplate<double, NumSpectralSamples>;
//    template SLR_API const uint32_t SampledSpectrumTemplate<double, NumSpectralSamples>::NumComponents;
//    template SLR_API const SampledSpectrumTemplate<double, NumSpectralSamples> SampledSpectrumTemplate<double, NumSpectralSamples>::Zero;
//    template SLR_API const SampledSpectrumTemplate<double, NumSpectralSamples> SampledSpectrumTemplate<double, NumSpectralSamples>::One;
//    template SLR_API const SampledSpectrumTemplate<double, NumSpectralSamples> SampledSpectrumTemplate<double, NumSpectralSamples>::Inf;
//    template SLR_API const SampledSpectrumTemplate<double, NumSpectralSamples> SampledSpectrumTemplate<double, NumSpectralSamples>::NaN;
    
    template <typename RealType, uint32_t N>
    SampledSpectrumTemplate<RealType, N> sqrt(const SampledSpectrumTemplate<RealType, N> &value) {
        SampledSpectrumTemplate<RealType, N> ret;
        for (int i = 0; i < N; ++i)
            ret[i] = std::sqrt(value[i]);
        return ret;
    }
    template SLR_API SampledSpectrumTemplate<float, NumSpectralSamples> sqrt(const SampledSpectrumTemplate<float, NumSpectralSamples> &value);
    template SLR_API SampledSpectrumTemplate<double, NumSpectralSamples> sqrt(const SampledSpectrumTemplate<double, NumSpectralSamples> &value);
    
    template <typename RealType, uint32_t N>
    SampledSpectrumTemplate<RealType, N> exp(const SampledSpectrumTemplate<RealType, N> &value) {
        SampledSpectrumTemplate<RealType, N> ret;
        for (int i = 0; i < N; ++i)
            ret[i] = std::exp(value[i]);
        return ret;
    }
    template SLR_API SampledSpectrumTemplate<float, NumSpectralSamples> exp(const SampledSpectrumTemplate<float, NumSpectralSamples> &value);
    template SLR_API SampledSpectrumTemplate<double, NumSpectralSamples> exp(const SampledSpectrumTemplate<double, NumSpectralSamples> &value);

    
    template struct DiscretizedSpectrumTemplate<float, NumStrataForStorage>;
//    template SLR_API const uint32_t DiscretizedSpectrumTemplate<float, NumStrataForStorage>::NumStrata;
//    template SLR_API std::unique_ptr<float[]> DiscretizedSpectrumTemplate<float, NumStrataForStorage>::xbar;
//    template SLR_API std::unique_ptr<float[]> DiscretizedSpectrumTemplate<float, NumStrataForStorage>::ybar;
//    template SLR_API std::unique_ptr<float[]> DiscretizedSpectrumTemplate<float, NumStrataForStorage>::zbar;
//    template SLR_API float DiscretizedSpectrumTemplate<float, NumStrataForStorage>::integralCMF;
//    template SLR_API const DiscretizedSpectrumTemplate<float, NumStrataForStorage> DiscretizedSpectrumTemplate<float, NumStrataForStorage>::Zero;
//    template SLR_API const DiscretizedSpectrumTemplate<float, NumStrataForStorage> DiscretizedSpectrumTemplate<float, NumStrataForStorage>::One;
//    template SLR_API const DiscretizedSpectrumTemplate<float, NumStrataForStorage> DiscretizedSpectrumTemplate<float, NumStrataForStorage>::Inf;
//    template SLR_API const DiscretizedSpectrumTemplate<float, NumStrataForStorage> DiscretizedSpectrumTemplate<float, NumStrataForStorage>::NaN;

    template struct DiscretizedSpectrumTemplate<double, NumStrataForStorage>;
//    template SLR_API const uint32_t DiscretizedSpectrumTemplate<double, NumStrataForStorage>::NumStrata;
//    template SLR_API std::unique_ptr<double[]> DiscretizedSpectrumTemplate<double, NumStrataForStorage>::xbar;
//    template SLR_API std::unique_ptr<double[]> DiscretizedSpectrumTemplate<double, NumStrataForStorage>::ybar;
//    template SLR_API std::unique_ptr<double[]> DiscretizedSpectrumTemplate<double, NumStrataForStorage>::zbar;
//    template SLR_API double DiscretizedSpectrumTemplate<double, NumStrataForStorage>::integralCMF;
//    template SLR_API const DiscretizedSpectrumTemplate<double, NumStrataForStorage> DiscretizedSpectrumTemplate<double, NumStrataForStorage>::Zero;
//    template SLR_API const DiscretizedSpectrumTemplate<double, NumStrataForStorage> DiscretizedSpectrumTemplate<double, NumStrataForStorage>::One;
//    template SLR_API const DiscretizedSpectrumTemplate<double, NumStrataForStorage> DiscretizedSpectrumTemplate<double, NumStrataForStorage>::Inf;
//    template SLR_API const DiscretizedSpectrumTemplate<double, NumStrataForStorage> DiscretizedSpectrumTemplate<double, NumStrataForStorage>::NaN;
    
    template struct SLR_API SpectrumStorageTemplate<float, NumStrataForStorage>;
    template struct SLR_API SpectrumStorageTemplate<double, NumStrataForStorage>;
}

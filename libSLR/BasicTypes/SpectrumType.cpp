//
//  SpectrumType.cpp
//
//  Created by 渡部 心 on 2016/01/26.
//  Copyright © 2016年 渡部 心. All rights reserved.
//

#include "SpectrumTypes.h"

namespace SLR {
    template struct WavelengthSamplesTemplate<float, NumSpectralSamples>;
    template struct WavelengthSamplesTemplate<double, NumSpectralSamples>;
    
    template struct ContinuousSpectrumTemplate<float, NumSpectralSamples>;
    template struct ContinuousSpectrumTemplate<double, NumSpectralSamples>;
    
    template struct RegularContinuousSpectrumTemplate<float, NumSpectralSamples>;
    template struct RegularContinuousSpectrumTemplate<double, NumSpectralSamples>;
    
    template struct IrregularContinuousSpectrumTemplate<float, NumSpectralSamples>;
    template struct IrregularContinuousSpectrumTemplate<double, NumSpectralSamples>;
    
    template struct UpsampledContinuousSpectrumTemplate<float, NumSpectralSamples>;
    template struct UpsampledContinuousSpectrumTemplate<double, NumSpectralSamples>;
    
    template struct SampledSpectrumTemplate<float, NumSpectralSamples>;
    template struct SampledSpectrumTemplate<double, NumSpectralSamples>;
    
    template struct DiscretizedSpectrumTemplate<float, NumStrataForStorage>;
    template struct DiscretizedSpectrumTemplate<double, NumStrataForStorage>;
    
    template struct SpectrumStorageTemplate<float, NumStrataForStorage>;
    template struct SpectrumStorageTemplate<double, NumStrataForStorage>;
}

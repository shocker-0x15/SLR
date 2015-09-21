//
//  directional_distribution_functions.cpp
//
//  Created by 渡部 心 on 2015/06/28.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "directional_distribution_functions.h"

SampledSpectrum FresnelNoOp::evaluate(float cosEnter) const {
    return SampledSpectrum::One;
}

SampledSpectrum FresnelConductor::evaluate(float cosEnter) const {
    cosEnter = std::fabs(cosEnter);
    float cosEnter2 = cosEnter * cosEnter;
    SampledSpectrum _2EtaCosEnter = 2.0f * m_eta * cosEnter;
    SampledSpectrum tmp_f = m_eta * m_eta + m_k * m_k;
    SampledSpectrum tmp = tmp_f * cosEnter2;
    SampledSpectrum Rparl2 = (tmp - _2EtaCosEnter + 1) / (tmp + _2EtaCosEnter + 1);
    SampledSpectrum Rperp2 = (tmp_f - _2EtaCosEnter + cosEnter2) / (tmp_f + _2EtaCosEnter + cosEnter2);
    return (Rparl2 + Rperp2) / 2.0f;
}

SampledSpectrum FresnelDielectric::evaluate(float cosEnter) const {
    cosEnter = std::clamp(cosEnter, -1.0f, 1.0f);
    
    bool entering = cosEnter > 0.0f;
    const SampledSpectrum &eEnter = entering ? m_etaExt : m_etaInt;
    const SampledSpectrum &eExit = entering ? m_etaInt : m_etaExt;
    
    SampledSpectrum sinExit = eEnter / eExit * std::sqrt(std::fmax(0.0f, 1.0f - cosEnter * cosEnter));
    SampledSpectrum ret = SampledSpectrum::Zero;
    cosEnter = std::fabs(cosEnter);
    for (int i = 0; i < SampledSpectrum::NumComponents; ++i) {
        if (sinExit[i] >= 1.0f) {
            ret[i] = 1.0f;
        }
        else {
            float cosExit = std::sqrt(std::fmax(0.0f, 1.0f - sinExit[i] * sinExit[i]));
            ret[i] = evalF(eEnter[i], eExit[i], cosEnter, cosExit);
        }
    }
    return ret;
}

float FresnelDielectric::evalF(float etaEnter, float etaExit, float cosEnter, float cosExit) {
    float Rparl = ((etaExit * cosEnter) - (etaEnter * cosExit)) / ((etaExit * cosEnter) + (etaEnter * cosExit));
    float Rperp = ((etaEnter * cosEnter) - (etaExit * cosExit)) / ((etaEnter * cosEnter) + (etaExit * cosExit));
    return (Rparl * Rparl + Rperp * Rperp) / 2.0f;
}

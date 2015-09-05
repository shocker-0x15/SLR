//
//  directional_distribution_functions.cpp
//
//  Created by 渡部 心 on 2015/06/28.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "directional_distribution_functions.h"

Spectrum FresnelNoOp::evaluate(float cosEnter) const {
    return Spectrum::One;
}

Spectrum FresnelConductor::evaluate(float cosEnter) const {
    cosEnter = std::fabs(cosEnter);
    float cosEnter2 = cosEnter * cosEnter;
    Spectrum _2EtaCosEnter = 2.0f * m_eta * cosEnter;
    Spectrum tmp_f = m_eta * m_eta + m_k * m_k;
    Spectrum tmp = tmp_f * cosEnter2;
    Spectrum Rparl2 = (tmp - _2EtaCosEnter + 1) / (tmp + _2EtaCosEnter + 1);
    Spectrum Rperp2 = (tmp_f - _2EtaCosEnter + cosEnter2) / (tmp_f + _2EtaCosEnter + cosEnter2);
    return (Rparl2 + Rperp2) / 2.0f;
}

Spectrum FresnelDielectric::evaluate(float cosEnter) const {
    cosEnter = std::clamp(cosEnter, -1.0f, 1.0f);
    
    bool entering = cosEnter > 0.0f;
    float eEnter = m_etaExt;
    float eExit = m_etaInt;
    if (!entering)
        std::swap(eEnter, eExit);
    
    float sinExit = eEnter / eExit * std::sqrt(std::fmax(0.0f, 1.0f - cosEnter * cosEnter));
    if (sinExit >= 1.0f) {
        return Spectrum::One;
    }
    else {
        float cosExit = std::sqrt(std::fmax(0.0f, 1.0f - sinExit * sinExit));
        cosEnter = std::fabs(cosEnter);
        
        Spectrum Rparl = ((eExit * cosEnter) - (eEnter * cosExit)) / ((eExit * cosEnter) + (eEnter * cosExit));
        Spectrum Rperp = ((eEnter * cosEnter) - (eExit * cosExit)) / ((eEnter * cosEnter) + (eExit * cosExit));
        return (Rparl * Rparl + Rperp * Rperp) / 2.0f;
    }
}

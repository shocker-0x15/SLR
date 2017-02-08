//
//  basic_phase_functions.cpp
//
//  Created by 渡部 心 on 2016/09/04.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "basic_phase_functions.h"
#include "../Core/distributions.h"

namespace SLR {
    SampledSpectrum IsotropicPhaseFunction::sample(const PFQuery &query, const PFSample &smp, PFQueryResult* result) const {
        result->dirLocal = uniformSampleSphere(smp.uDir[0], smp.uDir[1]);
        result->dirPDF = 1 / (4 * M_PI);
        result->sampledType = directionType();
        return SampledSpectrum(1 / (4 * M_PI));
    }
    
    SampledSpectrum IsotropicPhaseFunction::evaluate(const PFQuery &query, const Vector3D &dirIn) const {
        return SampledSpectrum(1 / (4 * M_PI));
    }
    
    float IsotropicPhaseFunction::evaluatePDF(const PFQuery &query, const Vector3D &dirIn) const {
        return 1 / (4 * M_PI);
    }
    
    
    
    SampledSpectrum HenyeyGreensteinPhaseFunction::sample(const PFQuery &query, const PFSample &smp, PFQueryResult* result) const {
        float sqTerm = (1 - m_g * m_g) / (1 - m_g + 2 * m_g * smp.uDir[1]);
        float cosTheta = std::clamp(0.5f / m_g * (1 + m_g * m_g - sqTerm * sqTerm), -1.0f, 1.0f);
        float phi = 2 * M_PI * smp.uDir[0];
        
        float value = (1 - m_g * m_g) / (4 * M_PI * std::pow(1 + m_g * m_g - 2 * m_g * cosTheta, 1.5f));
        float sinTheta = std::sqrt(1 - cosTheta * cosTheta);
        result->dirLocal = Vector3D(std::cos(phi) * sinTheta, std::sin(phi) * sinTheta, cosTheta);
        result->dirPDF = value;
        result->sampledType = directionType();
        return SampledSpectrum(value);
    }
    
    SampledSpectrum HenyeyGreensteinPhaseFunction::evaluate(const PFQuery &query, const Vector3D &dirIn) const {
        float cosTheta = dirIn.z;
        float value = (1 - m_g * m_g) / (4 * M_PI * std::pow(1 + m_g * m_g - 2 * m_g * cosTheta, 1.5f));
        return SampledSpectrum(value);
    }
    
    float HenyeyGreensteinPhaseFunction::evaluatePDF(const PFQuery &query, const Vector3D &dirIn) const {
        float cosTheta = dirIn.z;
        return (1 - m_g * m_g) / (4 * M_PI * std::pow(1 + m_g * m_g - 2 * m_g * cosTheta, 1.5f));
    }
    
    
    
    SampledSpectrum SchlickPhaseFunction::sample(const PFQuery &query, const PFSample &smp, PFQueryResult* result) const {
        float cosTheta = std::clamp((2 * smp.uDir[1] + m_k - 1) / (2 * m_k * smp.uDir[1] - m_k + 1), -1.0f, 1.0f);
        float phi = 2 * M_PI * smp.uDir[0];
        
        float dTerm = (1 + m_k * cosTheta);
        float value = (1 - m_k * m_k) / (4 * M_PI * dTerm * dTerm);
        float sinTheta = std::sqrt(1 - cosTheta * cosTheta);
        result->dirLocal = Vector3D(std::cos(phi) * sinTheta, std::sin(phi) * sinTheta, cosTheta);
        result->dirPDF = value;
        result->sampledType = directionType();
        return SampledSpectrum(value);
    }
    
    SampledSpectrum SchlickPhaseFunction::evaluate(const PFQuery &query, const Vector3D &dirIn) const {
        float cosTheta = dirIn.z;
        float dTerm = (1 + m_k * cosTheta);
        float value = (1 - m_k * m_k) / (4 * M_PI * dTerm * dTerm);
        return SampledSpectrum(value);
    }
    
    float SchlickPhaseFunction::evaluatePDF(const PFQuery &query, const Vector3D &dirIn) const {
        float cosTheta = dirIn.z;
        float dTerm = (1 + m_k * cosTheta);
        return (1 - m_k * m_k) / (4 * M_PI * dTerm * dTerm);
    }
}

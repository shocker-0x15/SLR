//
//  basic_phase_functions.h
//
//  Created by 渡部 心 on 2016/09/04.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef basic_phase_functions_h
#define basic_phase_functions_h

#include "../defines.h"
#include "../references.h"
#include "../Core/directional_distribution_functions.h"

namespace SLR {
    class SLR_API IsotropicPhaseFunction : public PhaseFunction {
    public:
        SampledSpectrum sample(const Vector3D &dirOut, const PFSample &smp, PFQueryResult* result) const override;
        SampledSpectrum evaluate(const Vector3D &dirOut, const Vector3D &dirIn) const override;
        float evaluatePDF(const Vector3D &dirOut, const Vector3D &dirIn) const override;
    };
    
    
    
    class SLR_API HenyeyGreensteinPhaseFunction : public PhaseFunction {
        float m_g;
    public:
        HenyeyGreensteinPhaseFunction(float g) : m_g(g) { SLRAssert(m_g != 0.0f, "g must not be equal to zero."); }
        
        SampledSpectrum sample(const Vector3D &dirOut, const PFSample &smp, PFQueryResult* result) const override;
        SampledSpectrum evaluate(const Vector3D &dirOut, const Vector3D &dirIn) const override;
        float evaluatePDF(const Vector3D &dirOut, const Vector3D &dirIn) const override;
    };
    
    
    
    class SLR_API SchlickPhaseFunction : public PhaseFunction {
        float m_k;
    public:
        SchlickPhaseFunction(float k) : m_k(k) { SLRAssert(k > -1.0f && k < 1.0f, "k must be in the range(-1, 1). %g", m_k); }
        
        SampledSpectrum sample(const Vector3D &dirOut, const PFSample &smp, PFQueryResult* result) const override;
        SampledSpectrum evaluate(const Vector3D &dirOut, const Vector3D &dirIn) const override;
        float evaluatePDF(const Vector3D &dirOut, const Vector3D &dirIn) const override;
    };
}

#endif /* basic_phase_functions_h */

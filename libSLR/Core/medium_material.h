//
//  medium_material.h
//
//  Created by 渡部 心 on 2017/01/04.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef __SLR_medium_material__
#define __SLR_medium_material__

#include "../defines.h"
#include "../references.h"
#include "../BasicTypes/RGBTypes.h"
#include "../BasicTypes/SpectrumTypes.h"

namespace SLR {
    class SLR_API MediumMaterial {
    public:
        MediumMaterial() { }
        virtual ~MediumMaterial() { }
        
        virtual PhaseFunction* getPhaseFunction(const MediumPoint &medPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const = 0;
        virtual SampledSpectrum emittance(const MediumPoint &medPt, const WavelengthSamples &wls) const { return SampledSpectrum::Zero; }
        virtual bool isEmitting() const { return false; }
        virtual EDF* getEDF(const MediumPoint &medPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const { SLRAssert(false, "Not implemented."); return nullptr; }
    };
    
    
    
    class SLR_API EmitterMediumProperty {
    public:
        virtual ~EmitterMediumProperty() { }
        
        virtual SampledSpectrum emittance(const MediumPoint &medPt, const WavelengthSamples &wls) const = 0;
        virtual EDF* getEDF(const MediumPoint &medPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const = 0;
    };
    
    
    
    class SLR_API EmitterMediumMaterial : public MediumMaterial {
        const MediumMaterial* m_mat;
        const EmitterMediumProperty* m_emit;
    public:
        EmitterMediumMaterial(const MediumMaterial* mat, const EmitterMediumProperty* emit) :
        m_mat(mat), m_emit(emit) {
            if (m_mat)
                SLRAssert(!m_mat->isEmitting(), "EmitterMediumMaterial can not have an emitting MediumMaterial.");
        }
        
        PhaseFunction* getPhaseFunction(const MediumPoint &medPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const override { return m_mat->getPhaseFunction(medPt, wls, mem); }
        SampledSpectrum emittance(const MediumPoint &medPt, const WavelengthSamples &wls) const override { return m_emit->emittance(medPt, wls); }
        bool isEmitting() const override { return true; }
        EDF* getEDF(const MediumPoint &medPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const override { return m_emit->getEDF(medPt, wls, mem); }
    };
}

#endif /* __SLR_medium_material__ */

//
//  geometry.cpp
//
//  Created by 渡部 心 on 2015/08/01.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "geometry.h"

#include "../MemoryAllocators/ArenaAllocator.h"
#include "transform.h"
#include "surface_object.h"
#include "medium_object.h"

namespace SLR {
    void SurfaceInteraction::calculateSurfacePoint(SurfacePoint* surfPt) const {
        m_obj->calculateSurfacePoint(*this, surfPt);
    }
    
    InteractionPoint* SurfaceInteraction::createInteractionPoint(ArenaAllocator &mem) const {
        SurfacePoint* surfPt = mem.create<SurfacePoint>();
        calculateSurfacePoint(surfPt);
        return surfPt;
    }
    
    
    
    void MediumInteraction::calculateMediumPoint(SLR::MediumPoint *medPt) const {
        m_obj->calculateMediumPoint(*this, medPt);
    }
    
    InteractionPoint* MediumInteraction::createInteractionPoint(ArenaAllocator &mem) const {
        MediumPoint* medPt = mem.create<MediumPoint>();
        calculateMediumPoint(medPt);
        return medPt;
    }
    
    
    
    void InteractionPoint::applyTransform(const SLR::StaticTransform &transform) {
        m_p = transform * m_p;
        m_shadingFrame.x = normalize(transform * m_shadingFrame.x);
        m_shadingFrame.y = normalize(transform * m_shadingFrame.y);
        m_shadingFrame.z = normalize(transform * m_shadingFrame.z);
    }
    
    
    
    SampledSpectrum SurfacePoint::emittance(const WavelengthSamples &wls) const {
        return m_obj->emittance(*this, wls);
    }
    
    float SurfacePoint::evaluateAreaPDF() const {
        return m_obj->evaluateAreaPDF(*this);
    }
    
    BSDF* SurfacePoint::createBSDF(const WavelengthSamples &wls, ArenaAllocator &mem) const {
        return m_obj->createBSDF(*this, wls, mem);
    }
    
    bool SurfacePoint::isEmitting() const {
        return m_obj->isEmitting();
    }
    
    EDF* SurfacePoint::createEDF(const WavelengthSamples &wls, ArenaAllocator &mem) const {
        return m_obj->createEDF(*this, wls, mem);
    }
    
    ABDFQuery* SurfacePoint::createABDFQuery(const Vector3D &dirLocal, int16_t selectedWL, DirectionType filter, bool reqRev, bool adjoint, ArenaAllocator &mem) const {
        return mem.create<BSDFQuery>(dirLocal, toLocal(m_gNormal), selectedWL, filter, reqRev, adjoint);
    }
    
    void SurfacePoint::applyTransform(const SLR::StaticTransform &transform) {
        InteractionPoint::applyTransform(transform);
        m_gNormal = normalize(transform * m_gNormal);
        m_texCoord0Dir = normalize(transform * m_texCoord0Dir);
    }
    
    
    
    SampledSpectrum MediumPoint::emittance(const WavelengthSamples &wls) const {
        SLRAssert_NotImplemented();
//        return m_obj->emittance(*this, wls) / evaluateExtinctionCoefficient(); // need to consider dividing by zero.
        return SampledSpectrum::Zero;
    }
    
    float MediumPoint::evaluateVolumePDF() const {
        SLRAssert_NotImplemented();
        return 0.0f;
    }
    
    BSDF* MediumPoint::createPhaseFunction(const WavelengthSamples &wls, ArenaAllocator &mem) const {
        SLRAssert_NotImplemented();
        return nullptr;
    }
    
    bool MediumPoint::isEmitting() const {
        return m_obj->isEmitting();
    }
    
    EDF* MediumPoint::createEDF(const WavelengthSamples &wls, SLR::ArenaAllocator &mem) const {
//        return m_obj->createEDF(*this, wls, mem);
        SLRAssert_NotImplemented();
        return nullptr;
    }
    
    ABDFQuery* MediumPoint::createABDFQuery(const Vector3D &dirLocal, int16_t selectedWL, DirectionType filter, bool reqRev, bool adjoint, ArenaAllocator &mem) const {
        return mem.create<VolumetricBSDFQuery>(dirLocal, selectedWL, filter, reqRev);
    }
    
    AbstractBDF* MediumPoint::createAbstractBDF(const WavelengthSamples &wls, SLR::ArenaAllocator &mem) const {
        return m_obj->createAbstractBDF(*this, wls, mem);
    }
    
    SampledSpectrum MediumPoint::evaluateExtinctionCoefficient(const WavelengthSamples &wls) const {
        return m_obj->evaluateExtinctionCoefficient(*this, wls);
    }
    
    void MediumPoint::applyTransform(const SLR::StaticTransform &transform) {
        InteractionPoint::applyTransform(transform);
    }
}

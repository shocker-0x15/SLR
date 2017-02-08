//
//  GridMedium.cpp
//
//  Created by 渡部 心 on 2017/01/04.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "GridMedium.h"
#include "light_path_samplers.h"

namespace SLR {
    bool GridMedium::subdivide(Allocator* mem, Medium** fragments, uint32_t* numFragments) const {
        SLRAssert_NotImplemented();
        return true;
    }
    
    bool GridMedium::interact(const Ray &ray, float distanceLimit, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                              MediumInteraction *mi, SampledSpectrum *medThroughput, bool* singleWavelength) const {
        FreePathSampler &sampler = pathSampler.getFreePathSampler();
        Point3D queryPoint = ray.org + ray.distMin * ray.dir;
        FloatSum sampledDistance = ray.distMin;
        
        float majorant = majorantExtinctionCoefficient();
        *singleWavelength = true;
        sampledDistance += -std::log(sampler.getSample()) / majorant;
        while (sampledDistance < distanceLimit) {
            queryPoint = ray.org + sampledDistance * ray.dir;
            SampledSpectrum extCoeff = extinctionCoefficient(queryPoint, wls);
            float probRealCollision = extCoeff[wls.selectedLambda] / majorant;
            if (sampler.getSample() < probRealCollision) {
                Point3D param;
                m_region.localCoordinates(queryPoint, &param);
                SampledSpectrum scatCoeff;
                SLRAssert_NotImplemented();
                *mi = MediumInteraction(ray.time, sampledDistance, queryPoint,
                                        scatCoeff, extCoeff, normalize(ray.dir), param.x, param.y, param.z);
                *medThroughput = SampledSpectrum::Zero;
                (*medThroughput)[wls.selectedLambda] = 1.0f / extCoeff[wls.selectedLambda];
                return true;
            }
            sampledDistance += -std::log(sampler.getSample()) / majorant;
            
            // TODO: handle out of boundary.
        }
        
        *medThroughput = SampledSpectrum::Zero;
        (*medThroughput)[wls.selectedLambda] = 1.0f;
        return false;
    }
    
    SampledSpectrum GridMedium::evaluateTransmittance(Ray &ray, float distanceLimit, const WavelengthSamples &wls, SLR::LightPathSampler &pathSampler, bool *singleWavelength) const {
        FreePathSampler sampler = pathSampler.getFreePathSampler();
        Point3D queryPoint = ray.org + ray.distMin * ray.dir;
        FloatSum sampledDistance = ray.distMin;
        
        SampledSpectrum transmittance = SampledSpectrum::Zero;
        transmittance[wls.selectedLambda] = 1.0f;
        
        float majorant = majorantExtinctionCoefficient();
        *singleWavelength = true;
        sampledDistance += -std::log(sampler.getSample()) / majorant;
        while (sampledDistance < distanceLimit) {
            queryPoint = ray.org + sampledDistance * ray.dir;
            SampledSpectrum extCoeff = extinctionCoefficient(queryPoint, wls);
            float probRealCollision = extCoeff[wls.selectedLambda] / majorant;
            transmittance[wls.selectedLambda] *= (1.0f - probRealCollision);
            sampledDistance += -std::log(sampler.getSample()) / majorant;
        }
        
        return transmittance;
    }
    
    SampledSpectrum GridMedium::extinctionCoefficient(const Point3D &p, const WavelengthSamples &wls) const {
        if (!m_region.contains(p))
            return SampledSpectrum::Zero;
        Point3D param;
        m_region.localCoordinates(p, &param);
        
        SLRAssert_NotImplemented();
        
        return SampledSpectrum::Zero;
    }
    
    void GridMedium::getMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const {
        SLRAssert_NotImplemented();
    }
    
    void GridMedium::sample(float u0, float u1, float u2, MediumPoint *medPt, float *volumePDF) const {
        SLRAssert_NotImplemented();
    }
    
    float GridMedium::evaluateVolumePDF(const MediumPoint &medPt) const {
        return 1.0f / volume();
    }
    
    
    
    bool SubGridMedium::interact(const Ray &ray, float distanceLimit, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                                 MediumInteraction *mi, SampledSpectrum *medThroughput, bool* singleWavelength) const {
        SLRAssert_NotImplemented();
        return false;
    }
    
    SampledSpectrum SubGridMedium::evaluateTransmittance(Ray &ray, float distanceLimit, const WavelengthSamples &wls, SLR::LightPathSampler &pathSampler, bool *singleWavelength) const {
        SLRAssert_NotImplemented();
        return SampledSpectrum::Zero;
    }
    
    SampledSpectrum SubGridMedium::extinctionCoefficient(const Point3D &p, const WavelengthSamples &wls) const {
        SLRAssert_NotImplemented();
        return SampledSpectrum::Zero;
    }
    
    void SubGridMedium::getMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const {
        SLRAssert_NotImplemented();
    }
    
    void SubGridMedium::sample(float u0, float u1, float u2, MediumPoint *medPt, float *volumePDF) const {
        SLRAssert_NotImplemented();
    }
    
    float SubGridMedium::evaluateVolumePDF(const MediumPoint &medPt) const {
        return 1.0f / volume();
    }
}

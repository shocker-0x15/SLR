//
//  GridMedium.cpp
//
//  Created by 渡部 心 on 2017/01/04.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "GridMedium.h"
#include "light_path_samplers.h"

namespace SLR {
    bool GridMediumDistribution::subdivide(Allocator* mem, MediumDistribution** fragments, uint32_t* numFragments) const {
        SLRAssert_NotImplemented();
        return true;
    }
    
    bool GridMediumDistribution::interact(const Ray &ray, float distanceLimit, const WavelengthSamples &wls, LightPathSampler &pathSampler,
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
    
    SampledSpectrum GridMediumDistribution::evaluateTransmittance(Ray &ray, float distanceLimit, const WavelengthSamples &wls, SLR::LightPathSampler &pathSampler, bool *singleWavelength) const {
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
    
    SampledSpectrum GridMediumDistribution::extinctionCoefficient(const Point3D &p, const WavelengthSamples &wls) const {
        if (!m_region.contains(p))
            return SampledSpectrum::Zero;
        Point3D param;
        m_region.localCoordinates(p, &param);
        
        SLRAssert_NotImplemented();
        
        return SampledSpectrum::Zero;
    }
    
    void GridMediumDistribution::getMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const {
        SLRAssert_NotImplemented();
    }
    
    void GridMediumDistribution::sample(float u0, float u1, float u2, MediumPoint *medPt, float *volumePDF) const {
        SLRAssert_NotImplemented();
    }
    
    float GridMediumDistribution::evaluateVolumePDF(const MediumPoint &medPt) const {
        return 1.0f / volume();
    }
    
    
    
    bool SubGridMediumDistribution::interact(const Ray &ray, float distanceLimit, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                                 MediumInteraction *mi, SampledSpectrum *medThroughput, bool* singleWavelength) const {
        SLRAssert_NotImplemented();
        return false;
    }
    
    SampledSpectrum SubGridMediumDistribution::evaluateTransmittance(Ray &ray, float distanceLimit, const WavelengthSamples &wls, SLR::LightPathSampler &pathSampler, bool *singleWavelength) const {
        SLRAssert_NotImplemented();
        return SampledSpectrum::Zero;
    }
    
    SampledSpectrum SubGridMediumDistribution::extinctionCoefficient(const Point3D &p, const WavelengthSamples &wls) const {
        SLRAssert_NotImplemented();
        return SampledSpectrum::Zero;
    }
    
    void SubGridMediumDistribution::getMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const {
        SLRAssert_NotImplemented();
    }
    
    void SubGridMediumDistribution::sample(float u0, float u1, float u2, MediumPoint *medPt, float *volumePDF) const {
        SLRAssert_NotImplemented();
    }
    
    float SubGridMediumDistribution::evaluateVolumePDF(const MediumPoint &medPt) const {
        return 1.0f / volume();
    }
}

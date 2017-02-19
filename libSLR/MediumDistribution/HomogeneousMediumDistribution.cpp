//
//  HomogeneousMediumDistribution.cpp
//
//  Created by 渡部 心 on 2016/12/31.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "HomogeneousMediumDistribution.h"
#include "light_path_sampler.h"

namespace SLR {
    bool HomogeneousMediumDistribution::interact(const Ray &ray, float distanceLimit, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                                                 MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const {
        FreePathSampler &sampler = pathSampler.getFreePathSampler();
        Point3D queryPoint = ray.org + ray.distMin * ray.dir;
        FloatSum sampledDistance = ray.distMin;
        
        SampledSpectrum extCoeff = extinctionCoefficient(queryPoint, wls);
        float distance = -std::log(sampler.getSample()) / extCoeff[wls.selectedLambda];
        SampledSpectrum transmittance = exp(-extCoeff * std::min(distance, distanceLimit - ray.distMin));
        *singleWavelength = false;
        sampledDistance += distance;
        if (sampledDistance < distanceLimit) {
            Point3D p = ray.org + sampledDistance * ray.dir;
            Point3D param;
            m_region.localCoordinates(p, &param);
            *mi = MediumInteraction(ray.time, sampledDistance, p,
                                    m_sigma_s->evaluate(wls), extCoeff, normalize(ray.dir), param.x, param.y, param.z);
            float distPDF = extCoeff[wls.selectedLambda] * transmittance[wls.selectedLambda];
            *medThroughput = transmittance / distPDF;
            return true;
        }
        
        *medThroughput = transmittance / transmittance[wls.selectedLambda];
        return false;
    }
    
    SampledSpectrum HomogeneousMediumDistribution::evaluateTransmittance(Ray &ray, float distanceLimit, const WavelengthSamples &wls, SLR::LightPathSampler &pathSampler, bool *singleWavelength) const {
        Point3D queryPoint = ray.org + ray.distMin * ray.dir;
        
        SampledSpectrum extCoeff = extinctionCoefficient(queryPoint, wls);
        float distance = distanceLimit - ray.distMin;
        SampledSpectrum transmittance = exp(-extCoeff * distance);
        *singleWavelength = false;
        
        return transmittance;
    }
    
    void HomogeneousMediumDistribution::getMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const {
        ReferenceFrame shadingFrame;
        shadingFrame.z = mi.getIncomingDirection();
        shadingFrame.z.makeCoordinateSystem(&shadingFrame.x, &shadingFrame.y);
        *medPt = MediumPoint(mi, false, shadingFrame);
    }
    
    void HomogeneousMediumDistribution::sample(float u0, float u1, float u2, MediumPoint* medPt, float* volumePDF) const {
        SLRAssert_NotImplemented();
    }
    
    float HomogeneousMediumDistribution::evaluateVolumePDF(const MediumPoint &medPt) const {
        return 1.0f / volume();
    }
}

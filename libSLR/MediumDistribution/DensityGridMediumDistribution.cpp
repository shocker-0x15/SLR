//
//  DensityGridMediumDistribution.cpp
//
//  Created by 渡部 心 on 2017/02/14.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "DensityGridMediumDistribution.h"

#include "../Core/light_path_sampler.h"

namespace SLR {
    float DensityGridMediumDistribution::calcDensity(const Point3D &param) const {
        if (param.x < 0 || param.y < 0 || param.z < 0 ||
            param.x >= 1 || param.y >= 1 || param.z >= 1)
            return 0.0f;
        
        uint32_t lx = std::min((uint32_t)(param.x * (m_numX - 1)), m_numX - 1);
        uint32_t ux = std::min(lx + 1, m_numX - 1);
        uint32_t ly = std::min((uint32_t)(param.y * (m_numY - 1)), m_numY - 1);
        uint32_t uy = std::min(ly + 1, m_numY - 1);
        uint32_t lz = std::min((uint32_t)(param.z * (m_numZ - 1)), m_numZ - 1);
        uint32_t uz = std::min(lz + 1, m_numZ - 1);
        float wlx = 1 - param.x;
        float wux = param.x;
        float wly = 1 - param.y;
        float wuy = param.y;
        float wlz = 1 - param.z;
        float wuz = param.z;
        float density = (m_density_grid[lz][m_numX * ly + lx] * wlz * wly * wlx +
                         m_density_grid[lz][m_numX * ly + ux] * wlz * wly * wux +
                         m_density_grid[lz][m_numX * uy + lx] * wlz * wuy * wlx +
                         m_density_grid[lz][m_numX * uy + ux] * wlz * wuy * wux +
                         m_density_grid[uz][m_numX * ly + lx] * wuz * wly * wlx +
                         m_density_grid[uz][m_numX * ly + ux] * wuz * wly * wux +
                         m_density_grid[uz][m_numX * uy + lx] * wuz * wuy * wlx +
                         m_density_grid[uz][m_numX * uy + ux] * wuz * wuy * wux);
        return density;
    }
    
    bool DensityGridMediumDistribution::subdivide(Allocator* mem, MediumDistribution** fragments, uint32_t* numFragments) const {
        SLRAssert_NotImplemented();
        return true;
    }
    
    bool DensityGridMediumDistribution::interact(const Ray &ray, const RaySegment &segment, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                                                 MediumInteraction *mi, SampledSpectrum *medThroughput, bool* singleWavelength) const {
        SLRAssert(std::isfinite(segment.distMax), "distanceLimit must be a finite value.");
        FreePathSampler &sampler = pathSampler.getFreePathSampler();
        
        SampledSpectrum base_sigma_e = m_base_sigma_e->evaluate(wls);
        
        // delta tracking to sample free path.
        float majorantSelected = majorantExtinctionCoefficientAtWavelength(wls.selectedWavelength());
        *singleWavelength = false;
        bool hit = false;
        float extCoeffSelected = 0.0f;
        float hitDistance = segment.distMax;
        FloatSum sampledDistance = segment.distMin;
        sampledDistance += -std::log(sampler.getSample()) / majorantSelected;
        while (sampledDistance < segment.distMax) {
            Point3D queryPoint = ray.org + sampledDistance * ray.dir;
            Point3D param;
            m_region.calculateLocalCoordinates(queryPoint, &param);
            float density = calcDensity(param);
            SampledSpectrum extCoeff = base_sigma_e * density;
            float probRealCollision = extCoeff[wls.selectedLambdaIndex] / majorantSelected;
            if (sampler.getSample() < probRealCollision) {
                *mi = MediumInteraction(ray.time, sampledDistance, queryPoint, normalize(ray.dir), param.x, param.y, param.z);
                hit = true;
                extCoeffSelected = extCoeff[wls.selectedLambdaIndex];
                hitDistance = sampledDistance;
                break;
            }
            sampledDistance += -std::log(sampler.getSample()) / majorantSelected;
            
            // TODO: handle out of boundary.
        }
        
        // estimate Monte Carlo throughput T(s, wl_j)/p(s, wl_i) by ratio tracking.
        if (wls.wavelengthSelected()) {
            *medThroughput = SampledSpectrum::Zero;
            (*medThroughput)[wls.selectedLambdaIndex] = 1.0f;
        }
        else {
            SampledSpectrum trDiff = SampledSpectrum::One;
            for (int wl = 0; wl < WavelengthSamples::NumComponents; ++wl) {
                if (wl == wls.selectedLambdaIndex)
                    continue;
                float majorantWL = majorantExtinctionCoefficientAtWavelength(wls[wl]);
                sampledDistance = segment.distMin;
                sampledDistance += -std::log(sampler.getSample()) / majorantWL;
                while (sampledDistance < hitDistance) {
                    Point3D queryPoint = ray.org + sampledDistance * ray.dir;
                    Point3D param;
                    m_region.calculateLocalCoordinates(queryPoint, &param);
                    float density = calcDensity(param);
                    SampledSpectrum extCoeff = base_sigma_e * density;
                    float probRealCollision = (extCoeff[wl] - extCoeff[wls.selectedLambdaIndex]) / majorantWL;
                    trDiff[wl] *= (1.0f - probRealCollision);
                    sampledDistance += -std::log(sampler.getSample()) / majorantWL;
                }
            }
            *medThroughput = trDiff;
        }
        if (hit)
            *medThroughput /= extCoeffSelected;
        
        return hit;
    }
    
    SampledSpectrum DensityGridMediumDistribution::evaluateTransmittance(const Ray &ray, const RaySegment &segment, const WavelengthSamples &wls, SLR::LightPathSampler &pathSampler, 
                                                                         bool* singleWavelength) const {
        SLRAssert(std::isfinite(segment.distMax), "distanceLimit must be a finite value.");
        FreePathSampler sampler = pathSampler.getFreePathSampler();
        
        SampledSpectrum base_sigma_e = m_base_sigma_e->evaluate(wls);
        
        *singleWavelength = false;
        
        // estimate transmittance by ratio tracking.
        SampledSpectrum transmittance = SampledSpectrum::One;
        if (wls.wavelengthSelected()) {
            transmittance = SampledSpectrum::Zero;
            transmittance[wls.selectedLambdaIndex] = 1.0f;
            
            float majorantSelected = majorantExtinctionCoefficientAtWavelength(wls.selectedWavelength());
            FloatSum sampledDistance = segment.distMin;
            sampledDistance += -std::log(sampler.getSample()) / majorantSelected;
            while (sampledDistance < segment.distMax) {
                Point3D queryPoint = ray.org + sampledDistance * ray.dir;
                Point3D param;
                m_region.calculateLocalCoordinates(queryPoint, &param);
                float density = calcDensity(param);
                SampledSpectrum extCoeff = base_sigma_e * density;
                float probRealCollision = extCoeff[wls.selectedLambdaIndex] / majorantSelected;
                transmittance[wls.selectedLambdaIndex] *= (1.0f - probRealCollision);
                sampledDistance += -std::log(sampler.getSample()) / majorantSelected;
            }
        }
        else {
            for (int wl = 0; wl < WavelengthSamples::NumComponents; ++wl) {
                float majorantWL = majorantExtinctionCoefficientAtWavelength(wls[wl]);
                FloatSum sampledDistance = segment.distMin;
                sampledDistance += -std::log(sampler.getSample()) / majorantWL;
                while (sampledDistance < segment.distMax) {
                    Point3D queryPoint = ray.org + sampledDistance * ray.dir;
                    Point3D param;
                    m_region.calculateLocalCoordinates(queryPoint, &param);
                    float density = calcDensity(param);
                    SampledSpectrum extCoeff = base_sigma_e * density;
                    float probRealCollision = extCoeff[wl] / majorantWL;
                    transmittance[wl] *= (1.0f - probRealCollision);
                    sampledDistance += -std::log(sampler.getSample()) / majorantWL;
                }
            }
        }
        
        return transmittance;
    }
    
    void DensityGridMediumDistribution::calculateMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const {
        ReferenceFrame shadingFrame(mi.getIncomingDirection());
        *medPt = MediumPoint(mi, false, shadingFrame);
    }
    
    SampledSpectrum DensityGridMediumDistribution::evaluateExtinctionCoefficient(const Point3D &param, const WavelengthSamples &wls) const {
        float density = calcDensity(param);
        return density * m_base_sigma_e->evaluate(wls);
    }
    
    SampledSpectrum DensityGridMediumDistribution::evaluateAlbedo(const Point3D &param, const WavelengthSamples &wls) const {
        if (param.x < 0 || param.y < 0 || param.z < 0 ||
            param.x >= 1 || param.y >= 1 || param.z >= 1)
            return SampledSpectrum::Zero;
        
        return m_base_sigma_s->evaluate(wls).safeDivide(m_base_sigma_e->evaluate(wls));
    }
    
    void DensityGridMediumDistribution::sample(float u0, float u1, float u2, MediumPoint *medPt, float *volumePDF) const {
        SLRAssert_NotImplemented();
    }
    
    float DensityGridMediumDistribution::evaluateVolumePDF(const MediumPoint &medPt) const {
        return 1.0f / volume();
    }
}

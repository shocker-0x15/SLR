//
//  CloudMediumDistribution.cpp
//
//  Created by 渡部 心 on 2017/03/24.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "CloudMediumDistribution.h"

#include "../Core/light_path_sampler.h"
#include "../Helper/bmp_exporter.h"

namespace SLR {
    template <typename RealType>
    class LayeredWorleyNoiseGeneratorTemplate {
        WorleyNoise3DGeneratorTemplate<RealType> m_primaryNoiseGen;
        uint32_t m_numOctaves;
        RealType m_initialFrequency;
        RealType m_initialVariation;
        RealType m_clipValue;
        RealType m_persistence;
        RealType m_frequencyMultiplier;
        
    public:
        LayeredWorleyNoiseGeneratorTemplate(uint32_t numOctaves, RealType initialFrequency, RealType initialVariation, RealType clipValue, 
                                            RealType frequencyMultiplier, RealType persistence) : 
        m_numOctaves(numOctaves), m_initialFrequency(initialFrequency), m_initialVariation(initialVariation), m_clipValue(clipValue),
        m_persistence(persistence), m_frequencyMultiplier(frequencyMultiplier) { } 
        
        RealType evaluate(const Point3DTemplate<RealType> &p) const {
            RealType total = 1;
            RealType frequency = m_initialFrequency;
            RealType variation = m_initialVariation;
            for (int i = 0; i < m_numOctaves; ++i) {
                RealType closestSqDistance;
                uint32_t hashOfClosest;
                uint32_t closestFPIdx;
                m_primaryNoiseGen.evaluate(frequency * p, &closestSqDistance, &hashOfClosest, &closestFPIdx);
                total *= std::fmax(1 - closestSqDistance / 3.0 * variation, 0.0);
                
                variation *= m_persistence;
                frequency *= m_frequencyMultiplier;
            }
            
            return std::fmin(total, m_clipValue);
        }
    };
    
    
    
    template <typename RealType>
    class PerlinWorleyNoiseGeneratorTemplate {
        MultiOctavePerlinNoise3DGeneratorTemplate<RealType> m_perlinGen;
        LayeredWorleyNoiseGeneratorTemplate<RealType> m_worleyGen;
        
    public:
        PerlinWorleyNoiseGeneratorTemplate(uint32_t numOctaves, RealType initialFrequency, RealType initialAmplitude, RealType initialVariation, 
                                           RealType frequencyMultiplier, RealType persistence) :
        m_perlinGen(numOctaves, initialFrequency, initialAmplitude, true, frequencyMultiplier, persistence, -1), 
        m_worleyGen(numOctaves, initialFrequency, initialVariation, 1.0f, frequencyMultiplier, persistence) {
        }
        
        RealType evaluate(const Point3DTemplate<RealType> &p) const {
            return 0.75f * m_perlinGen.evaluate(p) + 0.25f * m_worleyGen.evaluate(p);
        }
    };
    
    
    
    // For debug visualization.
    void CloudMediumDistribution::exportBMP(const std::string &filename, float gamma) const {
//        LayeredWorleyNoiseGeneratorTemplate<float> layeredWorleyGen(3, 10.0f, 8.0f, 1.0f, 2.0f, 0.5f);
        PerlinWorleyNoiseGeneratorTemplate<float> perlinWorleyGen(3, 10.0f, 1.0f, 8.0f, 2.0f, 0.5f);
        
        const uint32_t width = 512;
        const uint32_t height = 512;
        uint32_t byteWidth = width * 3 + width % 4;
        uint8_t* data = (uint8_t*)malloc(height * byteWidth);
        
        printf("start to generate a test noise image. ...");
        for (int i = 0; i < height; ++i) {
            for (int j = 0; j < width; ++j) {
                Point3D p(Point3D((j + 0.5f) / width, (i + 0.5f) / height, 0.0f));
                float value = perlinWorleyGen.evaluate(p);
//                float value = layeredWorleyGen.evaluate(p);
                uint8_t pixVal = uint8_t(std::pow(std::fmin(value, 1.0f), 1.0f / gamma) * 255);
                
                uint32_t idx = (height - i - 1) * byteWidth + 3 * j;
                data[idx + 0] = pixVal;
                data[idx + 1] = pixVal;
                data[idx + 2] = pixVal;
            }
        }
        saveBMP(filename.c_str(), data, width, height);
        free(data);
        printf("done\n");
    };
    
    void CloudMediumDistribution::saveToFile(const char* fileName, uint32_t resX, uint32_t resY, uint32_t resZ) const {
        printf("write to %s ...", fileName);
        fflush(stdout);
        FILE* fp = fopen(fileName, "wb");
        
        fwrite(&resX, sizeof(resX), 1, fp);
        fwrite(&resY, sizeof(resY), 1, fp);
        fwrite(&resZ, sizeof(resZ), 1, fp);
        
        float* zSlice = new float[resX * resY];
        for (int iz = 0; iz < resZ; ++iz) {
            float pz = (float)iz / (resZ - 1);
            for (int iy = 0; iy < resY; ++iy) {
                float py = (float)iy / (resY - 1);
                for (int ix = 0; ix < resX; ++ix) {
                    float px = (float)ix / (resX - 1);
                    zSlice[resX * iy + ix] = calcDensity(Point3D(px, py, pz));
                }
            }
            fwrite(zSlice, sizeof(float), resX * resY, fp);
        }
        delete[] zSlice;
        
        fclose(fp);
        printf("done.\n");
    }
    
    float CloudMediumDistribution::calcDensity(const Point3D &param) const {
        if (param.x < 0 || param.y < 0 || param.z < 0 ||
            param.x > 1 || param.y > 1 || param.z > 1)
            return 0.0f;
        
        Point3D position = m_aspect * param;
        
        float distributionValue = m_distributionGenerator.evaluate(position);
        float distributionThreshold = 0.525f;
        if (distributionValue > distributionThreshold)
            return 1.0;//m_densityGenerator.evaluate(position.x, position.y, position.z);
        return 0.0f;
    }
    
    bool CloudMediumDistribution::subdivide(Allocator* mem, MediumDistribution** fragments, uint32_t* numFragments) const {
        SLRAssert_NotImplemented();
        return true;
    }
    
    bool CloudMediumDistribution::interact(const Ray &ray, const RaySegment &segment, const WavelengthSamples &wls, LightPathSampler &pathSampler,
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
        *medThroughput = SampledSpectrum::One;
        if (hit)
            *medThroughput /= extCoeffSelected;
        
        return hit;
    }
    
    SampledSpectrum CloudMediumDistribution::evaluateTransmittance(const Ray &ray, const RaySegment &segment, const WavelengthSamples &wls, SLR::LightPathSampler &pathSampler, 
                                                                   bool* singleWavelength) const {
        SLRAssert(std::isfinite(segment.distMax), "distanceLimit must be a finite value.");
        FreePathSampler sampler = pathSampler.getFreePathSampler();
        
        SampledSpectrum base_sigma_e = m_base_sigma_e->evaluate(wls);
        
        // estimate transmittance by ratio tracking.
        *singleWavelength = false;
        SampledSpectrum transmittance = SampledSpectrum::One;
        
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
            if (probRealCollision >= 1.0f) {
                transmittance = SampledSpectrum::Zero;
                break;
            }
            transmittance *= (1.0f - probRealCollision);
            sampledDistance += -std::log(sampler.getSample()) / majorantSelected;
        }
        SLRAssert(transmittance.allFinite() && !transmittance.hasMinus(), "Invalid transmittance value.");
        
        return transmittance;
    }
    
    void CloudMediumDistribution::calculateMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const {
        ReferenceFrame shadingFrame;
        shadingFrame.z = mi.getIncomingDirection();
        shadingFrame.z.makeCoordinateSystem(&shadingFrame.x, &shadingFrame.y);
        *medPt = MediumPoint(mi, false, shadingFrame);
    }
    
    SampledSpectrum CloudMediumDistribution::evaluateExtinctionCoefficient(const Point3D &param, const WavelengthSamples &wls) const {
        float density = calcDensity(param);
        return density * m_base_sigma_e->evaluate(wls);
    }
    
    SampledSpectrum CloudMediumDistribution::evaluateAlbedo(const Point3D &param, const WavelengthSamples &wls) const {
        if (param.x < 0 || param.y < 0 || param.z < 0 ||
            param.x > 1 || param.y > 1 || param.z > 1)
            return SampledSpectrum::Zero;
        
        return m_albedo->evaluate(wls);
    }
    
    void CloudMediumDistribution::sample(float u0, float u1, float u2, MediumPoint *medPt, float *volumePDF) const {
        SLRAssert_NotImplemented();
    }
    
    float CloudMediumDistribution::evaluateVolumePDF(const MediumPoint &medPt) const {
        return 1.0f / volume();
    }
}

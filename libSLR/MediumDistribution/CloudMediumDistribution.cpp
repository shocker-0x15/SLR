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
    static float enhance(float x, float threshold, float power) {
        if (x < threshold)
            return threshold * std::pow(x / threshold, power);
        else
            return 1 - (1 - threshold) * std::pow((1 - x) / (1 - threshold), power);
    }
    
    static float enhanceLower(float x, float threshold, float power) {
        if (x < threshold)
            return threshold * std::pow(x / threshold, power);
        else
            return x;
    }
    
    
    
    template <typename RealType>
    RealType LayeredWorleyNoiseGeneratorTemplate<RealType>::evaluate(const Point3DTemplate<RealType> &p) const {
//        RealType total = 1;
//        RealType frequency = m_initialFrequency;
//        RealType variation = m_initialVariation;
//        for (int i = 0; i < m_numOctaves; ++i) {
//            RealType closestSqDistance;
//            uint32_t hashOfClosest;
//            uint32_t closestFPIdx;
//            m_primaryNoiseGen.evaluate(p, frequency, &closestSqDistance, &hashOfClosest, &closestFPIdx);
//            total *= std::fmax(1 - closestSqDistance / 3.0 * variation, 0.0);
//            
//            variation *= m_persistence;
//            frequency *= m_frequencyMultiplier;
//        }
//        
//        return total;
        
        RealType total = 0;
        RealType frequency = m_initialFrequency;
        RealType variation = m_initialVariation;
        RealType sumAmp = 0.0f;
        for (int i = 0; i < m_numOctaves; ++i) {
            RealType closestSqDistance;
            uint32_t hashOfClosest;
            uint32_t closestFPIdx;
            m_primaryNoiseGen.evaluate(p, frequency, &closestSqDistance, &hashOfClosest, &closestFPIdx);
            total += (1 - closestSqDistance) * variation;
            
            sumAmp += variation;
            variation *= m_persistence;
            frequency *= m_frequencyMultiplier;
        }
        
        return std::max<RealType>(total / sumAmp, 0.0f);
    }
    
    template class SLR_API LayeredWorleyNoiseGeneratorTemplate<float>;
    template class SLR_API LayeredWorleyNoiseGeneratorTemplate<double>;
    
    
    
    template class SLR_API PerlinWorleyNoiseGeneratorTemplate<float>;
    template class SLR_API PerlinWorleyNoiseGeneratorTemplate<double>;
    
    
    
    const float CloudMediumDistribution::CloudBaseAltitude = 1500;
    const float CloudMediumDistribution::CloudTopAltitude = 5000;
    
    // For debug visualization.
    void CloudMediumDistribution::exportBMPs() const {
//        LayeredWorleyNoiseGeneratorTemplate<float> layeredWorleyGen(3, 10.0f, 8.0f, 1.0f, 2.0f, 0.5f);
//        PerlinWorleyNoiseGeneratorTemplate<float> perlinWorleyGen(3, 10.0f, 1.0f, 8.0f, 2.0f, 0.5f);
//        CurlNoise3DGeneratorTemplate<float> curlGen(3, 10.0f);
        
        const auto produceSliceImage = [this](const std::string &filename, std::function<float(const Point3D &)> func, uint32_t resX, uint32_t resZ, float gamma) {
            uint32_t byteWidth = resX * 3 + resX % 4;
            uint8_t* data = (uint8_t*)malloc(resZ * byteWidth);
            
            printf("start to generate a test image: %s. ...", filename.c_str());
            for (int i = 0; i < resZ; ++i) {
                for (int j = 0; j < resX; ++j) {
                    Point3D p((j + 0.5f) / resX, 0.0f, (i + 0.5f) / resZ);
                    
                    float value = func(p);
                    
                    uint8_t pixVal = uint8_t(std::pow(std::clamp(value, 0.0f, 1.0f), 1.0f / gamma) * 255);
                    
                    uint32_t idx = (resZ - i - 1) * byteWidth + 3 * j;
                    data[idx + 0] = pixVal;
                    data[idx + 1] = pixVal;
                    data[idx + 2] = pixVal;
                }
            }
            saveBMP(filename.c_str(), data, resX, resZ);
            free(data);
            printf("done\n");
        };
        
        produceSliceImage("cloud_coverage.bmp", std::bind(&CloudMediumDistribution::calcCoverage, this, std::placeholders::_1), WeatherNumX, WeatherNumZ, 1.0f);
        produceSliceImage("cloud_cloud_type.bmp", std::bind(&CloudMediumDistribution::calcCloudType, this, std::placeholders::_1), WeatherNumX, WeatherNumZ, 1.0f);
        produceSliceImage("cloud_base_shape.bmp", std::bind(&CloudMediumDistribution::calcBaseShape, this, std::placeholders::_1), LowResNumX, LowResNumZ, 1.0f);
        produceSliceImage("cloud_erosion.bmp", std::bind(&CloudMediumDistribution::calcErosion, this, std::placeholders::_1), HighResNumX, HighResNumZ, 1.0f);
//        exit(0);
    };
    
    void CloudMediumDistribution::saveToFile(const std::string &name) const {
        const auto export2DData = [this](const std::string &fileName, std::function<float(const Point3D &p)> func, uint32_t resX, uint32_t resZ) {
            printf("write to %s\n", fileName.c_str());
            fflush(stdout);
            FILE* fp = fopen(fileName.c_str(), "wb");
            
            uint32_t numDims = 2;
            fwrite(&numDims, sizeof(numDims), 1, fp);
            fwrite(&resX, sizeof(resX), 1, fp);
            fwrite(&resZ, sizeof(resZ), 1, fp);
            
            float* ySlice = new float[resX * resZ];
            for (int iz = 0; iz < resZ; ++iz) {
                float pz = (float)iz / (resZ - 1);
                for (int ix = 0; ix < resX; ++ix) {
                    float px = (float)ix / (resX - 1);
                    Point3D p(px, 0, pz);
                    ySlice[resX * iz + ix] = func(p);
                }
            }
            fwrite(ySlice, sizeof(float), resX * resZ, fp);
            delete[] ySlice;
            
            fclose(fp);
            printf("done.\n");
        };
        
        const auto export3DData = [this](const std::string &fileName, std::function<float(const Point3D &p)> func, uint32_t resX, uint32_t resY, uint32_t resZ) {
            printf("write to %s\n", fileName.c_str());
            fflush(stdout);
            FILE* fp = fopen(fileName.c_str(), "wb");
            
            uint32_t numDims = 3;
            fwrite(&numDims, sizeof(numDims), 1, fp);
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
                        Point3D p(px, py, pz);
                        zSlice[resX * iy + ix] = func(p);
                    }
                }
                fwrite(zSlice, sizeof(float), resX * resY, fp);
                
                printf("%3u/%3u\n", iz + 1, resZ);
            }
            delete[] zSlice;
            
            fclose(fp);
            printf("done.\n");
        };
        
        exportBMPs();
        
        export2DData(name + "_coverage", std::bind(&CloudMediumDistribution::calcCoverage, this, std::placeholders::_1), WeatherNumX, WeatherNumZ);
        export2DData(name + "_cloud_type", std::bind(&CloudMediumDistribution::calcCloudType, this, std::placeholders::_1), WeatherNumX, WeatherNumZ);
        export3DData(name + "_base_shape", std::bind(&CloudMediumDistribution::calcBaseShape, this, std::placeholders::_1), LowResNumX, LowResNumY, LowResNumZ);
        export3DData(name + "_erosion", std::bind(&CloudMediumDistribution::calcErosion, this, std::placeholders::_1), HighResNumX, HighResNumY, HighResNumZ);
    }
    
    float CloudMediumDistribution::calcCoverage(const Point3D &position) const {
        float sum = 0.0f;
        const uint32_t NumOverlays = 3;
        for (int i = 0; i < NumOverlays; ++i) {
            float value = m_coverageGenerator.evaluate(Point3D(position.x, 0.0f, position.z) + i * Vector3D(0.5f, 0.3f, 0.2f));
            const float threshold = 0.5f;
            value = smoothstep(threshold - 0.1f, threshold + 0.1f, value);
            
            sum += value;
        }
        float value = m_coverageGenerator.evaluate(Point3D(position.x, 0.0f, position.z) + Vector3D(-0.27f, 0.15f, -0.35f));
        const float threshold = 0.5f;
        value = smoothstep(threshold - 0.05f, threshold + 0.05f, value);
        sum *= value;
        
        return sum / NumOverlays;
    }
    
    float CloudMediumDistribution::calcCloudType(const Point3D &position) const {
        float value = m_coverageGenerator.evaluate(Point3D(position.x, 0.0f, position.z) + Vector3D(-0.5f, 0.7f, 0.4f));
        const float threshold = 0.4f;
        
        return saturate(remap(value, threshold, 1.0f, 0.0f, 0.55f) / 0.25f);
    }
    
    float CloudMediumDistribution::calcBaseShape(const Point3D &position) const {
//        // 3D checker board pattern for Debug
//        int32_t px, py, pz;
//        px = std::min<int32_t>(position.x * LowResNumX, LowResNumX - 1);
//        py = std::min<int32_t>(position.y * LowResNumY, LowResNumY - 1);
//        pz = std::min<int32_t>(position.z * LowResNumZ, LowResNumZ - 1);
//        return (px + py + pz) % 2 == 0 ? 1.0f : 0.0f;
        
        float baseShapeValue = m_baseShapeGenerator.evaluate(position);
        float floorValue = m_baseShapeExtraGenerator.evaluate(position + Vector3D(3, -5, 2));
        baseShapeValue = remap(baseShapeValue, -floorValue, 1.0f, 0.0f, 1.0f);
//        baseShapeValue = remap(baseShapeValue, -(1 - floorValue), 1.0f, 0.0f, 1.0f);
//        baseShapeValue = remap(baseShapeValue, (1 - floorValue), 1.0f, 0.0f, 1.0f);
        
        return baseShapeValue;
    }
    
    float CloudMediumDistribution::calcErosion(const Point3D &position) const {
        float erosion = m_erosionGenerator.evaluate(position + Vector3D(5.0f, 3.0f, 2.0f));
        
        return erosion;
    }
    
    float CloudMediumDistribution::calcHeightGradient(float cloudType, float h) const {
        // [0, 1] => [1500, 2750]
        const float stratusGrad[] = {1525, 1562.5, 1612.5, 1637.5}; // 0.02f, 0.05f, 0.09f, 0.11f
        const float stratocumulusGrad[] = {1525, 1750, 2100, 2281.25}; // 0.02f, 0.2f, 0.48f, 0.625f
        const float cumulusGrad[] = {1512.5, 1578.125, 2475, 2750}; // 0.01f, 0.0625f, 0.78f, 1.0f
        const float cumulonimbusGrad[] = {1512.5, 1578.125, 5000, 5500};
        
        float pStratus = 1 - saturate(cloudType * 3);//1 - saturate(std::abs(cloudType - 0.0f / 3.0f) * 3);
        float pStratocumulus = 1 - saturate(std::abs(cloudType - 1.0f / 3.0f) * 3);
        float pCumulus = 1 - saturate(std::abs(cloudType - 2.0f / 3.0f) * 3);
        float pCumulonimbus = saturate((cloudType - 2.0f / 3.0f) * 3);//1 - saturate(std::abs(cloudType - 3.0f / 3.0f) * 3);
        
        float grad[4];
        for (int i = 0; i < 4; ++i)
            grad[i] = (pStratus * stratusGrad[i] + 
                       pStratocumulus * stratocumulusGrad[i] + 
                       pCumulus * cumulusGrad[i] + 
                       pCumulonimbus * cumulonimbusGrad[i]);
        
        //           ________________
        //          /                \
        // ________/                  \________
        return smoothstep(grad[0], grad[1], h) - smoothstep(grad[2], grad[3], h);
    }
    
    float CloudMediumDistribution::calcDensity(const Point3D &param) const {
        if (param.x < 0 || param.y < 0 || param.z < 0 ||
            param.x > 1 || param.y > 1 || param.z > 1)
            return 0.0f;
        
        Point3D position = (m_region.minP + (m_region.maxP - m_region.minP) * param) / m_featureScale;
        
#if defined(CLOUD_GENERATION)
        SLRAssert_NotImplemented();
        return 1.0f;
#else
        const float minimumDensity = 1.0f / (8e+4f * 0.1f);
        
        // JP: x:100km x z:100km を基準とする。
        const float DefaultWeatherScale = 1e-5f;
        const float DefaultBaseShapeScale = 5.314f * 1e-5f * 10;
        const float DefaultErosionScale = 5.739f * 1e-4f * 10;
        
        float baseShape = m_baseShape.evaluate(DefaultBaseShapeScale * position.x, 
                                               DefaultBaseShapeScale * position.y, 
                                               DefaultBaseShapeScale * position.z);
        if (baseShape <= 0.0f)
            return minimumDensity;
        
        float cloudType = m_cloudType.evaluate(DefaultWeatherScale * position.x, 
                                               DefaultWeatherScale * position.z);
        float heightGradient = calcHeightGradient(cloudType * 0.666f, position.y);
        baseShape *= heightGradient;
        if (baseShape <= 0.0f)
            return minimumDensity;
        
        float coverage = m_coverage.evaluate(DefaultWeatherScale * position.x, 
                                             DefaultWeatherScale * position.z);
        float ret = saturate(remap(baseShape, 1 - coverage, 1.0f, 0.0f, 1.0f)) * coverage;
        if (ret <= 0.0f)
            return minimumDensity;
        
        float erosion = m_erosion.evaluate(DefaultErosionScale * position.x, 
                                           DefaultErosionScale * position.y, 
                                           DefaultErosionScale * position.z);
        float heightFrac = saturate((position.y - CloudBaseAltitude) / (CloudTopAltitude - CloudBaseAltitude));
        erosion = erosion * (1 - heightFrac) + (1 - erosion) * heightFrac;
        ret = saturate(remap(ret, erosion * 0.2f, 1.0f, 0.0f, 1.0f));
        
        return std::max(enhanceLower(ret, 0.5f, 0.2f), minimumDensity);
        
//        const auto produceSliceImage = [this](const std::string &filename, std::function<float(float, float, float)> func, uint32_t iy, uint32_t resX, uint32_t resY, uint32_t resZ, float gamma) {
//            uint32_t byteWidth = resX * 3 + resX % 4;
//            uint8_t* data = (uint8_t*)malloc(resZ * byteWidth);
//            
//            printf("start to generate a test image: %s. ...", filename.c_str());
//            for (int i = 0; i < resZ; ++i) {
//                for (int j = 0; j < resX; ++j) {
//                    Point3D p((j + 0.5f) / resX + resX, (iy + 0.5f) / resY + resY, (i + 0.5f) / resZ + resZ);
//                    
//                    float value = func(p.x, p.y, p.z);
//                    
//                    uint8_t pixVal = uint8_t(std::pow(std::clamp(value, 0.0f, 1.0f), 1.0f / gamma) * 255);
//                    
//                    uint32_t idx = (resZ - i - 1) * byteWidth + 3 * j;
//                    data[idx + 0] = pixVal;
//                    data[idx + 1] = pixVal;
//                    data[idx + 2] = pixVal;
//                }
//            }
//            saveBMP(filename.c_str(), data, resX, resZ);
//            free(data);
//            printf("done\n");  
//        };
//        for (int i = 0; i < LowResNumY; ++i) {
//            char filename[256];
//            sprintf(filename, "cloud_base_shape_%03u.bmp", i);
//            produceSliceImage(filename, std::bind(&CloudMediumDistribution::Float3DGrid::evaluate, &m_baseShape, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3), 
//                              i, LowResNumX, LowResNumY, LowResNumZ, 1.0f);
//        }
//        exit(-1);

//        float density = 1.0f;
//        float coverage = m_coverage.evaluate(param.x + 0.3f, param.z + 0.3f);
//        density *= coverage;
//        if (density <= 0.0f)
//            return 0.0f;
//        
//        float cloudType = m_cloudType.evaluate(param.x, param.z);
//        float heightSignal = calcHeightGradient(cloudType * 0.666f, position.y);
//        density *= heightSignal;
//        if (density <= 0.0f)
//            return 0.0f;
//        
//        float baseShape = m_baseShape.evaluate(DefaultBaseShapeScale * position.x, 
//                                               DefaultBaseShapeScale * position.y, 
//                                               DefaultBaseShapeScale * position.z);
////        if (density > 0)
////            printf("(%g, %g, %g): %g\n", 
////                   DefaultBaseShapeScale * position.x, 
////                   DefaultBaseShapeScale * position.y, 
////                   DefaultBaseShapeScale * position.z, 
////                   baseShape);
//        density *= baseShape;
//        if (density <= 0.0f)
//            return 0.0f;
//        
////        float erosion = m_erosion.evaluate(DefaultErosionScale * position.x, 
////                                           DefaultErosionScale * position.y, 
////                                           DefaultErosionScale * position.z);
////        float heightFrac = saturate((position.y - CloudBaseAltitude) / (CloudTopAltitude - CloudBaseAltitude));
////        erosion = erosion * (1 - heightFrac) + (1 - erosion) * heightFrac;
////        density -= 0.2f * erosion;
//        
//        return enhanceLower(saturate(density), 0.5f, 0.2f);
#endif
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
        SLRAssert(transmittance.allFinite() && !transmittance.hasNegative(), "Invalid transmittance value.");
        
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

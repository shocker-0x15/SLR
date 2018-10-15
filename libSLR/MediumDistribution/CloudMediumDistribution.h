//
//  CloudMediumDistribution.h
//
//  Created by 渡部 心 on 2017/03/24.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef __SLR_CloudMediumDistribution__
#define __SLR_CloudMediumDistribution__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/geometry.h"
#include "../Core/distributions.h"

namespace SLR {
    template <typename RealType>
    class LayeredWorleyNoiseGeneratorTemplate {
        WorleyNoise3DGeneratorTemplate<RealType> m_primaryNoiseGen;
        uint32_t m_numOctaves;
        RealType m_initialFrequency;
        RealType m_initialVariation;
        RealType m_persistence;
        RealType m_frequencyMultiplier;
        int32_t m_repeat;
        
    public:
        LayeredWorleyNoiseGeneratorTemplate(uint32_t numOctaves, RealType initialFrequency, RealType initialVariation, 
                                            RealType frequencyMultiplier, RealType persistence, int32_t repeat) : 
        m_primaryNoiseGen(repeat), m_numOctaves(numOctaves), m_initialFrequency(initialFrequency), m_initialVariation(initialVariation),
        m_persistence(persistence), m_frequencyMultiplier(frequencyMultiplier) { } 
        
        RealType evaluate(const Point3DTemplate<RealType> &p) const;
    };
    
    
    
    template <typename RealType>
    class PerlinWorleyNoiseGeneratorTemplate {
        MultiOctavePerlinNoise3DGeneratorTemplate<RealType> m_perlinGen;
        LayeredWorleyNoiseGeneratorTemplate<RealType> m_worleyGen;
        
    public:
        PerlinWorleyNoiseGeneratorTemplate(uint32_t numOctaves, RealType initialFrequency, RealType initialAmplitude, RealType initialVariation, 
                                           RealType frequencyMultiplier, RealType persistence, int32_t repeat) :
        m_perlinGen(numOctaves, initialFrequency, initialAmplitude, true, frequencyMultiplier, persistence, repeat), 
        m_worleyGen(numOctaves, initialFrequency, initialVariation, frequencyMultiplier, persistence, repeat) {
        }
        
        RealType evaluate(const Point3DTemplate<RealType> &p) const {
//            return m_perlinGen.evaluate(p);
//            return m_worleyGen.evaluate(p);
//            return 0.75f * m_perlinGen.evaluate(p) + 0.25f * m_worleyGen.evaluate(p);
//            return saturate<RealType>(remap<RealType>(m_perlinGen.evaluate(p), 0.0f, 1.0, m_worleyGen.evaluate(p), 1.0)); // Frostbite-like
            
            // From Nubis slide:
            // Our approach, In principle, was to subtract the web like shapes of Worley noise from the low density regions of perlin noise in order to introduce round shapes there.
            return remap<RealType>(m_perlinGen.evaluate(p), 1.0f * (1 - m_worleyGen.evaluate(p)), 1.0, 0.0, 1.0); // Nubis-like
        }
        
        RealType getSupValue() const {
//            return 0.75f * m_perlinGen.getSupValue() + 0.25f;
            return m_perlinGen.getSupValue();
        }
    };
    
    
    
//#define CLOUD_GENERATION
    class CloudMediumDistribution : public MediumDistribution {
        std::array<float, NumStrataForStorage> m_majorantExtinctionCoefficient;
        BoundingBox3D m_region;
        float m_featureScale;
        const AssetSpectrum* m_base_sigma_e;
        const AssetSpectrum* m_albedo;
        MultiOctavePerlinNoise3DGeneratorTemplate<float> m_coverageGenerator;
        PerlinWorleyNoiseGeneratorTemplate<float> m_baseShapeGenerator;
        LayeredWorleyNoiseGeneratorTemplate<float> m_baseShapeExtraGenerator;
        LayeredWorleyNoiseGeneratorTemplate<float> m_erosionGenerator;
        
        static const float CloudBaseAltitude;
        static const float CloudTopAltitude;
        static const uint32_t WeatherNumX = 256;
        static const uint32_t WeatherNumZ = 256;
        static const uint32_t LowResNumX = 128;
        static const uint32_t LowResNumY = 128;
        static const uint32_t LowResNumZ = 128;
        static const uint32_t HighResNumX = 32;
        static const uint32_t HighResNumY = 32;
        static const uint32_t HighResNumZ = 32;
#if !defined(CLOUD_GENERATION)
        struct Float2DGrid {
            uint32_t numX, numZ;
            float* data;
            
            Float2DGrid() : data(nullptr) { }
            ~Float2DGrid() {
                if (data)
                    delete[] data;
            }
            
            void open(const std::string &filename) {
                FILE* fp = fopen(filename.c_str(), "rb");
                if (fp == nullptr) {
                    printf("failed to read grid file: %s.", filename.c_str());
                    return;
                }
                
                uint32_t numDims;
                fread(&numDims, sizeof(uint32_t), 1, fp);
                if (numDims != 2) {
                    fclose(fp);
                    data = nullptr;
                    return;
                }
                fread(&numX, sizeof(uint32_t), 1, fp);
                fread(&numZ, sizeof(uint32_t), 1, fp);
                
                data = new float[numX * numZ];
                fread(data, sizeof(float), numX * numZ, fp);
                
                fclose(fp);
            }
            
            float evaluate(float x, float z) const {
                x = std::fmod(x, 1.0f);
                z = std::fmod(z, 1.0f);
                if (x < 0)
                    x += 1.0f;
                if (z < 0)
                    z += 1.0f;
                
                uint32_t lx = std::min<uint32_t>((uint32_t)(x * (numX - 1)), numX - 1);
                uint32_t ux = std::min<uint32_t>(lx + 1, numX - 1);
                uint32_t lz = std::min<uint32_t>((uint32_t)(z * (numZ - 1)), numZ - 1);
                uint32_t uz = std::min<uint32_t>(lz + 1, numZ - 1);
                float wux = x * (numX - 1) - lx;
                float wlx = 1 - wux;
                float wuz = z * (numZ - 1) - lz;
                float wlz = 1 - wuz;
                
                return (wlz * (wlx * data[numX * lz + lx] + 
                               wux * data[numX * lz + ux]) + 
                        wuz * (wlx * data[numX * uz + lx] + 
                               wux * data[numX * uz + ux])); 
            }
        };
        struct Float3DGrid {
            uint32_t numX, numY, numZ;
            float** data;
            
            Float3DGrid() : data(nullptr) { }
            ~Float3DGrid() {
                if (data) {
                    for (int iz = numZ - 1; iz >= 0; --iz)
                        delete[] data[iz];
                    delete[] data;
                }
            }
            
            void open(const std::string &filename) {
                FILE* fp = fopen(filename.c_str(), "rb");
                if (fp == nullptr) {
                    printf("failed to read grid file: %s.", filename.c_str());
                    return;
                }
                
                uint32_t numDims;
                fread(&numDims, sizeof(uint32_t), 1, fp);
                if (numDims != 3) {
                    fclose(fp);
                    data = nullptr;
                    return;
                }
                fread(&numX, sizeof(uint32_t), 1, fp);
                fread(&numY, sizeof(uint32_t), 1, fp);
                fread(&numZ, sizeof(uint32_t), 1, fp);
                
                data = new float*[numZ];
                for (int iz = 0; iz < numZ; ++iz) {
                    float* &ySlice = data[iz];
                    ySlice = new float[numX * numY];
                    fread(ySlice, sizeof(float), numX * numY, fp);
                }
                
                fclose(fp);
            }
            
            float evaluate(float x, float y, float z) const {
                x = std::fmod(x, 1.0f);
                y = std::fmod(y, 1.0f);
                z = std::fmod(z, 1.0f);
                if (x < 0)
                    x += 1.0f;
                if (y < 0)
                    y += 1.0f;
                if (z < 0)
                    z += 1.0f;
                
                uint32_t lx = std::min<uint32_t>((uint32_t)(x * (numX - 1)), numX - 1);
                uint32_t ux = std::min<uint32_t>(lx + 1, numX - 1);
                uint32_t ly = std::min<uint32_t>((uint32_t)(y * (numY - 1)), numY - 1);
                uint32_t uy = std::min<uint32_t>(ly + 1, numY - 1);
                uint32_t lz = std::min<uint32_t>((uint32_t)(z * (numZ - 1)), numZ - 1);
                uint32_t uz = std::min<uint32_t>(lz + 1, numZ - 1);
                float wux = x * (numX - 1) - lx;
                float wlx = 1 - wux;
                float wuy = y * (numY - 1) - ly;
                float wly = 1 - wuy;
                float wuz = z * (numZ - 1) - lz;
                float wlz = 1 - wuz;
                
                return (wlz * (wly * (wlx * data[lz][numX * ly + lx] + 
                                      wux * data[lz][numX * ly + ux]) + 
                               wuy * (wlx * data[lz][numX * uy + lx] + 
                                      wux * data[lz][numX * uy + ux])) + 
                        wuz * (wly * (wlx * data[uz][numX * ly + lx] + 
                                      wux * data[uz][numX * ly + ux]) + 
                               wuy * (wlx * data[uz][numX * uy + lx] + 
                                      wux * data[uz][numX * uy + ux]))); 
            }
        };
        
        float m_maxDensity;
        Float2DGrid m_coverage;
        Float2DGrid m_cloudType;
        Float3DGrid m_baseShape;
        Float3DGrid m_erosion;
#endif
        
        void exportBMPs() const;
        void saveToFile(const std::string &name) const;
        float calcCoverage(const Point3D &param) const;
        float calcCloudType(const Point3D &param) const;
        float calcBaseShape(const Point3D &param) const;
        float calcErosion(const Point3D &param) const;
        float calcHeightGradient(float cloudType, float h) const;
        float calcDensity(const Point3D &param) const;
    public:
        CloudMediumDistribution(const BoundingBox3D &region, float featureScale, float density, uint32_t rngSeed) : 
        m_region(Point3D(region.minP.x, 0 * featureScale, region.minP.z), Point3D(region.maxP.x, 5000 * featureScale, region.maxP.z)),
        m_featureScale(featureScale), 
        m_coverageGenerator(5, 10.0f, 1.0f, true, 2.0f, 0.5f, 1), 
        m_baseShapeGenerator(3, 8.0f, 1.0f, 1.0f, 2.0f, 0.5f, 1),
        m_baseShapeExtraGenerator(5, 8.0f, 1.0f, 2.0f, 0.5f, 1),
        m_erosionGenerator(3, 4.0f, 1.0f, 2.0f, 0.5f, 1) {
            float sigma_e_values[] = {0.1f * density, 0.1f * density};
            m_base_sigma_e = new RegularContinuousSpectrum(WavelengthLowBound, WavelengthHighBound, sigma_e_values, 2);
            float albedo_values[] = {0.999f, 0.999f};
//            float albedo_values[] = {0.0f, 0.0f};
            m_albedo = new RegularContinuousSpectrum(WavelengthLowBound, WavelengthHighBound, albedo_values, 2);

#if defined(CLOUD_GENERATION)
            m_base_sigma_e->calcBounds(NumStrataForStorage, m_majorantExtinctionCoefficient.data());
            for (int i = 0; i < NumStrataForStorage; ++i)
                m_majorantExtinctionCoefficient[i] *= m_baseShapeGenerator.getSupValue();
            
            saveToFile("cloud");
            exit(0);
#else       
            m_coverage.open("cloud_coverage");
            m_cloudType.open("cloud_cloud_type");
//            m_baseShape.openTGA("noiseShapePacked.tga");
//            m_erosion.openTGA("noiseErosionPacked.tga");
            m_baseShape.open("cloud_base_shape");
            m_erosion.open("cloud_erosion");

            m_maxDensity = 1.0f;
//            m_maxDensity = -INFINITY;
//            for (int z = 0; z < m_baseShape.numZ; ++z) {
//                float pz = (z + 0.5f) / m_baseShape.numZ; 
//                for (int y = 0; y < m_baseShape.numY; ++y) {
//                    float py = (y + 0.5f) / m_baseShape.numY;
//                    for (int x = 0; x < m_baseShape.numX; ++x) {
//                        float px = (x + 0.5f) / m_baseShape.numX;
//                        m_maxDensity = calcDensity(Point3D(px, py, pz));
//                    }
//                }
//            }
            m_base_sigma_e->calcBounds(NumStrataForStorage, m_majorantExtinctionCoefficient.data());
            for (int i = 0; i < NumStrataForStorage; ++i)
                m_majorantExtinctionCoefficient[i] *= m_maxDensity;
#endif
        }
        ~CloudMediumDistribution() {
            delete m_albedo;
            delete m_base_sigma_e;
        }
        
        float majorantExtinctionCoefficientAtWavelength(float wl) const override {
            int index = (wl - WavelengthLowBound) / (WavelengthHighBound - WavelengthLowBound) * NumStrataForStorage;
            index = std::clamp(index, 0, (int)NumStrataForStorage - 1);
            return m_majorantExtinctionCoefficient[index];
        }
        
        bool subdivide(Allocator* mem, MediumDistribution** fragments, uint32_t* numFragments) const override;
        
        BoundingBox3D bounds() const override { return m_region; }
        bool contains(const Point3D &p) const override { return m_region.contains(p); }
        bool intersectBoundary(const Ray &ray, const RaySegment &segment, float* distToBoundary, bool* enter) const override {
            return m_region.intersectBoundary(ray, segment, distToBoundary, enter);
        }
        bool interact(const Ray &ray, const RaySegment &segment, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                      MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const override;
        SampledSpectrum evaluateTransmittance(const Ray &ray, const RaySegment &segment, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                                              bool* singleWavelength) const override;
        void calculateMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const override;
        SampledSpectrum evaluateExtinctionCoefficient(const Point3D &param, const WavelengthSamples &wls) const override;
        SampledSpectrum evaluateAlbedo(const Point3D &param, const WavelengthSamples &wls) const override;
        float volume() const override { return m_region.volume(); }
        void sample(float u0, float u1, float u2, MediumPoint* medPt, float* volumePDF) const override;
        float evaluateVolumePDF(const MediumPoint& medPt) const override;
    };
}

#endif /* __SLR_CloudMediumDistribution__ */

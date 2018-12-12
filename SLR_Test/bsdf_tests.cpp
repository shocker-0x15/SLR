//
//  bsdf_tests.cpp
//
//  Created by 渡部 心 on 2017/05/24.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include <gtest/gtest.h>

#include <libSLR/MemoryAllocators/ArenaAllocator.h>
#include <libSLR/BasicTypes/spectrum_types.h>
#include <libSLR/BasicTypes/spectrum_library.h>
#include <libSLR/Core/distributions.h>
#include <libSLR/RNG/XORShiftRNG.h>

#include <libSLR/BSDF/bsdf_headers.h>

#define EXPECT_SAMPLED_SPECTRUM_NEAR(val1, val2, abs_error) \
do { \
for (int i = 0; i < SLR::SampledSpectrum::NumComponents; ++i) \
EXPECT_NEAR((val1)[i], (val2)[i], (abs_error)[i]); \
} while (0) \

typedef std::shared_ptr<SLR::AssetSpectrum> AssetSpectrumRef;

static AssetSpectrumRef createSpectrumFromData(const SLR::SpectrumLibrary::Data &data) {
    using namespace SLR;
    
    AssetSpectrumRef spectrum;
    if (data.dType == SpectrumLibrary::DistributionType::Regular)
        spectrum = createShared<RegularContinuousSpectrum>(data.minLambdas, data.maxLambdas, data.values, data.numSamples);
    else if (data.dType == SpectrumLibrary::DistributionType::Irregular)
        spectrum = createShared<IrregularContinuousSpectrum>(data.lambdas, data.values, data.numSamples);
    
    return spectrum;
}

struct CommonDataForBSDFTest {
    AssetSpectrumRef FlatReflectance;
    AssetSpectrumRef ColorCheckerReflectances[24];
    
    CommonDataForBSDFTest() {
        using namespace SLR;
        // FlatReflectance
        {
            float values[2] = {1.0f, 1.0f};
            FlatReflectance = createShared<RegularContinuousSpectrum>(WavelengthLowBound, WavelengthHighBound, values, lengthof(values));
        }
        // Color Checker
        for (int i = 0; i < 24; ++i) {
            SpectrumLibrary::Data data;
            SpectrumLibrary::queryReflectanceSpectrum("Color Checker", i, &data);
            
            ColorCheckerReflectances[i] = createSpectrumFromData(data);
        }
    }
};

static CommonDataForBSDFTest s_commonDataForBSDFTest;

static void BSDFTest_Function(const SLR::BSDF* bsdf, const SLR::WavelengthSamples &wls, SLR::RandomNumberGenerator* rng, bool adjoint, bool lowerIncidence) {
    using namespace SLR;
    
    const uint32_t NumZenithSectors = 64;
    const uint32_t NumAzimuthSectors = 128;
    const float ZenithSectorAngle = M_PI / NumZenithSectors;
    const float AzimuthSectorAngle = 2 * M_PI / NumAzimuthSectors;
    uint32_t* sampleCounts = new uint32_t[NumZenithSectors * NumAzimuthSectors];
    
    bool hasDelta = bsdf->hasDelta();
    
    Vector3D incidentVectors[] = {
        Vector3D(0, 0, 1), // normal incident
        normalize(Vector3D(1, 1, 1)), // ordinary incident
        normalize(Vector3D(1, 0, 0.001)), // grazing angle in tangent direction
        normalize(Vector3D(0, 1, 0.001)), // grazing angle in bitangent direction
        //        Vector3D(1, 0, 0) // extremely grazing angle
    };
    Normal3D geometricNormals[] = {
        Normal3D(0, 0, 1), // no perturbation
        normalize(Normal3D(0.25f, 0.25f, 1)), // weak perturbation
        normalize(Normal3D(-1, -1, 1)), // strong perturbation
    };
    for (int vIdx = 0; vIdx < lengthof(incidentVectors); ++vIdx) {
        Vector3D dirOutLocal = incidentVectors[vIdx];
        if (lowerIncidence)
            dirOutLocal.z *= -1;
        for (int gnIdx = 0; gnIdx < lengthof(geometricNormals); ++gnIdx) {
            Normal3D geomNormalLocal = geometricNormals[gnIdx];
            BSDFQuery query(dirOutLocal, geomNormalLocal, wls.selectedLambdaIndex, DirectionType::All, true, adjoint);
            
            std::fill(sampleCounts, sampleCounts + NumZenithSectors * NumAzimuthSectors, 0);
            
            SampledSpectrumSum accumEnergyThroughput = SampledSpectrum::Zero;
            
            FloatSum accumError_FwdSampled_vs_FwdEvaluated = 0;
            FloatSum accumError_FwdSampledPDF_vs_FwdEvaluatedPDF = 0;
            FloatSum accumError_FwdSampledRev_vs_FwdEvaluatedRev = 0;
            FloatSum accumError_FwdSampledPDFRev_vs_FwdEvaluatedPDFRev = 0;
            
            FloatSum accumError_RevEvaluated_vs_FwdEvaluatedRev = 0;
            FloatSum accumError_RevEvaluatedPDF_vs_FwdEvaluatedPDFRev = 0;
            FloatSum accumError_RevEvaluatedRev_vs_FwdEvaluated = 0;
            FloatSum accumError_RevEvaluatedPDFRev_vs_FwdEvaluatedPDF = 0;
            
            // Percentage Difference
            const auto calcErrorScore = [](float a, float b) {
                SLRAssert(a >= 0 && b >= 0, "Assuming a and b as positive values.");
                float avg = 0.5f * (a + b);
                return avg > 0 ? (a - b) * (a - b) / (avg * avg) : 0; 
            };
            
            uint32_t numValidSamples = 0;
            const uint32_t NumSamplesPerIncident = 100000;
            for (int j = 0; j < NumSamplesPerIncident; ++j) {
                BSDFSample sample(rng->getFloat0cTo1o(), rng->getFloat0cTo1o(), rng->getFloat0cTo1o());
                BSDFQueryResult result;
                
                SampledSpectrum fwdSampledValue = bsdf->sample(query, sample, &result);
                if (result.dirPDF == 0.0f)
                    continue;
                ++numValidSamples;
                
                EXPECT_TRUE(fwdSampledValue.allFinite() && !fwdSampledValue.hasNegative());
                
                accumEnergyThroughput += fwdSampledValue * (absDot(result.dirLocal, geomNormalLocal) / result.dirPDF);
                
                // JP: デルタ関数を含むBSDFの場合、明示的なサンプリング以外では値を評価できないためスキップする。
                // EN: Skip in the case for a BSDF containing the delta function because it is impossible to evaluate a value without explicit sampling.
                if (hasDelta)
                    continue;
                
                // JP: サンプルされた値がサンプル方向を使って再評価したBSDFの値と一致するかを確かめる。
                // EN: check if sampled values are consistent with values from re-evaluating BSDF using the sampled direction.
                SampledSpectrum fwdEvaluatedValue, fwdEvaluatedValueRev;
                float fwdEvaluatedPDFValue, fwdEvaluatedPDFValueRev;
                fwdEvaluatedValue = bsdf->evaluate(query, result.dirLocal, &fwdEvaluatedValueRev);
                fwdEvaluatedPDFValue = bsdf->evaluatePDF(query, result.dirLocal, &fwdEvaluatedPDFValueRev);
                
                for (int wl = 0; wl < SampledSpectrum::NumComponents; ++wl)
                    accumError_FwdSampled_vs_FwdEvaluated += calcErrorScore(fwdSampledValue[wl], fwdEvaluatedValue[wl]);
                accumError_FwdSampledPDF_vs_FwdEvaluatedPDF += calcErrorScore(result.dirPDF, fwdEvaluatedPDFValue);
                for (int wl = 0; wl < SampledSpectrum::NumComponents; ++wl)
                    accumError_FwdSampledRev_vs_FwdEvaluatedRev += calcErrorScore(result.reverse.value[wl], fwdEvaluatedValueRev[wl]);
                accumError_FwdSampledPDFRev_vs_FwdEvaluatedPDFRev += calcErrorScore(result.reverse.dirPDF, fwdEvaluatedPDFValueRev);
                
                BSDFQuery revQuery(result.dirLocal, geomNormalLocal, wls.selectedLambdaIndex, DirectionType::All, true, !adjoint);
                
                // JP: 逆方向から評価した値が、順方向のものと一致するかを確かめる。
                // EN: check if evaluated values from the reverse direction are consistent with the forward direction.
                SampledSpectrum revEvaluatedValue, revEvaluatedValueRev;
                float revEvaluatedPDFValue, revEvaluatedPDFValueRev;
                revEvaluatedValue = bsdf->evaluate(revQuery, dirOutLocal, &revEvaluatedValueRev);
                revEvaluatedPDFValue = bsdf->evaluatePDF(revQuery, dirOutLocal, &revEvaluatedPDFValueRev);
                
                for (int wl = 0; wl < SampledSpectrum::NumComponents; ++wl)
                    accumError_RevEvaluated_vs_FwdEvaluatedRev += calcErrorScore(revEvaluatedValue[wl], fwdEvaluatedValueRev[wl]);
                accumError_RevEvaluatedPDF_vs_FwdEvaluatedPDFRev += calcErrorScore(revEvaluatedPDFValue, fwdEvaluatedPDFValueRev);
                for (int wl = 0; wl < SampledSpectrum::NumComponents; ++wl)
                    accumError_RevEvaluatedRev_vs_FwdEvaluated += calcErrorScore(revEvaluatedValueRev[wl], fwdEvaluatedValue[wl]);
                accumError_RevEvaluatedPDFRev_vs_FwdEvaluatedPDF += calcErrorScore(revEvaluatedPDFValueRev, fwdEvaluatedPDFValue);
                
                // JP: サンプル方向の該当するビンの数をインクリメントする。
                // EN: increment the number of the bin for the sampled direction.
                float theta, phi;
                result.dirLocal.toPolarZUp(&theta, &phi);
                uint32_t zenithSectorIdx = std::min((uint32_t)(theta / M_PI * NumZenithSectors), NumZenithSectors - 1);
                uint32_t azimuthSectorIdx = (uint32_t)std::fmod(phi / (2 * M_PI) * NumAzimuthSectors, NumAzimuthSectors);
                ++sampleCounts[zenithSectorIdx * NumAzimuthSectors + azimuthSectorIdx];   
            }
            
//            // JP: BSDFのウェイトを評価する。
//            // EN: evaluate a BSDF weight.
//            SampledSpectrum energyThroughput = accumEnergyThroughput.result / numValidSamples;
//            float iEnergyThroughput = energyThroughput.importance(wls.selectedLambdaIndex);
//            float weight = bsdf->weight(query);
////            if (!(std::abs((iEnergyThroughput - weight) / weight) < 0.5 || weight == 0))
////                printf("");
//            EXPECT_TRUE(std::abs((iEnergyThroughput - weight) / weight) < 0.5 || weight == 0);
            
            if (hasDelta)
                continue;
            
            
            
            // RMSPD: Root Mean Squared Percentage Differences
            float RMSPD_FwdSampled_vs_FwdEvaluated = std::sqrt(accumError_FwdSampled_vs_FwdEvaluated.result / (numValidSamples * SampledSpectrum::NumComponents));
            float RMSPD_FwdSampledPDF_vs_FwdEvaluatedPDF = std::sqrt(accumError_FwdSampledPDF_vs_FwdEvaluatedPDF.result / numValidSamples);
            float RMSPD_FwdSampledRev_vs_FwdEvaluatedRev = std::sqrt(accumError_FwdSampledRev_vs_FwdEvaluatedRev.result / (numValidSamples * SampledSpectrum::NumComponents));
            float RMSPD_FwdSampledPDFRev_vs_FwdEvaluatedPDFRev = std::sqrt(accumError_FwdSampledPDFRev_vs_FwdEvaluatedPDFRev.result / numValidSamples);
            
            float RMSPD_RevEvaluated_vs_FwdEvaluatedRev = std::sqrt(accumError_RevEvaluated_vs_FwdEvaluatedRev.result / (numValidSamples * SampledSpectrum::NumComponents));
            float RMSPD_RevEvaluatedPDF_vs_FwdEvaluatedPDFRev = std::sqrt(accumError_RevEvaluatedPDF_vs_FwdEvaluatedPDFRev.result / numValidSamples);
            float RMSPD_RevEvaluatedRev_vs_FwdEvaluated = std::sqrt(accumError_RevEvaluatedRev_vs_FwdEvaluated.result / (numValidSamples * SampledSpectrum::NumComponents));
            float RMSPD_RevEvaluatedPDFRev_vs_FwdEvaluatedPDF = std::sqrt(accumError_RevEvaluatedPDFRev_vs_FwdEvaluatedPDF.result / numValidSamples);
            
            EXPECT_LT(RMSPD_FwdSampled_vs_FwdEvaluated, 0.01);
            EXPECT_LT(RMSPD_FwdSampledPDF_vs_FwdEvaluatedPDF, 0.01);
            EXPECT_LT(RMSPD_FwdSampledRev_vs_FwdEvaluatedRev, 0.01);
            EXPECT_LT(RMSPD_FwdSampledPDFRev_vs_FwdEvaluatedPDFRev, 0.01);
            
            EXPECT_LT(RMSPD_RevEvaluated_vs_FwdEvaluatedRev, 0.01);
            EXPECT_LT(RMSPD_RevEvaluatedPDF_vs_FwdEvaluatedPDFRev, 0.01);
            EXPECT_LT(RMSPD_RevEvaluatedRev_vs_FwdEvaluated, 0.01);
            EXPECT_LT(RMSPD_RevEvaluatedPDFRev_vs_FwdEvaluatedPDF, 0.01);
            
            
            
            // JP: 確率密度関数の正当性を確認する。
            // EN: validate the correctness of the probability density function.
            query.requestReverse = false;
            uint32_t sumCounts = 0;
            FloatSum accumError = 0;
            for (int z = 0; z < NumZenithSectors; ++z) {
                for (int a = 0; a < NumAzimuthSectors; ++a) {
                    float theta = M_PI * (z + 0.5f) / NumZenithSectors;
                    float phi = 2 * M_PI * (a + 0.5f) / NumAzimuthSectors;
                    Vector3D dirInLocal = Vector3D::fromPolarZUp(phi, theta);
                    
                    uint32_t count = sampleCounts[z * NumAzimuthSectors + a];
                    float sinThetaLow = std::sin(M_PI * z / NumZenithSectors);
                    float sinThetaHigh = std::sin(M_PI * (z + 1) / NumZenithSectors);
                    float approxSolidAngle = 0.5f * (sinThetaLow + sinThetaHigh) * AzimuthSectorAngle * ZenithSectorAngle;
                    float measuredProb = (float)count / NumSamplesPerIncident;
                    float analyticProb = bsdf->evaluatePDF(query, dirInLocal) * approxSolidAngle;
                    
                    accumError += (measuredProb - analyticProb) * (measuredProb - analyticProb);
                    
                    sumCounts += count;
                }
            }
            
            // JP: 誤差はサンプル数の平方根に反比例して減少する。
            // EN: Error decreases in inverse proportion to the square root of the sample count.
            const float AcceptableRelativeErrorScale = 0.3f;
            float RMSE = std::sqrt(accumError.result / sumCounts);
            EXPECT_LT(RMSE, AcceptableRelativeErrorScale / std::sqrt(NumSamplesPerIncident));
            ASSERT_EQ(sumCounts, numValidSamples);   
        }
    }
    
    delete[] sampleCounts;
}

TEST(BSDFTest, LambertianBRDF) {
    using namespace SLR;
    
    ArenaAllocator mem;
    RandomNumberGenerator* rng = new XORShiftRNG(1592814120);
    
    const uint32_t NumWavelengthSamples = 10;
    for (int j = 0; j < NumWavelengthSamples; ++j) {
        float wlPDF;
        WavelengthSamples wls = WavelengthSamples::createWithEqualOffsets((j + 0.5f) / NumWavelengthSamples, rng->getFloat0cTo1o(), &wlPDF);
        
        std::shared_ptr<BSDF> bsdf = createShared<LambertianBRDF>(s_commonDataForBSDFTest.FlatReflectance->evaluate(wls));
        BSDFTest_Function(bsdf.get(), wls, rng, false, false);
        BSDFTest_Function(bsdf.get(), wls, rng, true, false);
        BSDFTest_Function(bsdf.get(), wls, rng, false, true);
        BSDFTest_Function(bsdf.get(), wls, rng, true, true);
    }
    
    delete rng;
}

TEST(BSDFTest, SpecularBRDF) {
    using namespace SLR;
    
    ArenaAllocator mem;
    RandomNumberGenerator* rng = new XORShiftRNG(1592814120);
    
    AssetSpectrumRef eta, k;
    {
        SpectrumLibrary::Data data;
        
        SpectrumLibrary::queryIoRSpectrum("Aluminium", 0, &data);
        eta = createSpectrumFromData(data);
        SpectrumLibrary::queryIoRSpectrum("Aluminium", 1, &data);
        k = createSpectrumFromData(data);
    }
    
    const uint32_t NumWavelengthSamples = 10;
    for (int j = 0; j < NumWavelengthSamples; ++j) {
        float wlPDF;
        WavelengthSamples wls = WavelengthSamples::createWithEqualOffsets((j + 0.5f) / NumWavelengthSamples, rng->getFloat0cTo1o(), &wlPDF);
        
        std::shared_ptr<BSDF> bsdf = createShared<SpecularBRDF>(s_commonDataForBSDFTest.FlatReflectance->evaluate(wls), eta->evaluate(wls), k->evaluate(wls));
        BSDFTest_Function(bsdf.get(), wls, rng, false, false);
        BSDFTest_Function(bsdf.get(), wls, rng, true, false);
        BSDFTest_Function(bsdf.get(), wls, rng, false, true);
        BSDFTest_Function(bsdf.get(), wls, rng, true, true);
    }
    
    delete rng;
}

TEST(BSDFTest, SpecularBSDF) {
    using namespace SLR;
    
    ArenaAllocator mem;
    RandomNumberGenerator* rng = new XORShiftRNG(1592814120);
    
    AssetSpectrumRef etaExt, etaInt;
    {
        SpectrumLibrary::Data data;
        
        SpectrumLibrary::queryIoRSpectrum("Air", 0, &data);
        etaExt = createSpectrumFromData(data);
        SpectrumLibrary::queryIoRSpectrum("Diamond", 0, &data);
        etaInt = createSpectrumFromData(data);
    }
    
    const uint32_t NumWavelengthSamples = 10;
    for (int j = 0; j < NumWavelengthSamples; ++j) {
        float wlPDF;
        WavelengthSamples wls = WavelengthSamples::createWithEqualOffsets((j + 0.5f) / NumWavelengthSamples, rng->getFloat0cTo1o(), &wlPDF);
        
        std::shared_ptr<BSDF> bsdf = createShared<SpecularBSDF>(s_commonDataForBSDFTest.FlatReflectance->evaluate(wls), etaExt->evaluate(wls), etaInt->evaluate(wls), true);
        BSDFTest_Function(bsdf.get(), wls, rng, false, false);
        BSDFTest_Function(bsdf.get(), wls, rng, true, false);
        BSDFTest_Function(bsdf.get(), wls, rng, false, true);
        BSDFTest_Function(bsdf.get(), wls, rng, true, true);
    }
    
    delete rng;
}

TEST(BSDFTest, OrenNayerBRDF) {
    using namespace SLR;
    
    ArenaAllocator mem;
    RandomNumberGenerator* rng = new XORShiftRNG(1592814120);
    
    const uint32_t NumWavelengthSamples = 10;
    for (int j = 0; j < NumWavelengthSamples; ++j) {
        float wlPDF;
        WavelengthSamples wls = WavelengthSamples::createWithEqualOffsets((j + 0.5f) / NumWavelengthSamples, rng->getFloat0cTo1o(), &wlPDF);
        
        std::shared_ptr<BSDF> bsdf = createShared<OrenNayerBRDF>(s_commonDataForBSDFTest.FlatReflectance->evaluate(wls), 0.3f);
        BSDFTest_Function(bsdf.get(), wls, rng, false, false);
        BSDFTest_Function(bsdf.get(), wls, rng, true, false);
        BSDFTest_Function(bsdf.get(), wls, rng, false, true);
        BSDFTest_Function(bsdf.get(), wls, rng, true, true);
    }
    
    delete rng;
}

TEST(BSDFTest, ModifiedWardDurBRDF) {
    using namespace SLR;
    
    ArenaAllocator mem;
    RandomNumberGenerator* rng = new XORShiftRNG(1592814120);
    
    const uint32_t NumWavelengthSamples = 10;
    for (int j = 0; j < NumWavelengthSamples; ++j) {
        float wlPDF;
        WavelengthSamples wls = WavelengthSamples::createWithEqualOffsets((j + 0.5f) / NumWavelengthSamples, rng->getFloat0cTo1o(), &wlPDF);
        
        std::shared_ptr<BSDF> bsdf = createShared<ModifiedWardDurBRDF>(s_commonDataForBSDFTest.FlatReflectance->evaluate(wls), 0.4f, 0.05f);
        BSDFTest_Function(bsdf.get(), wls, rng, false, false);
        BSDFTest_Function(bsdf.get(), wls, rng, true, false);
        BSDFTest_Function(bsdf.get(), wls, rng, false, true);
        BSDFTest_Function(bsdf.get(), wls, rng, true, true);
    }
    
    delete rng;
}

TEST(BSDFTest, AshikhminShirleyBRDF) {
    using namespace SLR;
    
    ArenaAllocator mem;
    RandomNumberGenerator* rng = new XORShiftRNG(1592814120);
    
    const uint32_t NumWavelengthSamples = 10;
    for (int j = 0; j < NumWavelengthSamples; ++j) {
        float wlPDF;
        WavelengthSamples wls = WavelengthSamples::createWithEqualOffsets((j + 0.5f) / NumWavelengthSamples, rng->getFloat0cTo1o(), &wlPDF);
        
        SampledSpectrum Rs = s_commonDataForBSDFTest.FlatReflectance->evaluate(wls);
        SampledSpectrum Rd = s_commonDataForBSDFTest.FlatReflectance->evaluate(wls);
        std::shared_ptr<BSDF> bsdf = createShared<AshikhminShirleyBRDF>(Rs, Rd, 1000, 5);
        BSDFTest_Function(bsdf.get(), wls, rng, false, false);
        BSDFTest_Function(bsdf.get(), wls, rng, true, false);
        BSDFTest_Function(bsdf.get(), wls, rng, false, true);
        BSDFTest_Function(bsdf.get(), wls, rng, true, true);
    }
    
    delete rng;
}

TEST(BSDFTest, MicrofacetBRDF) {
    using namespace SLR;
    
    ArenaAllocator mem;
    RandomNumberGenerator* rng = new XORShiftRNG(1592814120);
    
    AssetSpectrumRef eta, k;
    {
        SpectrumLibrary::Data data;
        
        SpectrumLibrary::queryIoRSpectrum("Aluminium", 0, &data);
        eta = createSpectrumFromData(data);
        SpectrumLibrary::queryIoRSpectrum("Aluminium", 1, &data);
        k = createSpectrumFromData(data);
    }
    
    std::shared_ptr<MicrofacetDistribution> D = createShared<GGXMicrofacetDistribution>(0.05f, 0.2f); 
    
    const uint32_t NumWavelengthSamples = 10;
    for (int j = 0; j < NumWavelengthSamples; ++j) {
        float wlPDF;
        WavelengthSamples wls = WavelengthSamples::createWithEqualOffsets((j + 0.5f) / NumWavelengthSamples, rng->getFloat0cTo1o(), &wlPDF);
        
        std::shared_ptr<BSDF> bsdf = createShared<MicrofacetBRDF>(eta->evaluate(wls), k->evaluate(wls), D.get());
        BSDFTest_Function(bsdf.get(), wls, rng, false, false);
        BSDFTest_Function(bsdf.get(), wls, rng, true, false);
        BSDFTest_Function(bsdf.get(), wls, rng, false, true);
        BSDFTest_Function(bsdf.get(), wls, rng, true, true);
    }
    
    delete rng;
}

TEST(BSDFTest, MicrofacetBSDF) {
    using namespace SLR;
    
    ArenaAllocator mem;
    RandomNumberGenerator* rng = new XORShiftRNG(1592814120);
    
    AssetSpectrumRef etaExt, etaInt;
    {
        SpectrumLibrary::Data data;
        
        SpectrumLibrary::queryIoRSpectrum("Air", 0, &data);
        etaExt = createSpectrumFromData(data);
        SpectrumLibrary::queryIoRSpectrum("Diamond", 0, &data);
        etaInt = createSpectrumFromData(data);
    }
    
    std::shared_ptr<MicrofacetDistribution> D = createShared<GGXMicrofacetDistribution>(0.05f, 0.2f); 
    
    const uint32_t NumWavelengthSamples = 10;
    for (int j = 0; j < NumWavelengthSamples; ++j) {
        float wlPDF;
        WavelengthSamples wls = WavelengthSamples::createWithEqualOffsets((j + 0.5f) / NumWavelengthSamples, rng->getFloat0cTo1o(), &wlPDF);
        
        std::shared_ptr<BSDF> bsdf = createShared<MicrofacetBSDF>(etaExt->evaluate(wls), etaInt->evaluate(wls), D.get());
        BSDFTest_Function(bsdf.get(), wls, rng, false, false);
        BSDFTest_Function(bsdf.get(), wls, rng, true, false);
        BSDFTest_Function(bsdf.get(), wls, rng, false, true);
        BSDFTest_Function(bsdf.get(), wls, rng, true, true);
    }
    
    delete rng;
}

TEST(BSDFTest, DisneyBRDF) {
    using namespace SLR;
    
    ArenaAllocator mem;
    RandomNumberGenerator* rng = new XORShiftRNG(1592814120);
    
    AssetSpectrumRef etaExt, etaInt;
    {
        SpectrumLibrary::Data data;
        
        SpectrumLibrary::queryIoRSpectrum("Air", 0, &data);
        etaExt = createSpectrumFromData(data);
        SpectrumLibrary::queryIoRSpectrum("Diamond", 0, &data);
        etaInt = createSpectrumFromData(data);
    }
    
    std::shared_ptr<MicrofacetDistribution> D = createShared<GGXMicrofacetDistribution>(0.05f, 0.2f);
    
    float baseColorXYZ[3];
    s_commonDataForBSDFTest.FlatReflectance->convertToXYZ(baseColorXYZ);
    float baseColorLuminance = baseColorXYZ[1];
    
    const uint32_t NumWavelengthSamples = 10;
    for (int j = 0; j < NumWavelengthSamples; ++j) {
        float wlPDF;
        WavelengthSamples wls = WavelengthSamples::createWithEqualOffsets((j + 0.5f) / NumWavelengthSamples, rng->getFloat0cTo1o(), &wlPDF);
        
        SampledSpectrum baseColor = s_commonDataForBSDFTest.FlatReflectance->evaluate(wls);
        float subsurface = 0.5f;
        float metallic = 0.5f;
        float specular = 1.0f;
        float specularTint = 0.5f;
        float roughness = 0.5f;
        float anisotropic = 1.0f;
        float sheen = 1.0f;
        float sheenTint = 0.5f;
        float clearCoat = 1.0f;
        float clearCoatGloss = 0.5f;
        std::shared_ptr<BSDF> bsdf = createShared<DisneyBRDF>(baseColor, baseColorLuminance, 
                                                              subsurface, metallic, specular, specularTint, roughness, anisotropic, 
                                                              sheen, sheenTint, clearCoat, clearCoatGloss);
        BSDFTest_Function(bsdf.get(), wls, rng, false, false);
        BSDFTest_Function(bsdf.get(), wls, rng, true, false);
        BSDFTest_Function(bsdf.get(), wls, rng, false, true);
        BSDFTest_Function(bsdf.get(), wls, rng, true, true);
    }
    
    delete rng;
}

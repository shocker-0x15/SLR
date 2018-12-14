//
//  DensityGridMediumDistribution.cpp
//
//  Created by 渡部 心 on 2017/02/14.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "DensityGridMediumDistribution.h"

#include "../Core/light_path_sampler.h"

#define UseSuperVoxels
#define UseRatioTracking

namespace SLR {
    float DensityGridMediumDistribution::calcDensityInSuperVoxels(const Point3D &param) const {
        if (param.x < 0 || param.y < 0 || param.z < 0 ||
            param.x > 1 || param.y > 1 || param.z > 1)
            return 0.0f;
        
        uint32_t lx = std::min((uint32_t)(param.x * (m_svNumX - 1)), m_svNumX - 1);
        uint32_t ux = std::min(lx + 1, m_svNumX - 1);
        uint32_t ly = std::min((uint32_t)(param.y * (m_svNumY - 1)), m_svNumY - 1);
        uint32_t uy = std::min(ly + 1, m_svNumY - 1);
        uint32_t lz = std::min((uint32_t)(param.z * (m_svNumZ - 1)), m_svNumZ - 1);
        uint32_t uz = std::min(lz + 1, m_svNumZ - 1);
        float wux = param.x * (m_svNumX - 1) - lx;
        float wlx = 1 - wux;
        float wuy = param.y * (m_svNumY - 1) - ly;
        float wly = 1 - wuy;
        float wuz = param.z * (m_svNumZ - 1) - lz;
        float wlz = 1 - wuz;
        float density = (m_superVoxels[m_svNumX * m_svNumY * lz + m_svNumX * ly + lx] * wlz * wly * wlx +
                         m_superVoxels[m_svNumX * m_svNumY * lz + m_svNumX * ly + ux] * wlz * wly * wux +
                         m_superVoxels[m_svNumX * m_svNumY * lz + m_svNumX * uy + lx] * wlz * wuy * wlx +
                         m_superVoxels[m_svNumX * m_svNumY * lz + m_svNumX * uy + ux] * wlz * wuy * wux +
                         m_superVoxels[m_svNumX * m_svNumY * uz + m_svNumX * ly + lx] * wuz * wly * wlx +
                         m_superVoxels[m_svNumX * m_svNumY * uz + m_svNumX * ly + ux] * wuz * wly * wux +
                         m_superVoxels[m_svNumX * m_svNumY * uz + m_svNumX * uy + lx] * wuz * wuy * wlx +
                         m_superVoxels[m_svNumX * m_svNumY * uz + m_svNumX * uy + ux] * wuz * wuy * wux);
        uint32_t diffCellIdx = (m_svNumX - 1) * (m_svNumY - 1) * std::min(lz, m_svNumZ - 2) + (m_svNumX - 1) * std::min(ly, m_svNumY - 2) + std::min(lx, m_svNumX - 2);
        return density + m_maximumDifferences[diffCellIdx];
    };
    
    void DensityGridMediumDistribution::setupSuperVoxels() {
        m_svNumX = std::max(m_numX / 16, 4u);
        m_svNumY = std::max(m_numY / 16, 4u);
        m_svNumZ = std::max(m_numZ / 16, 4u);
        
        m_superVoxels = new float[m_svNumX * m_svNumY * m_svNumZ];
        const uint32_t LengthOfCompensationCells = (m_svNumX - 1) * (m_svNumY - 1) * (m_svNumZ - 1);
        m_maximumDifferences = new float[LengthOfCompensationCells];
        float* maximumDifferences = new float[LengthOfCompensationCells];
        std::fill(m_maximumDifferences, m_maximumDifferences + LengthOfCompensationCells, 0.0f);
        std::fill(maximumDifferences, maximumDifferences + LengthOfCompensationCells, 0.0f);
        m_superVoxelWidth = (m_region.maxP - m_region.minP) / Vector3D(m_svNumX - 1, m_svNumY - 1, m_svNumZ - 1);
        
        // JP: スーパーボクセルの角においてオリジナルの値を評価してスーパーボクセルの値を初期化する。
        // EN: initialize super voxel values by original values evaluated at the corners of super voxels.
        for (int siz = 0; siz < m_svNumZ; ++siz) {
            float pz = (float)siz / (m_svNumZ - 1);
            for (int siy = 0; siy < m_svNumY; ++siy) {
                float py = (float)siy / (m_svNumY - 1); 
                for (int six = 0; six < m_svNumX; ++six) {
                    float px = (float)six / (m_svNumX - 1);
                    uint32_t idx = siz * m_svNumY * m_svNumX + siy * m_svNumX + six;
                    m_superVoxels[idx] = calcDensity(Point3D(px, py, pz));
                }
            }
        }

        // JP: 各スーパーボクセルにおける最大の差を計算して、スーパーボクセル中のmajorant extinctionが
        //     upper-boundingになるようにする。
        // EN: calculate maximum difference for each super voxel to make sure that 
        //     the majorant extinction in the super voxel is upper bounding.
        for (int svCellZ = 0; svCellZ < m_svNumZ - 1; ++svCellZ) {
            float lpz = (float)svCellZ / (m_svNumZ - 1);
            float upz = (float)(svCellZ + 1) / (m_svNumZ - 1);
            uint32_t liz = (uint32_t)std::floor(lpz * (m_numZ - 1));
            uint32_t uiz = (uint32_t)std::ceil(upz * (m_numZ - 1));
            
            for (int svCellY = 0; svCellY < m_svNumY - 1; ++svCellY) {
                float lpy = (float)svCellY / (m_svNumY - 1);
                float upy = (float)(svCellY + 1) / (m_svNumY - 1);
                uint32_t liy = (uint32_t)std::floor(lpy * (m_numY - 1));
                uint32_t uiy = (uint32_t)std::ceil(upy * (m_numY - 1));
                
                for (int svCellX = 0; svCellX < m_svNumX - 1; ++svCellX) {
                    float lpx = (float)svCellX / (m_svNumX - 1);
                    float upx = (float)(svCellX + 1) / (m_svNumX - 1);
                    uint32_t lix = (uint32_t)std::floor(lpx * (m_numX - 1));
                    uint32_t uix = (uint32_t)std::ceil(upx * (m_numX - 1));
                    
                    float maxDiff = 0;
                    for (int iz = liz; iz <= uiz; ++iz) {
                        float pz = (float)iz / (m_numZ - 1);
                        for (int iy = liy; iy <= uiy; ++iy) {
                            float py = (float)iy / (m_numY - 1);
                            for (int ix = lix; ix <= uix; ++ix) {
                                float px = (float)ix / (m_numX - 1);
                                if (iz > liz && iz < uiz && iy > liy && iy < uiy && ix > lix && ix < uix) {
                                    // JP: オリジナルの値のサンプル点における差を計算する。
                                    // EN: calculate the difference at the sampling point of original values. 
                                    Point3D p(px, py, pz);
                                    float actualDensity = m_density_grid[iz][m_numX * iy + ix];
                                    float coarseDensity = calcDensityInSuperVoxels(p);
                                    float diff = actualDensity - coarseDensity;
                                    maxDiff = std::max(maxDiff, diff);   
                                }
                                else {
                                    // JP: 補間されたオリジナルの値がスーパーボクセル境界でスーパーボクセルの値を超える可能性がある。
                                    // EN: There is a possibility that an interpolated original value exceeds that of super voxel values at super voxel boundaries.
                                    int zBase = iz;
                                    int yBase = iy;
                                    int xBase = ix;
                                    int numZIterations = (iz > liz && iz < uiz) ? 1 : 2;
                                    int numYIterations = (iy > liy && iy < uiy) ? 1 : 2;
                                    int numXIterations = (ix > lix && ix < uix) ? 1 : 2;
                                    if (numZIterations == 2)
                                        zBase = iz == liz ? liz : uiz - 1;
                                    if (numYIterations == 2)
                                        yBase = iy == liy ? liy : uiy - 1;
                                    if (numXIterations == 2)
                                        xBase = ix == lix ? lix : uix - 1;
                                    for (int diz = 0; diz < numZIterations; ++diz) {
                                        int iiz = zBase + diz;
                                        float ppz = std::clamp((float)iiz / (m_numZ - 1), lpz, upz);
                                        for (int diy = 0; diy < numYIterations; ++diy) {
                                            int iiy = yBase + diy;
                                            float ppy = std::clamp((float)iiy / (m_numY - 1), lpy, upy);
                                            for (int dix = 0; dix < numXIterations; ++dix) {
                                                int iix = xBase + dix;
                                                float ppx = std::clamp((float)iix / (m_numX - 1), lpx, upx);
                                                Point3D p(ppx, ppy, ppz);
                                                float actualDensity = calcDensity(p);
                                                float coarseDensity = calcDensityInSuperVoxels(p);
                                                float diff = actualDensity - coarseDensity;
                                                maxDiff = std::max(maxDiff, diff);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                    uint32_t svCellIndex = (m_svNumX - 1) * (m_svNumY - 1) * svCellZ + (m_svNumX - 1) * svCellY + svCellX;
                    maximumDifferences[svCellIndex] = maxDiff;
                    SLRAssert(svCellIndex < LengthOfCompensationCells, "Invalid index.");
                    SLRAssert(maxDiff >= 0, "maximum difference must be greater than or equal to 0.");
                }
            }
        }

        std::copy(maximumDifferences, maximumDifferences + LengthOfCompensationCells, m_maximumDifferences);
        delete[] maximumDifferences;

//        // sanity check
//        slrDevPrintf("start sanity check.\n");
//        for (int siz = 0; siz < m_svNumZ; ++siz) {
//            float pz = (float)siz / (m_svNumZ - 1);
//            for (int siy = 0; siy < m_svNumY; ++siy) {
//                float py = (float)siy / (m_svNumY - 1); 
//                for (int six = 0; six < m_svNumX; ++six) {
//                    float px = (float)six / (m_svNumX - 1);
//                    float coarseDensity = calcDensityInSuperVoxels(Point3D(px, py, pz));
//                    if (coarseDensity < 0.0f)
//                        slrDevPrintf("%g at (%g, %g, %g)\n", coarseDensity, px, py, pz);
//                }
//            }
//        }
//        FloatSum avgNaiveRatio = 0;
//        FloatSum avgRatio = 0;
//        uint32_t numValidPoints = 0;
//        for (int iz = 0; iz < m_numZ * 4; ++iz) {
//            float pz = (float)iz / (m_numZ * 4 - 1);
//            for (int iy = 0; iy < m_numY * 4; ++iy) {
//                float py = (float)iy / (m_numY * 4 - 1);
//                for (int ix = 0; ix < m_numX * 4; ++ix) {
//                    float px = (float)ix / (m_numX * 4 - 1);
//                    Point3D param(px, py, pz);
//                    float actualDensity = calcDensity(param);
//                    float coarseDensity = calcDensityInSuperVoxels(param);
//                    if (actualDensity > coarseDensity || (coarseDensity == 0.0f && actualDensity > 0))
//                        slrDevPrintf("%g, %g / %g at (%g, %g, %g)\n", actualDensity / coarseDensity, actualDensity, coarseDensity, 
//                                     param.x * (m_svNumX - 1), param.y * (m_svNumY - 1), param.z * (m_svNumZ - 1));
//                    avgNaiveRatio += actualDensity / m_maxDensity;
//                    if (coarseDensity > 0.0f) {
//                        avgRatio += actualDensity / coarseDensity;
//                        ++numValidPoints;
//                    }
//                    else if (coarseDensity == 0.0f && actualDensity == 0.0f) {
//                        avgRatio += 1;
//                        ++numValidPoints;
//                    }
////                    SLRAssert(actualDensity <= coarseDensity, 
////                              "The value in the super voxel must be equal or greater than the actual value.\n"
////                              "%g / %g, %g", actualDensity, coarseDensity, actualDensity / coarseDensity);
//                }
//            }
//        }
//        slrDevPrintf("Avg Ratio: %g, (%g)\n", avgRatio / numValidPoints, avgNaiveRatio / (m_numZ * m_numY * m_numX * 64));
//        slrDevPrintf("");
    }
    
    // References:
    // Free Path Sampling in High Resolution Inhomogeneous Participating Media
    bool DensityGridMediumDistribution::traverseSuperVoxels(const Ray &ray, const RaySegment &segment, FreePathSampler &sampler, float base_sigma_e, 
                                                            const int32_t step[3], const float delta_t[3], const int32_t outsideIndices[3],   
                                                            float max_t[3], int32_t superVoxel[3], 
                                                            float* sampledDistance, float* majorantAtScattering) const {
        // JP: 多項式の係数を求める。
        // EN: compute polynomial coefficients.
        const auto polyCoeff = [this](const Point3D &entryPoint, const Point3D &exitPoint, const int32_t cellIndices[3], float coeffs[4]) {
            float density000 = m_superVoxels[m_svNumX * m_svNumY * (cellIndices[2] + 0) + m_svNumX * (cellIndices[1] + 0) + (cellIndices[0] + 0)];
            float density100 = m_superVoxels[m_svNumX * m_svNumY * (cellIndices[2] + 0) + m_svNumX * (cellIndices[1] + 0) + (cellIndices[0] + 1)];
            float density010 = m_superVoxels[m_svNumX * m_svNumY * (cellIndices[2] + 0) + m_svNumX * (cellIndices[1] + 1) + (cellIndices[0] + 0)];
            float density110 = m_superVoxels[m_svNumX * m_svNumY * (cellIndices[2] + 0) + m_svNumX * (cellIndices[1] + 1) + (cellIndices[0] + 1)];
            float density001 = m_superVoxels[m_svNumX * m_svNumY * (cellIndices[2] + 1) + m_svNumX * (cellIndices[1] + 0) + (cellIndices[0] + 0)];
            float density101 = m_superVoxels[m_svNumX * m_svNumY * (cellIndices[2] + 1) + m_svNumX * (cellIndices[1] + 0) + (cellIndices[0] + 1)];
            float density011 = m_superVoxels[m_svNumX * m_svNumY * (cellIndices[2] + 1) + m_svNumX * (cellIndices[1] + 1) + (cellIndices[0] + 0)];
            float density111 = m_superVoxels[m_svNumX * m_svNumY * (cellIndices[2] + 1) + m_svNumX * (cellIndices[1] + 1) + (cellIndices[0] + 1)];
            
            Point3D baseCoordinates = m_region.minP + Vector3D(cellIndices[0], cellIndices[1], cellIndices[2]) * m_superVoxelWidth;
            Point3D epInUnitCube = (entryPoint - baseCoordinates) / m_superVoxelWidth;
            Point3D opInUnitCube = (exitPoint - baseCoordinates) / m_superVoxelWidth;
            Vector3D delta = opInUnitCube - epInUnitCube;
            
            float dxyz = density111 - density011 - density101 - density110 + density001 + density010 + density100 - density000;
            float dxy = density000 - density100 - density010 + density110;
            float dxz = density000 - density100 - density001 + density101;
            float dyz = density000 - density010 - density001 + density011;
            float dx = density100 - density000;
            float dy = density010 - density000;
            float dz = density001 - density000;
            
            coeffs[3] = dxyz * delta.x * delta.y * delta.z;
            coeffs[2] = ((epInUnitCube.z * delta.x * delta.y + epInUnitCube.y * delta.x * delta.z + epInUnitCube.x * delta.y * delta.z) * dxyz + 
                         dxy * delta.x * delta.y + dxz * delta.x * delta.z + dyz * delta.y * delta.z);
            coeffs[1] = ((epInUnitCube.y * epInUnitCube.z * delta.x + epInUnitCube.x * epInUnitCube.z * delta.y + epInUnitCube.x * epInUnitCube.y * delta.z) * dxyz + 
                         dx * delta.x + dy * delta.y + dz * delta.z + 
                         (epInUnitCube.y * delta.x + epInUnitCube.x * delta.y) * dxy + 
                         (epInUnitCube.z * delta.x + epInUnitCube.x * delta.z) * dxz + 
                         (epInUnitCube.z * delta.y + epInUnitCube.y * delta.z) * dyz);
            coeffs[0] = (epInUnitCube.x * epInUnitCube.y * epInUnitCube.z * dxyz + 
                         epInUnitCube.x * epInUnitCube.y * dxy + 
                         epInUnitCube.x * epInUnitCube.z * dxz + 
                         epInUnitCube.y * epInUnitCube.z * dyz + 
                         epInUnitCube.x * dx + epInUnitCube.y * dy + epInUnitCube.z * dz + density000);
            float maximumDifference = m_maximumDifferences[(m_svNumX - 1) * (m_svNumY - 1) * cellIndices[2] + (m_svNumX - 1) * cellIndices[1] + cellIndices[0]];
            coeffs[0] += maximumDifference;
        };
        
        // JP: 散乱が起こる距離をregular falsi法を使って求める。
        // EN: find a distance at which scattering occurs by regula falsi method.
        const auto solvePolyInt = [](float base_sigma_e, const float coeffs[4], float entryDistance, float exitDistance, float entryOpticalDepth, float exitOpticalDepth, float sampledOpticalDepth, 
                                     float* sampledDistance, float* distParamInSuperVoxel) {
            float paramLow = 0.0f;
            float paramHigh = 1.0f;
            float opticalDepthLow = entryOpticalDepth;
            float opticalDepthHigh = exitOpticalDepth;
            uint32_t numIterations = 0;
            float prevParam = -INFINITY;
            float curParam = paramHigh;
            while (std::fabs(curParam - prevParam) >= 1e-3f * std::fabs(curParam) && numIterations < 100) {
                prevParam = curParam;
                curParam = paramLow + (paramHigh - paramLow) * (sampledOpticalDepth - opticalDepthLow) / (opticalDepthHigh - opticalDepthLow);
                float temp = base_sigma_e * 
                (coeffs[0] * curParam + 
                 coeffs[1] * curParam * curParam / 2 + 
                 coeffs[2] * curParam * curParam * curParam / 3 + 
                 coeffs[3] * curParam * curParam * curParam * curParam / 4);
                float opticalDepth = entryOpticalDepth + (exitDistance - entryDistance) * temp;
                
                ++numIterations;
                
                if (opticalDepth < sampledOpticalDepth) {
                    paramLow = curParam;
                    opticalDepthLow = opticalDepth;
                }
                else if (opticalDepth > sampledOpticalDepth) {
                    paramHigh = curParam;
                    opticalDepthHigh = opticalDepth;
                }
                else {
                    break;
                }
            }
            SLRAssert(std::isfinite(curParam), "Invalid value.\n"
                      "curParam: %g, paramLow/High: %g/%g, opticalDepthLow/High: %g/%g\n"
                      "base_sigma_e: %g, entry/exitDistance: %g/%g, entry/exit/sampledOpticalDepth: %g/%g/%g", 
                      curParam, paramLow, paramHigh, opticalDepthLow, opticalDepthHigh, 
                      base_sigma_e, entryDistance, exitDistance, entryOpticalDepth, exitOpticalDepth, sampledOpticalDepth);
            *distParamInSuperVoxel = curParam;
            *sampledDistance = entryDistance + (exitDistance - entryDistance) * curParam;
        };
        
        float sampledOpticalDepth = -std::log(1 - sampler.getSample());
        if (sampledOpticalDepth == 0.0f) {
            Point3D queryPoint = ray.org + *sampledDistance * ray.dir;
            Point3D param;
            m_region.calculateLocalCoordinates(queryPoint, &param);
            *majorantAtScattering = base_sigma_e * calcDensityInSuperVoxels(param);
            return true;
        }
        float entryDistance = 0.0f;
        float exitDistance = segment.distMin;
        Point3D entryPoint;
        Point3D exitPoint = ray.org + exitDistance * ray.dir;
        float entryOpticalDepth = 0.0f;
        float exitOpticalDepth = 0.0f;
        float coeffs[4];
        while (true) {
            entryPoint = exitPoint;
            entryDistance = exitDistance;
            entryOpticalDepth = exitOpticalDepth;
            
            const uint8_t axisMap[] = {2, 2, 2, 0, 1, 2, 1, 0};
            uint32_t condition = (((max_t[0] < max_t[1]) << 0) |
                                  ((max_t[0] < max_t[2]) << 1) |
                                  ((max_t[1] < max_t[2]) << 2));
            uint32_t steppingAxis = axisMap[condition];
            
            exitDistance = max_t[steppingAxis];
            exitPoint = ray.org + exitDistance * ray.dir;
            polyCoeff(entryPoint, exitPoint, superVoxel, coeffs);
            
            exitOpticalDepth += (exitDistance - entryDistance) * base_sigma_e * (coeffs[0] + coeffs[1] / 2 + coeffs[2] / 3 + coeffs[3] / 4);
            if (exitOpticalDepth >= sampledOpticalDepth)
                break;
            
            superVoxel[steppingAxis] += step[steppingAxis];
            
            // Out of Volume, No Scattering
            if (superVoxel[steppingAxis] == outsideIndices[steppingAxis] || exitDistance >= segment.distMax) {
                *sampledDistance = segment.distMax;
                *majorantAtScattering = NAN;
                return false;
            }
            SLRAssert(superVoxel[0] >= 0 && superVoxel[0] < m_svNumX - 1 && 
                      superVoxel[1] >= 0 && superVoxel[1] < m_svNumY - 1 && 
                      superVoxel[2] >= 0 && superVoxel[2] < m_svNumZ - 1, "Invalid Super Voxel.");
            
#if 1
            float nextPlaneCoord = m_region.minP[steppingAxis] + m_superVoxelWidth[steppingAxis] * (superVoxel[steppingAxis] + (step[steppingAxis] == 1 ? 1 : 0));
            max_t[steppingAxis] = (nextPlaneCoord - ray.org[steppingAxis]) / ray.dir[steppingAxis];
#else
            max_t[steppingAxis] += delta_t[steppingAxis];
#endif
        }
        
        float distParamInSuperVoxel;
        solvePolyInt(base_sigma_e, coeffs, entryDistance, exitDistance, entryOpticalDepth, exitOpticalDepth, sampledOpticalDepth, 
                     sampledDistance, &distParamInSuperVoxel);
        if (*sampledDistance >= segment.distMax) {
            *sampledDistance = segment.distMax;
            *majorantAtScattering = NAN;
            return false;
        }
        *majorantAtScattering = base_sigma_e * 
        (coeffs[0] + 
         coeffs[1] * distParamInSuperVoxel + 
         coeffs[2] * distParamInSuperVoxel * distParamInSuperVoxel + 
         coeffs[3] * distParamInSuperVoxel * distParamInSuperVoxel * distParamInSuperVoxel);
        
        return true;
    }

    float DensityGridMediumDistribution::calcDensity(const Point3D &param) const {
        if (param.x < 0 || param.y < 0 || param.z < 0 ||
            param.x > 1 || param.y > 1 || param.z > 1)
            return 0.0f;
        
        uint32_t lx = std::min((uint32_t)(param.x * (m_numX - 1)), m_numX - 1);
        uint32_t ux = std::min(lx + 1, m_numX - 1);
        uint32_t ly = std::min((uint32_t)(param.y * (m_numY - 1)), m_numY - 1);
        uint32_t uy = std::min(ly + 1, m_numY - 1);
        uint32_t lz = std::min((uint32_t)(param.z * (m_numZ - 1)), m_numZ - 1);
        uint32_t uz = std::min(lz + 1, m_numZ - 1);
        float wux = param.x * (m_numX - 1) - lx;
        float wlx = 1 - wux;
        float wuy = param.y * (m_numY - 1) - ly;
        float wly = 1 - wuy;
        float wuz = param.z * (m_numZ - 1) - lz;
        float wlz = 1 - wuz;
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
        
#if defined(UseSuperVoxels)
        // initialize 3D DDA process.
        const uint32_t MaxVoxelIndices[3] = {m_svNumX - 2, m_svNumY - 2, m_svNumZ - 2};
        Point3D initialPoint = ray.org + segment.distMin * ray.dir;
        Point3D localParam;
        m_region.calculateLocalCoordinates(initialPoint, &localParam);
        localParam = clamp(localParam, Point3D(0.0f), Point3D(1.0f));
        int32_t initSuperVoxel[3] = {
            std::min((int32_t)(localParam.x * (m_svNumX - 1)), (int32_t)m_svNumX - 2), 
            std::min((int32_t)(localParam.y * (m_svNumY - 1)), (int32_t)m_svNumY - 2), 
            std::min((int32_t)(localParam.z * (m_svNumZ - 1)), (int32_t)m_svNumZ - 2), 
        };
        int32_t step[3] = {0, 0, 0};
        float delta_t[3];
        int32_t outsideIndices[3] = {-1, -1, -1};
        float initMax_t[3];
        for (int i = 0; i < 3; ++i) {
            delta_t[i] = std::fabs(m_superVoxelWidth[i] / ray.dir[i]);
            if (ray.dir[i] > 0) {
                step[i] = 1;
                outsideIndices[i] = MaxVoxelIndices[i] + 1;
                float nextPlaneCoord = m_region.minP[i] + m_superVoxelWidth[i] * (initSuperVoxel[i] + 1);
                initMax_t[i] = (nextPlaneCoord - ray.org[i]) / ray.dir[i];
            }
            else {
                step[i] = -1;
                outsideIndices[i] = -1;
                float nextPlaneCoord = m_region.minP[i] + m_superVoxelWidth[i] * (initSuperVoxel[i] + 0);
                initMax_t[i] = (nextPlaneCoord - ray.org[i]) / ray.dir[i];
            }
        }
        
        // traverse super voxels
        bool hit = false;
        float extCoeffAtScattering = 0.0f;
        float sampledDistance = segment.distMin;
        int32_t superVoxel[3] = {initSuperVoxel[0], initSuperVoxel[1], initSuperVoxel[2]};
        float max_t[3] = {initMax_t[0], initMax_t[1], initMax_t[2]};
        while (sampledDistance < segment.distMax) {
            float majorantAtScattering;
            bool tentativeHit = traverseSuperVoxels(ray, RaySegment(sampledDistance, segment.distMax), sampler, base_sigma_e[wls.selectedLambdaIndex], 
                                                    step, delta_t, outsideIndices,  
                                                    max_t, superVoxel, 
                                                    &sampledDistance, &majorantAtScattering);
            if (tentativeHit) {
                Point3D queryPoint = ray.org + sampledDistance * ray.dir;
                Point3D param;
                m_region.calculateLocalCoordinates(queryPoint, &param);
                float extCoeff = base_sigma_e[wls.selectedLambdaIndex] * calcDensity(param);
                float probRealCollision = extCoeff / majorantAtScattering;
                SLRAssert(probRealCollision <= 1.0f, "Real extinction coefficient exceeds the majorant: %g > %g", extCoeff, majorantAtScattering);
                if (sampler.getSample() < probRealCollision) {
                    hit = true;
                    *mi = MediumInteraction(ray.time, sampledDistance, queryPoint, normalize(ray.dir), param.x, param.y, param.z);
                    extCoeffAtScattering = extCoeff;
                    break;
                }
            }
            else {
                hit = false;
                break;
            }
        }
        
        // estimate Monte Carlo throughput T(s, wl_j)/p(s, wl_i) by ratio tracking.
        *singleWavelength = false;
        if (wls.wavelengthSelected()) {
            *medThroughput = SampledSpectrum::Zero;
            (*medThroughput)[wls.selectedLambdaIndex] = 1.0f;
        }
        else {
            const float hitDistance = std::min(sampledDistance, segment.distMax);
            sampledDistance = segment.distMin;
            superVoxel[0] = initSuperVoxel[0];
            superVoxel[1] = initSuperVoxel[1];
            superVoxel[2] = initSuperVoxel[2];
            max_t[0] = initMax_t[0];
            max_t[1] = initMax_t[1];
            max_t[2] = initMax_t[2];
            
            SampledSpectrum mcThroughput = SampledSpectrum::One;
            for (int wl = 0; wl < WavelengthSamples::NumComponents; ++wl) {
                if (wl == wls.selectedLambdaIndex)
                    continue;
                
                while (sampledDistance < hitDistance) {
                    float majorantAtScattering;
                    bool tentativeHit = traverseSuperVoxels(ray, RaySegment(sampledDistance, hitDistance), sampler, base_sigma_e[wl], 
                                                            step, delta_t, outsideIndices,  
                                                            max_t, superVoxel, 
                                                            &sampledDistance, &majorantAtScattering);
                    if (tentativeHit) {
                        Point3D queryPoint = ray.org + sampledDistance * ray.dir;
                        Point3D param;
                        m_region.calculateLocalCoordinates(queryPoint, &param);
                        SampledSpectrum extCoeff = base_sigma_e * calcDensity(param);
                        float probRealCollision = (extCoeff[wl] - extCoeff[wls.selectedLambdaIndex]) / majorantAtScattering;
                        if (majorantAtScattering != 0.0f)
                            mcThroughput[wl] *= (1.0f - probRealCollision);
                    }
                }
            }
            *medThroughput = mcThroughput;
        }
        if (hit)
            *medThroughput /= extCoeffAtScattering;
        
        return hit;
#else
        // delta tracking to sample free path.
        float majorantSelected = majorantExtinctionCoefficientAtWavelength(wls.selectedWavelength());
        *singleWavelength = false;
        bool hit = false;
        float extCoeffAtScattering = 0.0f;
        FloatSum sampledDistance = segment.distMin;
        sampledDistance += -std::log(sampler.getSample()) / majorantSelected;
        while (sampledDistance < segment.distMax) {
            Point3D queryPoint = ray.org + sampledDistance * ray.dir;
            Point3D param;
            m_region.calculateLocalCoordinates(queryPoint, &param);
            float density = calcDensity(param);
            float extCoeff = base_sigma_e[wls.selectedLambdaIndex] * density;
            float probRealCollision = extCoeff / majorantSelected;
            if (sampler.getSample() < probRealCollision) {
                *mi = MediumInteraction(ray.time, sampledDistance, queryPoint, normalize(ray.dir), param.x, param.y, param.z);
                hit = true;
                extCoeffAtScattering = extCoeff;
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
            float hitDistance = std::min(sampledDistance.result, segment.distMax);
            SampledSpectrum mcThroughput = SampledSpectrum::One;
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
                    mcThroughput[wl] *= (1.0f - probRealCollision);
                    sampledDistance += -std::log(sampler.getSample()) / majorantWL;
                }
            }
            *medThroughput = mcThroughput;
        }
        if (hit)
            *medThroughput /= extCoeffAtScattering;
        
        return hit;
#endif
    }
    
    SampledSpectrum DensityGridMediumDistribution::evaluateTransmittance(const Ray &ray, const RaySegment &segment, const WavelengthSamples &wls, SLR::LightPathSampler &pathSampler, 
                                                                         bool* singleWavelength) const {
        SLRAssert(std::isfinite(segment.distMax), "distanceLimit must be a finite value.");
        FreePathSampler sampler = pathSampler.getFreePathSampler();
        
        SampledSpectrum base_sigma_e = m_base_sigma_e->evaluate(wls);
        
#if defined(UseSuperVoxels)
        // initialize 3D DDA process.
        const uint32_t MaxVoxelIndices[3] = {m_svNumX - 2, m_svNumY - 2, m_svNumZ - 2};
        Point3D initialPoint = ray.org + segment.distMin * ray.dir;
        Point3D localParam;
        m_region.calculateLocalCoordinates(initialPoint, &localParam);
        localParam = clamp(localParam, Point3D(0.0f), Point3D(1.0f));
        int32_t initSuperVoxel[3] = {
            std::min((int32_t)(localParam.x * (m_svNumX - 1)), (int32_t)m_svNumX - 2), 
            std::min((int32_t)(localParam.y * (m_svNumY - 1)), (int32_t)m_svNumY - 2), 
            std::min((int32_t)(localParam.z * (m_svNumZ - 1)), (int32_t)m_svNumZ - 2), 
        };
        int32_t step[3] = {0, 0, 0};
        float delta_t[3];
        int32_t outsideIndices[3] = {-1, -1, -1};
        float initMax_t[3];
        for (int i = 0; i < 3; ++i) {
            delta_t[i] = std::fabs(m_superVoxelWidth[i] / ray.dir[i]);
            if (ray.dir[i] > 0) {
                step[i] = 1;
                outsideIndices[i] = MaxVoxelIndices[i] + 1;
                float nextPlaneCoord = m_region.minP[i] + m_superVoxelWidth[i] * (initSuperVoxel[i] + 1);
                initMax_t[i] = (nextPlaneCoord - ray.org[i]) / ray.dir[i];
            }
            else {
                step[i] = -1;
                outsideIndices[i] = -1;
                float nextPlaneCoord = m_region.minP[i] + m_superVoxelWidth[i] * (initSuperVoxel[i] + 0);
                initMax_t[i] = (nextPlaneCoord - ray.org[i]) / ray.dir[i];
            }
        }

#   if defined(UseRatioTracking)
        // Ratio Tracking
        const auto estimateTransmittance = [&, this](int wl) {
            float transmittance = 1.0f;
            int32_t superVoxel[] = {initSuperVoxel[0], initSuperVoxel[1], initSuperVoxel[2]};
            float max_t[] = {initMax_t[0], initMax_t[1], initMax_t[2]};
            
            float sampledDistance = segment.distMin;
            while (sampledDistance < segment.distMax) {
                float majorantAtScattering;
                bool tentativeHit = traverseSuperVoxels(ray, RaySegment(sampledDistance, segment.distMax), sampler, base_sigma_e[wl], 
                                                        step, delta_t, outsideIndices,  
                                                        max_t, superVoxel, 
                                                        &sampledDistance, &majorantAtScattering);
                if (tentativeHit) {
                    Point3D queryPoint = ray.org + sampledDistance * ray.dir;
                    Point3D param;
                    m_region.calculateLocalCoordinates(queryPoint, &param);
                    float extCoeff = base_sigma_e[wl] * calcDensity(param);
                    float probRealCollision = extCoeff / majorantAtScattering;
                    if (majorantAtScattering != 0.0f)
                        transmittance *= (1.0f - probRealCollision);
                    const float RRThreshold = 0.1f; 
                    if (transmittance < RRThreshold) {
                        if (sampler.getSample() < transmittance)
                            transmittance = 1.0f;
                        else
                            return 0.0f;
                    }
                    SLRAssert(std::isfinite(transmittance), "Invalid Value.\n"
                              "extCoeff: %g, majorantAtScattering: %g", extCoeff, majorantAtScattering);
                }
            }
            
            return transmittance;
        };
#   else
        // Delta Tracking
        const auto estimateTransmittance = [&, this](int wl) {
            int32_t superVoxel[] = {initSuperVoxel[0], initSuperVoxel[1], initSuperVoxel[2]};
            float max_t[] = {initMax_t[0], initMax_t[1], initMax_t[2]};
            
            float sampledDistance = segment.distMin;
            while (sampledDistance < segment.distMax) {
                float majorantAtScattering;
                bool tentativeHit = traverseSuperVoxels(ray, RaySegment(sampledDistance, segment.distMax), sampler, base_sigma_e[wl], 
                                                        step, delta_t, outsideIndices,  
                                                        max_t, superVoxel, 
                                                        &sampledDistance, &majorantAtScattering);
                if (tentativeHit) {
                    Point3D queryPoint = ray.org + sampledDistance * ray.dir;
                    Point3D param;
                    m_region.calculateLocalCoordinates(queryPoint, &param);
                    SampledSpectrum extCoeff = base_sigma_e * calcDensity(param);
                    float probRealCollision = extCoeff[wl] / majorantAtScattering;
                    if (sampler.getSample() < probRealCollision)
                        return 0.0f;
                }
            }
            
            return 1.0f;
        };
#   endif
#else
#   if defined(UseRatioTracking)
        // Ratio Tracking
        const auto estimateTransmittance = [&, this](int wl) {
            float transmittance = 1.0f;
            float majorantWL = majorantExtinctionCoefficientAtWavelength(wls[wl]);
            FloatSum sampledDistance = segment.distMin;
            sampledDistance += -std::log(sampler.getSample()) / majorantWL;
            while (sampledDistance < segment.distMax) {
                Point3D queryPoint = ray.org + sampledDistance * ray.dir;
                Point3D param;
                m_region.calculateLocalCoordinates(queryPoint, &param);
                float density = calcDensity(param);
                float extCoeff = base_sigma_e[wl] * density;
                float probRealCollision = extCoeff / majorantWL;
                transmittance *= (1.0f - probRealCollision);
                
                const float RRThreshold = 0.1f; 
                if (transmittance < RRThreshold) {
                    if (sampler.getSample() < transmittance)
                        transmittance = 1.0f;
                    else
                        return 0.0f;
                }
                
                sampledDistance += -std::log(sampler.getSample()) / majorantWL;
            }
            
            return transmittance;
        };
#   else
        // Delta Tracking
        const auto estimateTransmittance = [&, this](int wl) {
            float majorantWL = majorantExtinctionCoefficientAtWavelength(wls[wl]);
            FloatSum sampledDistance = segment.distMin;
            sampledDistance += -std::log(sampler.getSample()) / majorantWL;
            while (sampledDistance < segment.distMax) {
                Point3D queryPoint = ray.org + sampledDistance * ray.dir;
                Point3D param;
                m_region.calculateLocalCoordinates(queryPoint, &param);
                float density = calcDensity(param);
                float extCoeff = base_sigma_e[wl] * density;
                float probRealCollision = extCoeff / majorantWL;
                if (sampler.getSample() < probRealCollision)
                    return 0.0f;
                sampledDistance += -std::log(sampler.getSample()) / majorantWL;
            }
            
            return 1.0f;
        };
#   endif
#endif
        
        // estimate transmittance by ratio tracking.
        *singleWavelength = false;
        SampledSpectrum transmittance = SampledSpectrum::Zero;
        if (wls.wavelengthSelected()) {
            transmittance[wls.selectedLambdaIndex] = estimateTransmittance(wls.selectedLambdaIndex);
        }
        else {
            for (int wl = 0; wl < WavelengthSamples::NumComponents; ++wl)
                transmittance[wl] = estimateTransmittance(wl);
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
            param.x > 1 || param.y > 1 || param.z > 1)
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

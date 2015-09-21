//
//  Spectrum.cpp
//
//  Created by 渡部 心 on 2015/09/21.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "Spectrum.h"
#include "../BasicTypes/CompensatedSum.h"

static SpectrumFloat integralCMF;

void initSpectrum() {
    DiscretizedSpectrum::init();
    
    CompensatedSum<SpectrumFloat> cum(0);
    for (int i = 1; i < NumCMFSamples; ++i)
        cum += (ybar_2deg[i - 1] + ybar_2deg[i]) * 1 * 0.5;
    integralCMF = cum;
}

#ifdef Use_Spectral_Representation
InputSpectrumRef createInputSpectrum(SpectrumType spType, ColorSpace space, SpectrumFloat e0, SpectrumFloat e1, SpectrumFloat e2) {
    return createShared<UpsampledContinuousSpectrum>(spType, space, e0, e1, e2);
}
InputSpectrumRef createInputSpectrum(SpectrumType spType, SpectrumFloat minLambda, SpectrumFloat maxLambda, const SpectrumFloat* values, uint32_t numSamples) {
    return createShared<RegularContinuousSpectrum>(minLambda, maxLambda, values, numSamples);
}
InputSpectrumRef createInputSpectrum(SpectrumType spType, const SpectrumFloat* lambdas, const SpectrumFloat* values, uint32_t numSamples) {
    return createShared<IrregularContinuousSpectrum>(lambdas, values, numSamples);
}
#else
static void spectrum_to_XYZ(SpectrumFloat minLambda, SpectrumFloat maxLambda, const SpectrumFloat* values, uint32_t numSamples, SpectrumFloat XYZ[3]) {
    const SpectrumFloat CMFBinWidth = (WavelengthHighBound - WavelengthLowBound) / (NumCMFSamples - 1);
    const SpectrumFloat binWidth = (maxLambda - minLambda) / (numSamples - 1);
    uint32_t curCMFIdx = 0;
    uint32_t baseIdx = 0;
    SpectrumFloat curWL = WavelengthLowBound;
    SpectrumFloat prev_xbarVal = 0, prev_ybarVal = 0, prev_zbarVal = 0;
    SpectrumFloat prevValue = 0;
    SpectrumFloat halfWidth = 0;
    CompensatedSum<SpectrumFloat> X(0), Y(0), Z(0);
    while (true) {
        SpectrumFloat xbarValue, ybarValue, zbarValue;
        if (curWL == WavelengthLowBound + curCMFIdx * CMFBinWidth) {
            xbarValue = xbar_2deg[curCMFIdx];
            ybarValue = ybar_2deg[curCMFIdx];
            zbarValue = zbar_2deg[curCMFIdx];
            ++curCMFIdx;
        }
        else {
            uint32_t idx = std::min(uint32_t((curWL - WavelengthLowBound) / CMFBinWidth), NumCMFSamples - 1);
            SpectrumFloat CMFBaseWL = WavelengthLowBound + idx * CMFBinWidth;
            SpectrumFloat t = (curWL - CMFBaseWL) / CMFBinWidth;
            xbarValue = (1 - t) * xbar_2deg[idx] + t * xbar_2deg[idx + 1];
            ybarValue = (1 - t) * ybar_2deg[idx] + t * ybar_2deg[idx + 1];
            zbarValue = (1 - t) * zbar_2deg[idx] + t * zbar_2deg[idx + 1];
        }
        
        SpectrumFloat value;
        if (curWL < minLambda) {
            value = values[0];
        }
        else if (curWL > maxLambda) {
            value = values[numSamples - 1];
        }
        else if (curWL == minLambda + baseIdx * binWidth) {
            value = values[baseIdx];
            ++baseIdx;
        }
        else {
            uint32_t idx = std::min(uint32_t((curWL - minLambda) / binWidth), numSamples - 1);
            SpectrumFloat baseWL = minLambda + idx * binWidth;
            SpectrumFloat t = (curWL - baseWL) / binWidth;
            value = (1 - t) * values[idx] + t * values[idx + 1];
        }
        
        SpectrumFloat avgValue = (prevValue + value) * 0.5f;
        X += avgValue * (prev_xbarVal + xbarValue) * halfWidth;
        Y += avgValue * (prev_ybarVal + ybarValue) * halfWidth;
        Z += avgValue * (prev_zbarVal + zbarValue) * halfWidth;
        
        prev_xbarVal = xbarValue;
        prev_ybarVal = ybarValue;
        prev_zbarVal = zbarValue;
        prevValue = value;
        SpectrumFloat prevWL = curWL;
        curWL = std::min(WavelengthLowBound + curCMFIdx * CMFBinWidth,
                         baseIdx < numSamples ? (minLambda + baseIdx * binWidth) : INFINITY);
        halfWidth = (curWL - prevWL) * 0.5f;
        
        if (curCMFIdx == NumCMFSamples)
            break;
    }
    XYZ[0] = X / integralCMF;
    XYZ[1] = Y / integralCMF;
    XYZ[2] = Z / integralCMF;
}

static void spectrum_to_XYZ(const SpectrumFloat* lambdas, const SpectrumFloat* values, uint32_t numSamples, SpectrumFloat XYZ[3]) {
    const SpectrumFloat CMFBinWidth = (WavelengthHighBound - WavelengthLowBound) / (NumCMFSamples - 1);
    uint32_t curCMFIdx = 0;
    uint32_t baseIdx = 0;
    SpectrumFloat curWL = WavelengthLowBound;
    SpectrumFloat prev_xbarVal = 0, prev_ybarVal = 0, prev_zbarVal = 0;
    SpectrumFloat prevValue = 0;
    SpectrumFloat halfWidth = 0;
    CompensatedSum<SpectrumFloat> X(0), Y(0), Z(0);
    while (true) {
        SpectrumFloat xbarValue, ybarValue, zbarValue;
        if (curWL == WavelengthLowBound + curCMFIdx * CMFBinWidth) {
            xbarValue = xbar_2deg[curCMFIdx];
            ybarValue = ybar_2deg[curCMFIdx];
            zbarValue = zbar_2deg[curCMFIdx];
            ++curCMFIdx;
        }
        else {
            uint32_t idx = std::min(uint32_t((curWL - WavelengthLowBound) / CMFBinWidth), NumCMFSamples - 1);
            SpectrumFloat CMFBaseWL = WavelengthLowBound + idx * CMFBinWidth;
            SpectrumFloat t = (curWL - CMFBaseWL) / CMFBinWidth;
            xbarValue = (1 - t) * xbar_2deg[idx] + t * xbar_2deg[idx + 1];
            ybarValue = (1 - t) * ybar_2deg[idx] + t * ybar_2deg[idx + 1];
            zbarValue = (1 - t) * zbar_2deg[idx] + t * zbar_2deg[idx + 1];
        }
        
        SpectrumFloat value;
        if (curWL < lambdas[0]) {
            value = values[0];
        }
        else if (curWL > lambdas[1]) {
            value = values[numSamples - 1];
        }
        else if (curWL == lambdas[baseIdx]) {
            value = values[baseIdx];
            ++baseIdx;
        }
        else {
            const SpectrumFloat* lb = std::lower_bound(lambdas + std::max(baseIdx - 1, 0u), lambdas + numSamples, curWL);
            uint32_t idx = std::max(int32_t(std::distance(lambdas, lb)) - 1, 0);
            SpectrumFloat t = (curWL - lambdas[idx]) / (lambdas[idx + 1] - lambdas[idx]);
            value = (1 - t) * values[idx] + t * values[idx + 1];
        }
        
        SpectrumFloat avgValue = (prevValue + value) * 0.5f;
        X += avgValue * (prev_xbarVal + xbarValue) * halfWidth;
        Y += avgValue * (prev_ybarVal + ybarValue) * halfWidth;
        Z += avgValue * (prev_zbarVal + zbarValue) * halfWidth;
        
        prev_xbarVal = xbarValue;
        prev_ybarVal = ybarValue;
        prev_zbarVal = zbarValue;
        prevValue = value;
        SpectrumFloat prevWL = curWL;
        curWL = std::min(WavelengthLowBound + curCMFIdx * CMFBinWidth, baseIdx < numSamples ? lambdas[baseIdx] : INFINITY);
        halfWidth = (curWL - prevWL) * 0.5f;
        
        if (curCMFIdx == NumCMFSamples)
            break;
    }
    XYZ[0] = X / integralCMF;
    XYZ[1] = Y / integralCMF;
    XYZ[2] = Z / integralCMF;
}

InputSpectrumRef createInputSpectrum(SpectrumType spType, ColorSpace space, SpectrumFloat e0, SpectrumFloat e1, SpectrumFloat e2) {
    SLRAssert(e0 >= 0.0 && e1 >= 0.0 && e2 >= 0.0, "Values should not be minus.");
    switch (space) {
        case ColorSpace::sRGB:
            return createShared<RGBInputSpectrum>(e0, e1, e2);
        case ColorSpace::sRGB_NonLinear: {
            e0 = sRGB_degamma(e0);
            e1 = sRGB_degamma(e1);
            e2 = sRGB_degamma(e2);
            return createShared<RGBInputSpectrum>(e0, e1, e2);
        }
        case ColorSpace::xyY: {
            SpectrumFloat xyY[3] = {e0, e1, e2};
            SpectrumFloat XYZ[3];
            xyY_to_XYZ(xyY, XYZ);
            e0 = XYZ[0];
            e1 = XYZ[1];
            e2 = XYZ[2];
        }
        case ColorSpace::XYZ: {
            SpectrumFloat XYZ[3] = {e0, e1, e2};
            SpectrumFloat RGB[3];
            switch (spType) {
                case SpectrumType::Reflectance:
                    XYZ_to_sRGB_E(XYZ, RGB);
                    break;
                case SpectrumType::Illuminant:
                    XYZ_to_sRGB(XYZ, RGB);
                    break;
                case SpectrumType::IndexOfRefraction:
                    XYZ_to_sRGB_E(XYZ, RGB);
                    break;
                default:
                    break;
            }
            RGB[0] = RGB[0] < 0.0f ? 0.0f : RGB[0];
            RGB[1] = RGB[1] < 0.0f ? 0.0f : RGB[1];
            RGB[2] = RGB[2] < 0.0f ? 0.0f : RGB[2];
            return createShared<RGBInputSpectrum>(RGB[0], RGB[1], RGB[2]);
        }
        default:
            SLRAssert(false, "Invalid color space.");
            return createShared<RGBInputSpectrum>();
    }
}
InputSpectrumRef createInputSpectrum(SpectrumType spType, SpectrumFloat minLambda, SpectrumFloat maxLambda, const SpectrumFloat* values, uint32_t numSamples) {
    SpectrumFloat XYZ[3];
    spectrum_to_XYZ(minLambda, maxLambda, values, numSamples, XYZ);
    SpectrumFloat RGB[3];
    switch (spType) {
        case SpectrumType::Reflectance:
            XYZ_to_sRGB_E(XYZ, RGB);
            break;
        case SpectrumType::Illuminant:
            XYZ_to_sRGB(XYZ, RGB);
            break;
        case SpectrumType::IndexOfRefraction:
            XYZ_to_sRGB_E(XYZ, RGB);
            break;
        default:
            break;
    }
    RGB[0] = RGB[0] < 0.0f ? 0.0f : RGB[0];
    RGB[1] = RGB[1] < 0.0f ? 0.0f : RGB[1];
    RGB[2] = RGB[2] < 0.0f ? 0.0f : RGB[2];
    return createShared<RGBInputSpectrum>(RGB[0], RGB[1], RGB[2]);
}
InputSpectrumRef createInputSpectrum(SpectrumType spType, const SpectrumFloat* lambdas, const SpectrumFloat* values, uint32_t numSamples) {
    SpectrumFloat XYZ[3];
    spectrum_to_XYZ(lambdas, values, numSamples, XYZ);
    SpectrumFloat RGB[3];
    switch (spType) {
        case SpectrumType::Reflectance:
            XYZ_to_sRGB_E(XYZ, RGB);
            break;
        case SpectrumType::Illuminant:
            XYZ_to_sRGB(XYZ, RGB);
            break;
        case SpectrumType::IndexOfRefraction:
            XYZ_to_sRGB_E(XYZ, RGB);
            break;
        default:
            break;
    }
    RGB[0] = RGB[0] < 0.0f ? 0.0f : RGB[0];
    RGB[1] = RGB[1] < 0.0f ? 0.0f : RGB[1];
    RGB[2] = RGB[2] < 0.0f ? 0.0f : RGB[2];
    return createShared<RGBInputSpectrum>(RGB[0], RGB[1], RGB[2]);
}
#endif

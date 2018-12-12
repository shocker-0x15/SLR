//
//  spectrum_base.h
//
//  Created by 渡部 心 on 2015/05/04.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_spectrum_base__
#define __SLR_spectrum_base__

#include "../defines.h"
#include "../declarations.h"

namespace SLR {
    //             |  RGB Rendering|  Spectral Rendering|   |            typedefs
    // ------------+---------------+--------------------+   +--------------------
    //        Input|            RGB|  ContinuousSpectrum| > |       AssetSpectrum
    //  Calculation|            RGB|     SampledSpectrum| > |     SampledSpectrum
    //             | (+ RGBSamples)| + WavelengthSamples| > |   WavelengthSamples
    // Contribution|            RGB| DiscretizedSpectrum| > | DiscretizedSpectrum
    //      Storage|     RGBStorage|     SpectrumStorage| > |     SpectrumStorage
    
    enum class ColorSpace {
        sRGB,
        sRGB_NonLinear,
        XYZ,
        xyY,
    };
    
    enum class RGBColorSpace {
        sRGB = 0,
    };
    
    enum class SpectrumType : uint32_t {
        Reflectance = 0,
        Illuminant,
        IndexOfRefraction,
    };
    
    template <typename RealType>
    SLR_API RealType sRGB_gamma(RealType value);
    
    template <typename RealType>
    SLR_API RealType sRGB_degamma(RealType value);
    
    // TODO: implement a method to generate arbitrary XYZ<->RGB matrices.
    // http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
    //template <typename RealType>
    //static void RGB_to_XYZ(RealType xR, RealType yR, RealType xG, RealType yG, RealType xB, RealType yB, RealType xW, RealType yW) {
    //    RealType XR = xR / yR, YR = 1, ZR = (1 - xR - yR) / yR;
    //    RealType XG = xG / yG, YG = 1, ZG = (1 - xG - yG) / yG;
    //    RealType XB = xB / yB, YB = 1, ZB = (1 - xB - yB) / yB;
    //}
    
    template <typename RealType>
    inline float sRGB_to_Luminance(const RealType r, const RealType g, const RealType b) {
        return 0.2126729 * r + 0.7151522 * g + 0.0721750 * b;
    }
    
    template <typename RealType>
    inline void sRGB_to_XYZ(const RealType rgb[3], RealType xyz[3]) {
        xyz[0] = 0.4124564 * rgb[0] + 0.3575761 * rgb[1] + 0.1804375 * rgb[2];
        xyz[1] = 0.2126729 * rgb[0] + 0.7151522 * rgb[1] + 0.0721750 * rgb[2];
        xyz[2] = 0.0193339 * rgb[0] + 0.1191920 * rgb[1] + 0.9503041 * rgb[2];
    }
    
    template <typename RealType>
    inline void XYZ_to_sRGB(const RealType xyz[3], RealType rgb[3]) {
        rgb[0] = 3.2404542 * xyz[0] - 1.5371385 * xyz[1] - 0.4985314 * xyz[2];
        rgb[1] = -0.9692660 * xyz[0] + 1.8760108 * xyz[1] + 0.0415560 * xyz[2];
        rgb[2] = 0.0556434 * xyz[0] - 0.2040259 * xyz[1] + 1.0572252 * xyz[2];
    }
    
    template <typename RealType>
    inline void sRGB_E_to_XYZ(const RealType rgb[3], RealType xyz[3]) {
        xyz[0] = 0.4969 * rgb[0] + 0.3391 * rgb[1] + 0.1640 * rgb[2];
        xyz[1] = 0.2562 * rgb[0] + 0.6782 * rgb[1] + 0.0656 * rgb[2];
        xyz[2] = 0.0233 * rgb[0] + 0.1130 * rgb[1] + 0.8637 * rgb[2];
    }
    
    template <typename RealType>
    inline void XYZ_to_sRGB_E(const RealType xyz[3], RealType rgb[3]) {
        rgb[0] = 2.6897 * xyz[0] - 1.2759 * xyz[1] - 0.4138 * xyz[2];
        rgb[1] = -1.0221 * xyz[0] + 1.9783 * xyz[1] + 0.0438 * xyz[2];
        rgb[2] = 0.0612 * xyz[0] - 0.2245 * xyz[1] + 1.1633 * xyz[2];
    }
    
    template <typename RealType>
    inline void XYZ_to_xyY(const RealType xyz[3], RealType xyY[3]) {
        RealType b = xyz[0] + xyz[1] + xyz[2];
        if (b == 0) {
            xyY[0] = xyY[1] = 1.0f / 3.0f;
            xyY[2] = 0.0f;
            return;
        }
        xyY[0] = xyz[0] / b;
        xyY[1] = xyz[1] / b;
        xyY[2] = xyz[1];
    }
    
    template <typename RealType>
    inline void xyY_to_XYZ(const RealType xyY[3], RealType xyz[3]) {
        RealType b = xyY[2] / xyY[1];
        xyz[0] = xyY[0] * b;
        xyz[1] = xyY[2];
        xyz[2] = (1.0f - xyY[0] - xyY[1]) * b;
    }
    
#if SLR_Color_System_is_based_on == SLR_Color_System_CIE_1931_2deg
    static const SpectrumFloat WavelengthLowBound = 360.0f;
    static const SpectrumFloat WavelengthHighBound = 830.0f;
    static const uint32_t NumCMFSamples = 471;
    
    extern const SpectrumFloat xbar_CIE1931_2deg[NumCMFSamples];
    extern const SpectrumFloat ybar_CIE1931_2deg[NumCMFSamples];
    extern const SpectrumFloat zbar_CIE1931_2deg[NumCMFSamples];
#elif SLR_Color_System_is_based_on == SLR_Color_System_CIE_1964_10deg

#elif SLR_Color_System_is_based_on == SLR_Color_System_CIE_2012_2deg

#elif SLR_Color_System_is_based_on == SLR_Color_System_CIE_2012_10deg
    
#endif

    
    
    extern SLR_API SpectrumFloat integralCMF;
    extern SLR_API std::unique_ptr<ContinuousSpectrum> ybarSpectrum;
    
    // JP: 色システムの初期化を行う。このライブラリを使う前に呼ぶ必要がある。
    // EN: initialize the color system. It is required to call this before using this library.
    SLR_API void initializeColorSystem();
}

#endif /* __SLR_spectrum_base__ */

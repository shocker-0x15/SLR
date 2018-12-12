//
//  directional_distribution_functions.h
//
//  Created by 渡部 心 on 2015/06/28.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_directional_distribution_functions__
#define __SLR_directional_distribution_functions__

#include "../defines.h"
#include "../declarations.h"
#include "../BasicTypes/spectrum_base.h"
#include "geometry.h"
#include <type_traits>

namespace SLR {
    struct SLR_API DirectionType {
        enum InternalEnum : uint32_t {
            IE_LowFreq = 1 << 0,
            IE_HighFreq = 1 << 1,
            IE_Delta0D = 1 << 2,
            IE_Delta1D = 1 << 3,
            IE_NonDelta = IE_LowFreq | IE_HighFreq,
            IE_Delta = IE_Delta0D | IE_Delta1D,
            IE_AllFreq = IE_NonDelta | IE_Delta,
            
            IE_Reflection = 1 << 4,
            IE_Transmission = 1 << 5,
            IE_Emission = IE_Reflection,
            IE_Acquisition = IE_Reflection,
            IE_WholeSphere = IE_Reflection | IE_Transmission,
            
            IE_All = IE_AllFreq | IE_WholeSphere,
            
            IE_Dispersive = 1 << 6,
            
            IE_LowFreqReflection = IE_LowFreq | IE_Reflection,
            IE_LowFreqTransmission = IE_LowFreq | IE_Transmission,
            IE_LowFreqScattering = IE_LowFreqReflection | IE_LowFreqTransmission,
            IE_HighFreqReflection = IE_HighFreq | IE_Reflection,
            IE_HighFreqTransmission = IE_HighFreq | IE_Transmission,
            IE_HighFreqScattering = IE_HighFreqReflection | IE_HighFreqTransmission,
            IE_Delta0DReflection = IE_Delta0D | IE_Reflection,
            IE_Delta0DTransmission = IE_Delta0D | IE_Transmission,
            IE_Delta0DScattering = IE_Delta0DReflection | IE_Delta0DTransmission,
        };
        static const DirectionType LowFreq;
        static const DirectionType HighFreq;
        static const DirectionType Delta0D;
        static const DirectionType Delta1D;
        static const DirectionType NonDelta;
        static const DirectionType Delta;
        static const DirectionType AllFreq;
        static const DirectionType Reflection;
        static const DirectionType Transmission;
        static const DirectionType Emission;
        static const DirectionType Acquisition;
        static const DirectionType WholeSphere;
        static const DirectionType All;
        static const DirectionType Dispersive;
        static const DirectionType LowFreqReflection;
        static const DirectionType LowFreqTransmission;
        static const DirectionType LowFreqScattering;
        static const DirectionType HighFreqReflection;
        static const DirectionType HighFreqTransmission;
        static const DirectionType HighFreqScattering;
        static const DirectionType Delta0DReflection;
        static const DirectionType Delta0DTransmission;
        static const DirectionType Delta0DScattering;
        
        InternalEnum value;
        
        DirectionType(InternalEnum v = (InternalEnum)0) : value(v) { }
        DirectionType operator&(const DirectionType &r) const { return (InternalEnum)(value & r.value); }
        DirectionType operator|(const DirectionType &r) const { return (InternalEnum)(value | r.value); }
        DirectionType &operator&=(const DirectionType &r) { value = (InternalEnum)(value & r.value); return *this; }
        DirectionType &operator|=(const DirectionType &r) { value = (InternalEnum)(value | r.value); return *this; }
        DirectionType flip() const { return InternalEnum(value ^ IE_WholeSphere); }
        operator bool() const { return value; }
        bool operator==(const DirectionType &r) const { return value == r.value; }
        bool operator!=(const DirectionType &r) const { return value != r.value; }
        
        bool matches(DirectionType t) const { uint32_t res = value & t.value; return (res & IE_WholeSphere) && (res & IE_AllFreq); }
        bool hasNonDelta() const { return value & IE_NonDelta; }
        bool hasDelta() const { return value & IE_Delta; }
        bool isDelta() const { return (value & IE_Delta) && !(value & IE_NonDelta); }
        bool isReflection() const { return (value & IE_Reflection) && !(value & IE_Transmission); }
        bool isTransmission() const { return !(value & IE_Reflection) && (value & IE_Transmission); }
        bool isDispersive() const { return value & IE_Dispersive; }
    };
    
    
    
    struct SLR_API EDFQuery {
        DirectionType flags;
        EDFQuery(DirectionType f = DirectionType::All) : flags(f) { };
    };
    
    struct SLR_API EDFSample {
        float uComponent;
        float uDir[2];
        EDFSample(float uComp, float uDir0, float uDir1) : uComponent(uComp), uDir{uDir0, uDir1} {}
    };
    
    struct SLR_API EDFQueryResult {
        Vector3D dir_sn;
        float dirPDF;
        DirectionType dirType;
    };
    
    
    
    struct SLR_API ABDFQuery {
        Vector3D dirLocal;
        int16_t wlHint;
        DirectionType dirTypeFilter;
        bool requestReverse;
        
        ABDFQuery(const Vector3D &dirL, int16_t wl, DirectionType filter, bool reqRev) : dirLocal(dirL), wlHint(wl), dirTypeFilter(filter), requestReverse(reqRev) { }
    };
    
    struct SLR_API ABDFSample {
        float uComponent;
        float uDir[2];
        
        ABDFSample() { }
        ABDFSample(float uComp, float uDir0, float uDir1) : uComponent(uComp), uDir{uDir0, uDir1} { }
    };
    
    struct SLR_API ABDFReverseInfo {
        SampledSpectrum value;
        float dirPDF;
    };
    
    struct SLR_API ABDFQueryResult {
        Vector3D dirLocal;
        float dirPDF;
        DirectionType sampledType;
        ABDFReverseInfo reverse;
        ABDFQueryResult() { }
        
        SampledSpectrum reverseValue() const {
            return reverse.value;
        }
        
        float reverseDirPDF() const {
            return reverse.dirPDF;
        }
    };
    
    
    
    struct SLR_API BSDFQuery : public ABDFQuery {
        Normal3D gNormalLocal;
        bool adjoint;
        
        BSDFQuery(const Vector3D &dirSN, const Normal3D &gNormalSN, int16_t wl, DirectionType filter = DirectionType::All, bool reqRev = false, bool adj = false) :
        ABDFQuery(dirSN, wl, filter, reqRev), gNormalLocal(gNormalSN), adjoint(adj) { }
    };
    
    struct SLR_API BSDFSample : public ABDFSample {
        BSDFSample() : ABDFSample() {}
        BSDFSample(float uComp, float uDir0, float uDir1) : ABDFSample(uComp, uDir0, uDir1) { }
    };
    
    struct SLR_API BSDFReverseInfo : public ABDFReverseInfo {
    };
    
    struct SLR_API BSDFQueryResult : public ABDFQueryResult {
    };
    
    
    
    struct SLR_API PFQuery {
        Vector3D dirLocal;
        int16_t wlHint;
        DirectionType dirTypeFilter;
        bool requestReverse;
        
        PFQuery(const Vector3D &dirL, int16_t wl, DirectionType filter = DirectionType::All, bool reqRev = false) :
        dirLocal(dirL), wlHint(wl), dirTypeFilter(filter), requestReverse(reqRev) { }
    };
    
    struct SLR_API PFSample {
        float uDir[2];
        PFSample(float uDir0, float uDir1) : uDir{uDir0, uDir1} { }
    };
    
    struct SLR_API PFQueryResult {
        Vector3D dirLocal;
        float dirPDF;
        DirectionType sampledType;
        float revDirPDF;
    };
    
    
    
    struct SLR_API VolumetricBSDFQuery : public ABDFQuery {
        VolumetricBSDFQuery(const Vector3D &dirSN, int16_t wl, DirectionType filter, bool reqRev) :
        ABDFQuery(dirSN, wl, filter, reqRev) { }
    };
    
    struct SLR_API VolumetricBSDFSample : public ABDFSample {
        VolumetricBSDFSample() : ABDFSample() {}
        VolumetricBSDFSample(float uComp, float uDir0, float uDir1) : ABDFSample(uComp, uDir0, uDir1) { }
    };
    
    struct SLR_API VolumetricBSDFQueryResult : public ABDFQueryResult {
        VolumetricBSDFQueryResult() { }
        VolumetricBSDFQueryResult(const PFQueryResult &pfResult) {
            dirLocal = pfResult.dirLocal;
            dirPDF = pfResult.dirPDF;
            sampledType = pfResult.sampledType;
        }
    };
    

    
    struct SLR_API IDFSample {
        float uDir[2];
        IDFSample(float uDir0, float uDir1) : uDir{uDir0, uDir1} { }
    };
    
    struct SLR_API IDFQueryResult {
        Vector3D dirLocal;
        float dirPDF;
        DirectionType dirType;
    };
    
    
    
    class SLR_API EDF {
    protected:
        const DirectionType m_type;
        
        bool matches(DirectionType flags, EDFQueryResult* result) const {
            if (matches(flags))
                return true;
            result->dirPDF = 0.0f;
            result->dirType = DirectionType();
            return false;
        }
        friend class MultiEDF;
    public:
        EDF(DirectionType type) : m_type(type) { };
        virtual ~EDF() { };
        
        virtual SampledSpectrum sample(const EDFQuery &query, const EDFSample &smp, EDFQueryResult* result) const = 0;
        virtual SampledSpectrum evaluate(const EDFQuery &query, const Vector3D &dir) const = 0;
        virtual float evaluatePDF(const EDFQuery &query, const Vector3D &dir) const = 0;
        virtual float weight(const EDFQuery &query) const = 0;
        
        virtual bool matches(DirectionType flags) const { return m_type.matches(flags); }
        bool hasNonDelta() const { return matches(DirectionType::WholeSphere | DirectionType::NonDelta); }
    };
    
    
    
    // Abstract Bidirectional Distribution Function
    class SLR_API AbstractBDF {
    protected:
        DirectionType m_type;
    public:
        AbstractBDF(DirectionType type) : m_type(type) { }
        virtual ~AbstractBDF() { }
        
        virtual SampledSpectrum sample(const ABDFQuery* query, LightPathSampler &sampler, ArenaAllocator &mem, ABDFQueryResult** result) const = 0;
        virtual SampledSpectrum evaluate(const ABDFQuery* query, const Vector3D &dir, SampledSpectrum* rev_fs = nullptr) const = 0;
        virtual float evaluatePDF(const ABDFQuery* query, const Vector3D &dir, float* revPDF = nullptr) const = 0;
        
        virtual bool matches(DirectionType flags) const { return m_type.matches(flags); };
        bool hasDelta() const { return matches(DirectionType::WholeSphere | DirectionType::Delta); }
        bool hasNonDelta() const { return matches(DirectionType::WholeSphere | DirectionType::NonDelta); }
    };
    
    
    
    class SLR_API BSDF : public AbstractBDF {
    protected:        
        DirectionType sideTest(const Normal3D& ng, const Vector3D& d0, const Vector3D &d1) const {
            bool reflect = dot(Vector3D(ng), d0) * dot(Vector3D(ng), d1) > 0;
            return DirectionType::AllFreq | (reflect ? DirectionType::Reflection : DirectionType::Transmission);
        }
        virtual SampledSpectrum sampleInternal(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult* result) const = 0;
        virtual SampledSpectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const = 0;
        virtual float evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, float* revPDF) const = 0;
        virtual float weightInternal(const BSDFQuery &query) const = 0;
        virtual SampledSpectrum getBaseColorInternal(DirectionType flags) const = 0;
        friend class MultiBSDF;
        friend class FlippedBSDF;
    public:
        BSDF(DirectionType type) : SLR::AbstractBDF(type) { }
        
        SampledSpectrum sample(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const {
            if (!matches(query.dirTypeFilter)) {
                result->dirPDF = 0.0f;
                result->sampledType = DirectionType();
                return SampledSpectrum::Zero;
            }
            SampledSpectrum fs_sn = sampleInternal(query, smp.uComponent, smp.uDir, result);
            float snCorrection = (query.adjoint ?
                                  std::fabs(query.dirLocal.z / dot(query.dirLocal, query.gNormalLocal)) :
                                  std::fabs(result->dirLocal.z / dot(result->dirLocal, query.gNormalLocal)));
            //float snCorrection = std::fabs(result->dirLocal.z / dot(result->dirLocal, query.gNormalLocal));
            if (query.requestReverse)
                result->reverse.value *= snCorrection;
            SLRAssert(result->dirPDF == 0 || (fs_sn.allFinite() && !fs_sn.hasNegative() && std::isfinite(snCorrection) &&
                                              std::isfinite(result->dirPDF)),
                      "fs_sn: %s, snCorrection: %g, PDF: %g, wlIdx: %u, qDir: %s, sample: (%g, %g, %g), rDir: %s, gNormal: %s",
                      fs_sn.toString().c_str(), snCorrection, result->dirPDF, query.wlHint,
                      query.dirLocal.toString().c_str(), smp.uComponent, smp.uDir[0], smp.uDir[1], result->dirLocal.toString().c_str(), query.gNormalLocal.toString().c_str());
            return fs_sn * snCorrection;
        }
        SampledSpectrum evaluate(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs = nullptr) const {
//            BSDFQuery mQuery = query;
//            mQuery.dirTypeFilter &= sideTest(query.gNormalLocal, query.dirLocal, dir);
//            if (!matches(mQuery.dirTypeFilter)) {
//                if (query.requestReverse)
//                    *rev_fs = SampledSpectrum::Zero;
//                return SampledSpectrum::Zero;
//            }
//            SampledSpectrum fs_sn = evaluateInternal(mQuery, dir, rev_fs);
            SampledSpectrum fs_sn = evaluateInternal(query, dir, rev_fs);
            float snCorrection = (query.adjoint ?
                                  std::fabs(query.dirLocal.z / dot(query.dirLocal, query.gNormalLocal)) :
                                  std::fabs(dir.z / dot(dir, query.gNormalLocal)));
            //float snCorrection = std::fabs(dir.z / dot(dir, query.gNormalLocal));
            SLRAssert(fs_sn.allFinite() && !fs_sn.hasNegative() && std::isfinite(snCorrection),
                      "fs_sn: %s, snCorrection: %g, wlIdx: %u, qDir: %s, rDir: %s, gNormal: %s",
                      fs_sn.toString().c_str(), snCorrection, query.wlHint,
                      query.dirLocal.toString().c_str(), dir.toString().c_str(), query.gNormalLocal.toString().c_str());
            
            if (query.requestReverse)
                *rev_fs *= snCorrection;
            return fs_sn * snCorrection;
        }
        float evaluatePDF(const BSDFQuery &query, const Vector3D &dir, float* revPDF = nullptr) const {
            if (!matches(query.dirTypeFilter)) {
                if (query.requestReverse)
                    *revPDF = 0;
                return 0;
            }
            float ret = evaluatePDFInternal(query, dir, revPDF);
            SLRAssert(std::isfinite(ret) && ret >= 0,
                      "PDF: %g, wlIdx: %u, qDir: %s, rDir: %s",
                      ret, query.wlHint, query.dirLocal.toString().c_str(), dir.toString().c_str());
            return ret;
        }
        float weight(const BSDFQuery &query) const {
            if (!matches(query.dirTypeFilter))
                return 0;
            float weight_sn = weightInternal(query);
            float snCorrection = query.adjoint ? std::fabs(query.dirLocal.z / dot(query.dirLocal, query.gNormalLocal)) : 1;
            SLRAssert(std::isfinite(weight_sn) && weight_sn >= 0.0f && std::isfinite(snCorrection),
                      "weight_sn: %g, snCorrection: %g, wlIdx: %u, qDir: %s",
                      weight_sn, snCorrection, query.wlHint, query.dirLocal.toString().c_str());
            return weight_sn * snCorrection;
        }
        
        SampledSpectrum rho(uint32_t numSamples, BSDFSample* samples, float* uDir0, float* uDir1, float* uWl, DirectionType flags = DirectionType::All, bool fromUpper = true) const;
        
        SampledSpectrum getBaseColor(DirectionType flags) const {
            if (!matches(flags))
                return SampledSpectrum::Zero;
            return getBaseColorInternal(flags);
        }
        
        SampledSpectrum sample(const ABDFQuery* query, LightPathSampler &sampler, ArenaAllocator &mem, ABDFQueryResult** result) const override;
        SampledSpectrum evaluate(const ABDFQuery* query, const Vector3D &dir, SampledSpectrum* rev_fs = nullptr) const override;
        float evaluatePDF(const ABDFQuery* query, const Vector3D &dir, float* revPDF = nullptr) const override;
    };
    
    
    
    class SLR_API PhaseFunction {
    public:
        PhaseFunction() {}
        virtual ~PhaseFunction() {}
        
        virtual DirectionType directionType() const = 0;
        
        virtual SampledSpectrum sample(const PFQuery &query, const PFSample &smp, PFQueryResult* result) const = 0;
        virtual SampledSpectrum evaluate(const PFQuery &query, const Vector3D &dirIn) const = 0;
        virtual float evaluatePDF(const PFQuery &query, const Vector3D &dirIn, float* revPDF = nullptr) const = 0;
    };
    
    
    
    class SLR_API VolumetricBSDF : public AbstractBDF {
        const PhaseFunction* m_pf;
        SampledSpectrum m_albedo;
    public:
        VolumetricBSDF(const SampledSpectrum &albedo, const PhaseFunction* pf) : AbstractBDF(pf->directionType()), m_pf(pf), m_albedo(albedo) { }
        
        SampledSpectrum sample(const ABDFQuery* query, LightPathSampler &sampler, ArenaAllocator &mem, ABDFQueryResult** result) const override;
        SampledSpectrum evaluate(const ABDFQuery* query, const Vector3D &dir, SampledSpectrum* rev_fs = nullptr) const override;
        float evaluatePDF(const ABDFQuery* query, const Vector3D &dir, float* revPDF = nullptr) const override;
    };
    
    
    
    class SLR_API IDF {
    protected:
        const DirectionType m_type;
        
    public:
        IDF(DirectionType type) : m_type(type) { }
        virtual ~IDF() { }
        
        virtual SampledSpectrum sample(const IDFSample &smp, IDFQueryResult* result) const = 0;
        virtual SampledSpectrum evaluate(const Vector3D &dirIn) const = 0;
        virtual float evaluatePDF(const Vector3D &dirIn) const = 0;
        virtual void calculatePixel(const Vector3D &dirIn, float* hitPx, float* hitPy) const = 0;
        
        virtual bool matches(DirectionType flags) const { return m_type.matches(flags); }
        bool hasNonDelta() const { return matches(DirectionType::WholeSphere | DirectionType::NonDelta); }
    };
    
    
    
    class SLR_API Fresnel {
    public:
        virtual ~Fresnel() { };
        virtual SampledSpectrum evaluate(float cosEnter) const = 0;
        virtual float evaluate(float cosEnter, uint32_t wlIdx) const = 0;
    };
    
    class SLR_API FresnelNoOp : public Fresnel {
    public:
        SampledSpectrum evaluate(float cosEnter) const override;
        float evaluate(float cosEnter, uint32_t wlIdx) const override;
    };
    
    class SLR_API FresnelConductor : public Fresnel {
        SampledSpectrum m_eta;
        SampledSpectrum m_k;
    public:
        FresnelConductor(const SampledSpectrum &eta, const SampledSpectrum &k) : m_eta(eta), m_k(k) { }
        
        SampledSpectrum evaluate(float cosEnter) const override;
        float evaluate(float cosEnter, uint32_t wlIdx) const override;
    };
    
    class SLR_API FresnelDielectric : public Fresnel {
        SampledSpectrum m_etaExt;
        SampledSpectrum m_etaInt;
    public:
        FresnelDielectric(const SampledSpectrum &etaExt, const SampledSpectrum &etaInt) : m_etaExt(etaExt), m_etaInt(etaInt) { }
        
        SampledSpectrum etaExt() const { return m_etaExt; }
        SampledSpectrum etaInt() const { return m_etaInt; }
        
        SampledSpectrum evaluate(float cosEnter) const override;
        float evaluate(float cosEnter, uint32_t wlIdx) const override;
        
        static float evalF(float etaEnter, float etaExit, float cosEnter, float cosExit);
    };
    
    
    
    class SLR_API MicrofacetDistribution {
    public:
        virtual float sample(float u0, float u1, Normal3D* m, float* normalPDF) const = 0;
        virtual float evaluate(const Normal3D &m) const = 0;
        virtual float evaluatePDF(const Normal3D &m) const = 0;
        
        virtual float sample(const Vector3D &v, float u0, float u1, Normal3D* m, float* normalPDF) const {
            return sample(u0, u1, m, normalPDF);
        }
        virtual float evaluatePDF(const Vector3D &v, const Normal3D &m) const {
            return evaluatePDF(m);
        }
        
        virtual float evaluateSmithG1(const Vector3D &v, const Normal3D &m) const = 0;
    };
    
    class SLR_API BerryMicrofacetDistribution : public MicrofacetDistribution {
        float m_alpha_g;
    public:
        BerryMicrofacetDistribution(float alpha_g) : m_alpha_g(alpha_g) { }
        
        float sample(float u0, float u1, Normal3D* m, float* normalPDF) const override;
        float evaluate(const Normal3D &m) const override;
        float evaluatePDF(const Normal3D &m) const override;
        
        // JP: これらはクラス変数が派生クラスとして宣言される場合に必要となる。
        // EN: These are needed for the case where a class variable is declared as a derived class.
        float sample(const Vector3D &v, float u0, float u1, Normal3D* m, float* normalPDF) const override {
            return MicrofacetDistribution::sample(v, u0, u1, m, normalPDF);
        }
        float evaluatePDF(const Vector3D &v, const Normal3D &m) const override {
            return MicrofacetDistribution::evaluatePDF(v, m);
        }
        
        float evaluateSmithG1(const Vector3D &v, const Normal3D &m) const override;
    };
    
    // GGX: Ground Glass X (a.k.a. Trowbridge-Reitz)
    class SLR_API GGXMicrofacetDistribution : public MicrofacetDistribution {
        float m_alpha_gx, m_alpha_gy;
    public:
        GGXMicrofacetDistribution(float alpha_gx, float alpha_gy) : m_alpha_gx(alpha_gx), m_alpha_gy(alpha_gy) { }
        
        float sample(float u0, float u1, Normal3D* m, float* normalPDF) const override;
        float evaluate(const Normal3D &m) const override;
        float evaluatePDF(const Normal3D &m) const override;
        
        float sample(const Vector3D &v, float u0, float u1, Normal3D* m, float* normalPDF) const override;
        float evaluatePDF(const Vector3D &v, const Normal3D &m) const override;
        
        float evaluateSmithG1(const Vector3D &v, const Normal3D &m) const override;
    };
    
    // GTR: Generalized Trowbridge-Reitz
    // gamma 1 is equivalent to the Berry distribution.
    // gamma 2 is equivalent to the GGX distribution. 
    class SLR_API GTRMicrofacetDistribution : public MicrofacetDistribution {
        float m_gamma;
        float m_alpha_g;
    public:
        GTRMicrofacetDistribution(float gamma, float alpha_g) : m_gamma(gamma), m_alpha_g(alpha_g) { }
        
        float sample(float u0, float u1, Normal3D* m, float* normalPDF) const override;
        float evaluate(const Normal3D &m) const override;
        float evaluatePDF(const Normal3D &m) const override;
        
        // JP: これらはクラス変数が派生クラスとして宣言される場合に必要となる。
        // EN: These are needed for the case where a class variable is declared as a derived class.
        float sample(const Vector3D &v, float u0, float u1, Normal3D* m, float* normalPDF) const override {
            return MicrofacetDistribution::sample(v, u0, u1, m, normalPDF);
        }
        float evaluatePDF(const Vector3D &v, const Normal3D &m) const override {
            return MicrofacetDistribution::evaluatePDF(v, m);
        }
        
        float evaluateSmithG1(const Vector3D &v, const Normal3D &m) const override;
    };
}

#endif /* __SLR_directional_distribution_functions__ */

//
//  directional_distribution_functions.h
//
//  Created by 渡部 心 on 2015/06/28.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__directional_distribution_functions__
#define __SLR__directional_distribution_functions__

#include "../defines.h"
#include "../references.h"
#include "../BasicTypes/Spectrum.h"
#include "geometry.h"
#include <type_traits>

#define BSDF_SAMPLE_ASSERT \
SLRAssert(ret.hasInf() == false && ret.hasNaN() == false,\
"BSDF value is an expected value.\n"\
"uComponent: %f, uDir: (%f, %f), dirIn: (%f, %f, %f), dirOut: (%f, %f, %f)",\
smp.uComponent, smp.uDir[0], smp.uDir[1], \
query.dir_sn.x, query.dir_sn.y, query.dir_sn.z, \
result->dir_sn.x, result->dir_sn.y, result->dir_sn.z)

#define BSDF_EVALUATE_ASSERT \
SLRAssert(ret.hasInf() == false && ret.hasNaN() == false,\
"BSDF value is an expected value.\n"\
"dirIn: (%f, %f, %f), dirOut: (%f, %f, %f)",\
query.dir_sn.x, query.dir_sn.y, query.dir_sn.z, \
dir.x, dir.y, dir.z)

#define BSDF_EVALUATE_PDF_ASSERT \
SLRAssert(std::isinf(ret) == false && std::isnan(ret) == false,\
"PDF value is an expected value. : dirIn: (%f, %f, %f)",\
query.dir_sn.x, query.dir_sn.y, query.dir_sn.z);

template <typename EnumType>
inline EnumType operator|(EnumType e0, EnumType e1) {
    typedef typename std::underlying_type<EnumType>::type T;
    return (EnumType)(static_cast<T>(e0) | static_cast<T>(e1));
}
template <typename EnumType>
inline EnumType operator&(EnumType e0, EnumType e1) {
    typedef typename std::underlying_type<EnumType>::type T;
    return (EnumType)(static_cast<T>(e0) & static_cast<T>(e1));
}

struct DirectionType {
    enum InternalEnum : uint32_t {
        LowFreq = 1 << 0,
        HighFreq = 1 << 1,
        Delta0D = 1 << 2,
        Delta1D = 1 << 3,
        NonDelta = LowFreq | HighFreq,
        Delta = Delta0D | Delta1D,
        AllFreq = NonDelta | Delta,
        
        Reflection = 1 << 4,
        Transmission = 1 << 5,
        WholeSphere = Reflection | Transmission,
        
        All = AllFreq | WholeSphere,
        
        LowFreqReflection = LowFreq | Reflection,
        LowFreqTransmission = LowFreq | Transmission,
        LowFreqScattering = LowFreqReflection | LowFreqTransmission,
        HighFreqReflection = HighFreq | Reflection,
        HighFreqTransmission = HighFreq | Transmission,
        HighFreqScattering = HighFreqReflection | HighFreqTransmission,
        Delta0DReflection = Delta0D | Reflection,
        Delta0DTransmission = Delta0D | Transmission,
        Delta0DScattering = Delta0DReflection | Delta0DTransmission,
    };
    
    InternalEnum value;

    DirectionType(InternalEnum v = (InternalEnum)0) : value(v) { };
    DirectionType operator&(const DirectionType &r) const { return (InternalEnum)(value & r.value); };
    DirectionType operator|(const DirectionType &r) const { return (InternalEnum)(value | r.value); };
    DirectionType &operator&=(const DirectionType &r) { value = (InternalEnum)(value & r.value); return *this; };
    DirectionType &operator|=(const DirectionType &r) { value = (InternalEnum)(value | r.value); return *this; };
    
    bool matches(DirectionType t) const { return (value & t.value) == t.value; };
    bool isDelta() const { return (value & (Delta | WholeSphere)) == value; };
    bool hasNonDelta() const { return value & NonDelta; };
    bool hasDelta() const { return value & Delta; };
};


struct EDFQuery {
    DirectionType flags;
    EDFQuery(DirectionType f = DirectionType::All) : flags(f) { };
};

struct EDFSample {
    float uComponent;
    float uDir[2];
    EDFSample(float uComp, float uDir0, float uDir1) : uComponent(uComp), uDir{uDir0, uDir1} {};
};

struct EDFQueryResult {
    Vector3D dir_sn;
    float dirPDF;
    DirectionType dirType;
};

struct BSDFQuery {
    Vector3D dir_sn;
    Normal3D gNormal_sn;
    DirectionType flags;
    bool adjoint;
    
    BSDFQuery(const Vector3D &dirSN, const Normal3D &gNormalSN, DirectionType f = DirectionType::All, bool adj = false) :
    dir_sn(dirSN), gNormal_sn(gNormalSN), flags(f), adjoint(adj) { };
};

struct BSDFSample {
    float uComponent;
    float uDir[2];
    BSDFSample(float uComp, float uDir0, float uDir1) : uComponent(uComp), uDir{uDir0, uDir1} { }
};

struct BSDFQueryResult {
    Vector3D dir_sn;
    float dirPDF;
    DirectionType dirType;
};

struct IDFSample {
    float uDir[2];
    IDFSample(float uDir0, float uDir1) : uDir{uDir0, uDir1} { };
};

struct IDFQueryResult {
    Vector3D dirLocal;
    float dirPDF;
};


class EDF {
protected:
    DirectionType m_type;
    
    bool matches(DirectionType flags, EDFQueryResult* result) const {
        if (flags.matches(m_type))
            return true;
        result->dirPDF = 0.0f;
        result->dirType = DirectionType();
        return false;
    };
    friend class MultiEDF;
public:
    EDF(DirectionType t) : m_type(t) { };
    virtual ~EDF() { };
    
    virtual Spectrum sample(const EDFQuery &query, const EDFSample &smp, EDFQueryResult* result) const = 0;
    virtual Spectrum evaluate(const EDFQuery &query, const Vector3D &dir) const = 0;
    virtual float evaluatePDF(const EDFQuery &query, const Vector3D &dir) const = 0;
    virtual float weight(const EDFQuery &query) const = 0;
};

class BSDF {
protected:
    DirectionType m_type;
    
    DirectionType sideTest(const Normal3D& ng, const Vector3D& d0, const Vector3D &d1) const {
        bool reflect = dot(Vector3D(ng), d0) * dot(Vector3D(ng), d1) > 0;
        return DirectionType::AllFreq | (reflect ? DirectionType::Reflection : DirectionType::Transmission);
    };
    bool matches(DirectionType flags, BSDFQueryResult* result) const {
        if (flags.matches(m_type))
            return true;
        result->dirPDF = 0.0f;
        result->dirType = DirectionType();
        return false;
    };
    virtual Spectrum evaluateInternal(const BSDFQuery &query, const Vector3D &dir) const = 0;
    friend class MultiBSDF;
public:
    BSDF(DirectionType t) : m_type(t) { };
    virtual ~BSDF() { };
    
    virtual Spectrum sample(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult* result) const = 0;
    Spectrum evaluate(const BSDFQuery &query, const Vector3D &dir) const {
        BSDFQuery mQuery = query;
        mQuery.flags &= sideTest(query.gNormal_sn, query.dir_sn, dir);
        return evaluateInternal(mQuery, dir);
    };
    virtual float evaluatePDF(const BSDFQuery &query, const Vector3D &dir) const = 0;
    virtual float weight(const BSDFQuery &query) const = 0;
    
    virtual Spectrum getBaseColor(DirectionType flags) const = 0;

    virtual bool matches(DirectionType flags) const { return flags.matches(m_type); };
    bool hasNonDelta() const { return m_type.hasNonDelta(); };
};

class IDF {
protected:
    DirectionType m_type;
    
    bool matches(DirectionType flags, BSDFQueryResult* result) const {
        if (flags.matches(m_type))
            return true;
        result->dirPDF = 0.0f;
        result->dirType = DirectionType();
        return false;
    };
public:
    IDF(DirectionType t) : m_type(t) { };
    virtual ~IDF() { };
    
    virtual Spectrum sample(const IDFSample &smp, IDFQueryResult* result) const = 0;
    virtual Spectrum evaluate(const Vector3D &dirIn) const = 0;
    virtual float evaluatePDF(const Vector3D &dirIn) const = 0;
};


class Fresnel {
public:
    virtual Spectrum evaluate(float cosEnter) const = 0;
};

class FresnelNoOp : public Fresnel {
public:
    Spectrum evaluate(float cosEnter) const override;
};

class FresnelConductor : public Fresnel {
    Spectrum m_eta;
    Spectrum m_k;
public:
    FresnelConductor(const Spectrum &eta, const Spectrum &k) : m_eta(eta), m_k(k) { };
    
    Spectrum evaluate(float cosEnter) const override;
};

class FresnelDielectric : public Fresnel {
    float m_etaExt;
    float m_etaInt;
public:
    FresnelDielectric(float etaExt, float etaInt) : m_etaExt(etaExt), m_etaInt(etaInt) { };
    
    float etaExt() const { return m_etaExt; };
    float etaInt() const { return m_etaInt; };
    
    Spectrum evaluate(float cosEnter) const override;
};

#endif

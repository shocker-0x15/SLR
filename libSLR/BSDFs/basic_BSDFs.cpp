//
//  basic_BSDFs.cpp
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "basic_BSDFs.h"
#include "../Core/distributions.h"

namespace SLR {
    SampledSpectrum LambertianBRDF::sampleInternal(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult* result) const {
        result->dir_sn = cosineSampleHemisphere(uDir[0], uDir[1]);
        result->dirPDF = result->dir_sn.z / M_PI;
        result->dirType = m_type;
        result->dir_sn.z *= dot(query.dir_sn, query.gNormal_sn) > 0 ? 1 : -1;
        SampledSpectrum fs = m_R / M_PI;
//        if (dot(result->dir_sn, query.gNormal_sn) < 0)
//            fs = SampledSpectrum::Zero;
        if (result->reverse) {
            result->reverse->fs = fs;
            result->reverse->dirPDF = std::fabs(query.dir_sn.z) / M_PI;
        }
        return fs;
    }
    
    SampledSpectrum LambertianBRDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const {
        if (query.dir_sn.z * dir.z <= 0.0f) {
            SampledSpectrum fs = SampledSpectrum::Zero;
            if (rev_fs)
                *rev_fs = fs;
            return fs;
        }
        SampledSpectrum fs = m_R / M_PI;
        if (rev_fs)
            *rev_fs = fs;
        return fs;
    }
    
    SampledSpectrum LambertianBRDF::evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* revPDF) const {
        if (query.dir_sn.z * dir.z <= 0.0f) {
            if (revPDF)
                *revPDF = 0.0f;
            return 0.0f;
        }
        if (revPDF)
            *revPDF = std::fabs(query.dir_sn.z) / M_PI;
        return std::abs(dir.z) / M_PI;
    }
    
    SampledSpectrum LambertianBRDF::weightInternal(const SLR::BSDFQuery &query) const {
        return m_R;
    }
    
    SampledSpectrum LambertianBRDF::getBaseColorInternal(DirectionType flags) const {
        return m_R;
    }
    
    
    
    SampledSpectrum SpecularBRDF::sampleInternal(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult* result) const {
        result->dir_sn = Vector3D(-query.dir_sn.x, -query.dir_sn.y, query.dir_sn.z);
        result->dirPDF = 1.0f;
        result->dirType = m_type;
        SampledSpectrum fs = m_coeffR * m_fresnel.evaluate(query.dir_sn.z) / std::fabs(query.dir_sn.z);
        if (result->reverse) {
            result->reverse->fs = fs;
            result->reverse->dirPDF = 1.0f;
        }
        return fs;
    }
    
    SampledSpectrum SpecularBRDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const {
        if (rev_fs)
            *rev_fs = SampledSpectrum::Zero;
        return SampledSpectrum::Zero;
    }
    
    SampledSpectrum SpecularBRDF::evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* revPDF) const {
        if (revPDF)
            *revPDF = 0.0f;
        return 0.0f;
    }
    
    SampledSpectrum SpecularBRDF::weightInternal(const BSDFQuery &query) const {
        return m_coeffR * m_fresnel.evaluate(query.dir_sn.z);
    }
    
    SampledSpectrum SpecularBRDF::getBaseColorInternal(DirectionType flags) const {
        return m_coeffR;
    }
    
    
    
    SampledSpectrum SpecularBSDF::sampleInternal(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult* result) const {
        SampledSpectrum F = m_fresnel.evaluate(query.dir_sn.z);
        float reflectProb = F[query.heroIndex];
        if (query.flags.isReflection())
            reflectProb = 1.0f;
        if (query.flags.isTransmission())
            reflectProb = 0.0f;
        if (uComponent < reflectProb) {
            result->dir_sn = Vector3D(-query.dir_sn.x, -query.dir_sn.y, query.dir_sn.z);
            result->dirPDF = SampledSpectrum::Zero;
            result->dirPDF[query.heroIndex] = reflectProb;
            result->dirType = DirectionType::Reflection | DirectionType::Delta0D;
            SampledSpectrum fs = m_coeff * F / std::fabs(query.dir_sn.z);
            if (result->reverse) {
                result->reverse->fs = fs;
                result->reverse->dirPDF = SampledSpectrum::Zero;
                result->reverse->dirPDF[query.heroIndex] = reflectProb;
            }
            return fs;
        }
        else {
            bool entering = query.dir_sn.z > 0.0f;
            float eEnter = entering ? m_fresnel.etaExt()[query.heroIndex] : m_fresnel.etaInt()[query.heroIndex];
            float eExit = entering ? m_fresnel.etaInt()[query.heroIndex] : m_fresnel.etaExt()[query.heroIndex];
            
            float sinEnter2 = 1.0f - query.dir_sn.z * query.dir_sn.z;
            float rrEta = eEnter / eExit;// reciprocal of relative IOR.
            float sinExit2 = rrEta * rrEta * sinEnter2;
            
//            if (sinExit2 >= 1.0f) {
//                result->dirPDF = 0.0f;
//                return SampledSpectrum::Zero;
//            }
            float cosExit = std::sqrt(std::fmax(0.0f, 1.0f - sinExit2));
            if (entering)
                cosExit = -cosExit;
            result->dir_sn = Vector3D(rrEta * -query.dir_sn.x, rrEta * -query.dir_sn.y, cosExit);
            result->dirPDF = SampledSpectrum::Zero;
            result->dirPDF[query.heroIndex] = 1.0f - reflectProb;
            result->dirType = DirectionType::Transmission | DirectionType::Delta0D | (m_type.isDispersive() ? DirectionType::Dispersive : DirectionType());
            SampledSpectrum ret = SampledSpectrum::Zero;
            ret[query.heroIndex] = m_coeff[query.heroIndex] * (1.0f - F[query.heroIndex]) * std::pow(query.adjoint ? (eExit / eEnter) : (eEnter / eExit), 2);
            if (result->reverse) {
                result->reverse->fs = ret / std::fabs(query.dir_sn.z);
                result->reverse->dirPDF = SampledSpectrum::Zero;
                result->reverse->dirPDF[query.heroIndex] = reflectProb;
            }
            SampledSpectrum fs = ret / cosExit;
            return fs;
        }
    }
        
    SampledSpectrum SpecularBSDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const {
        if (rev_fs)
            *rev_fs = SampledSpectrum::Zero;
        return SampledSpectrum::Zero;
    }
    
    SampledSpectrum SpecularBSDF::evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* revPDF) const {
        if (revPDF)
            *revPDF = 0.0f;
        return 0.0f;
    }
    
    SampledSpectrum SpecularBSDF::weightInternal(const SLR::BSDFQuery &query) const {
        return m_coeff;
    }
    
    SampledSpectrum SpecularBSDF::getBaseColorInternal(DirectionType flags) const {
        return m_coeff;
    }
    
    
    
    SampledSpectrum InverseBSDF::sampleInternal(const BSDFQuery &query, float uComponent, const float uDir[2], BSDFQueryResult* result) const {
        BSDFQuery mQuery = query;
        mQuery.flags = mQuery.flags.flip();
        BSDFSample smp{uComponent, uDir[0], uDir[1]};
        SampledSpectrum ret = m_baseBSDF->sample(mQuery, smp, result);
        result->dirType = result->dirType.flip();
        result->dir_sn.z *= -1;
        return ret;
    }
        
    SampledSpectrum InverseBSDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const {
        BSDFQuery mQuery = query;
        mQuery.flags = mQuery.flags.flip();
        Vector3D mDir = dir;
        mDir.z *= -1;
        return m_baseBSDF->evaluate(mQuery, mDir, rev_fs);
    }
    
    SampledSpectrum InverseBSDF::evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* revPDF) const {
        BSDFQuery mQuery = query;
        mQuery.flags.flip();
        Vector3D mDir = dir;
        mDir.z *= -1;
        return m_baseBSDF->evaluatePDF(mQuery, mDir, revPDF);
    }
    
    SampledSpectrum InverseBSDF::weightInternal(const SLR::BSDFQuery &query) const {
        BSDFQuery mQuery = query;
        mQuery.flags = mQuery.flags.flip();
        return m_baseBSDF->weight(mQuery);
    }
    
    SampledSpectrum InverseBSDF::getBaseColorInternal(DirectionType flags) const {
        return m_baseBSDF->getBaseColor(flags.flip());
    }
}

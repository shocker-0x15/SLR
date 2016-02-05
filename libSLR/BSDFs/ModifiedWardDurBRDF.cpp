//
//  ModifiedWardDurBRDF.cpp
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "ModifiedWardDurBRDF.h"

namespace SLR {
    SampledSpectrum ModifiedWardDurBRDF::sampleInternal(const BSDFQuery &query, const BSDFSample &smp, BSDFQueryResult *result) const {
        float quad = 2 * M_PI * smp.uDir[1];
        float phi_h = std::atan2(m_anisoY * std::sin(quad), m_anisoX * std::cos(quad));
        float cosphi_ax = std::cos(phi_h) / m_anisoX;
        float sinphi_ay = std::sin(phi_h) / m_anisoY;
        float theta_h = std::atan(std::sqrt(-std::log(1 - smp.uDir[0]) / (cosphi_ax * cosphi_ax + sinphi_ay * sinphi_ay)));
        Vector3D halfv = Vector3D(std::sin(theta_h) * std::cos(phi_h), std::sin(theta_h) * std::sin(phi_h), std::cos(theta_h));
        halfv.z *= query.dir_sn.z > 0 ? 1 : -1;
        result->dir_sn = 2 * dot(query.dir_sn, halfv) * halfv - query.dir_sn;
        if (result->dir_sn.z * query.dir_sn.z <= 0) {
            result->dirPDF = 0.0f;
            return SampledSpectrum::Zero;
        }
        
        float hx_ax = halfv.x / m_anisoX;
        float hy_ay = halfv.y / m_anisoY;
        float dotHN = std::fabs(halfv.z);
        float dotHI = dot(halfv, result->dir_sn);
        float numerator = std::exp(-(hx_ax * hx_ax + hy_ay * hy_ay) / (dotHN * dotHN));
        float commonDenom = 4 * M_PI * m_anisoX * m_anisoY * dotHI * dotHN * dotHN * dotHN;
        
        result->dirPDF = numerator / commonDenom;
        result->dirType = m_type;
        SampledSpectrum fs = m_R * (numerator / (commonDenom * dotHI * dotHN));
        if (result->reverse) {
            result->reverse->fs = fs;
            result->reverse->dirPDF = result->dirPDF;
        }
        return fs;
    }
    
    SampledSpectrum ModifiedWardDurBRDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const {
        if (dir.z * query.dir_sn.z <= 0) {
            if (rev_fs)
                *rev_fs = SampledSpectrum::Zero;
            return SampledSpectrum::Zero;
        }
        Vector3D halfv = normalize(query.dir_sn + dir);
        float hx_ax = halfv.x / m_anisoX;
        float hy_ay = halfv.y / m_anisoY;
        float dotHN = std::fabs(halfv.z);
        float dotHI = dot(halfv, dir);
        float numerator = std::exp(-(hx_ax * hx_ax + hy_ay * hy_ay) / (dotHN * dotHN));
        float denominator = 4 * M_PI * m_anisoX * m_anisoY * dotHI * dotHI * dotHN * dotHN * dotHN * dotHN;
        SampledSpectrum fs = m_R * numerator / denominator;
        if (rev_fs)
            *rev_fs = fs;
        return fs;
    }
    
    float ModifiedWardDurBRDF::evaluatePDFInternal(const BSDFQuery &query, const Vector3D &dir, float* revPDF) const {
        if (dir.z * query.dir_sn.z <= 0) {
            if (revPDF)
                *revPDF = 0.0f;
            return 0.0f;
        }
        Vector3D halfv = normalize(query.dir_sn + dir);
        float hx_ax = halfv.x / m_anisoX;
        float hy_ay = halfv.y / m_anisoY;
        float dotHN = std::fabs(halfv.z);
        float dotHI = dot(halfv, dir);
        float numerator = std::exp(-(hx_ax * hx_ax + hy_ay * hy_ay) / (dotHN * dotHN));
        float denominator = 4 * M_PI * m_anisoX * m_anisoY * dotHI * dotHN * dotHN * dotHN;
        float ret = numerator / denominator;
        if (revPDF)
            *revPDF = ret;
        return ret;
    }
    
    float ModifiedWardDurBRDF::weightInternal(const BSDFQuery &query, const BSDFSample &smp) const {
        BSDFQueryResult result;
        float fs = sample(query, smp, &result)[query.wlHint];
        return result.dirPDF > 0.0f ? fs * std::fabs(result.dir_sn.z) / result.dirPDF : 0.0f;
    }
    
    float ModifiedWardDurBRDF::weightInternal(const BSDFQuery &query, const Vector3D &dir, float* revWeight) const {
        if (revWeight) {
            SampledSpectrum rev_fs;
            float revDirPDF;
            float fs = evaluate(query, dir, &rev_fs)[query.wlHint];
            float dirPDF = evaluatePDF(query, dir, &revDirPDF);
            if (dirPDF > 0.0f) {
                *revWeight = rev_fs[query.wlHint] * std::fabs(query.dir_sn.z) / revDirPDF;
                return fs * std::fabs(dir.z) / dirPDF;
            }
            else {
                *revWeight = 0.0f;
                return 0.0f;
            }
        }
        else {
            BSDFQueryResult result;
            float fs = evaluate(query, dir)[query.wlHint];
            float dirPDF = evaluatePDF(query, dir);
            return dirPDF > 0.0f ? fs * std::fabs(dir.z) / dirPDF : 0.0f;
        }
    }
    
    float ModifiedWardDurBRDF::weightInternal(const SLR::BSDFQuery &query) const {
        return m_R[query.wlHint];
    }
    
    SampledSpectrum ModifiedWardDurBRDF::getBaseColorInternal(DirectionType flags) const {
        return m_R;
    }    
}

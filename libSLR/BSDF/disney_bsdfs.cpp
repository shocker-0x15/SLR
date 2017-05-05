//
//  disney_bsdfs.cpp
//
//  Created by 渡部 心 on 2017/05/04.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "disney_bsdfs.h"

namespace SLR {
    SampledSpectrum DisneyBRDF::evaluateInternal(const BSDFQuery &query, const Vector3D &dir, SampledSpectrum* rev_fs) const {
        Vector3D halfVector = normalize(query.dirLocal + dir);
        float dotNH = halfVector.z;
        float dotLH = dot(dir, halfVector);
        
        float fresnelL = std::pow(1 - dir.z, 5);
        float fresnelV = std::pow(1 - query.dirLocal.z, 5);
        
        // Diffuse Term
        float F_D90 = 0.5f + 2 * m_roughness * dotLH * dotLH;
        float diffuseFresnel = (1 + (F_D90 - 1) * fresnelL) * (1 + (F_D90 - 1) * fresnelV);
        
        // Subsurface Term
        float F_SS90 = m_roughness * dotLH * dotLH;
        float ssFresnel = (1 + (F_SS90 - 1) * fresnelL) * (1 + (F_SS90 - 1) * fresnelV);
        float ssCoeff = 1.25f * (ssFresnel * (1.0f / (dir.z + query.dirLocal.z) - 0.5f) + 0.5f);
        
        // Specular Term
        float aspect = std::sqrt(1 - 0.9f * m_anisotropic);
        float alpha_x = std::max(0.001f, std::sqrt(m_roughness) / aspect);
        float alpha_y = std::max(0.001f, std::sqrt(m_roughness) * aspect);
        
        SLRAssert_NotImplemented();
        return SampledSpectrum::Zero;
    }
}

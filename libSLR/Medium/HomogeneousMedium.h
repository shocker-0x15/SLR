//
//  HomogeneousMedium.h
//
//  Created by 渡部 心 on 2016/12/31.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef __SLR_HomogeneousMedium__
#define __SLR_HomogeneousMedium__

#include "../defines.h"
#include "../references.h"
#include "../Core/geometry.h"

namespace SLR {
    class HomogeneousMedium : public Medium {
        BoundingBox3D m_region;
        const InputSpectrum* m_sigma_s;
        const InputSpectrum* m_sigma_e;
    public:
        HomogeneousMedium(const BoundingBox3D &region, const InputSpectrum* sigma_s, const InputSpectrum* sigma_e) :
        Medium(sigma_e->calcBounds()), m_region(region), m_sigma_s(sigma_s), m_sigma_e(sigma_e) { }
        
        bool subdivide(Allocator* mem, Medium** fragments, uint32_t* numFragments) const override { return false; }
        
        BoundingBox3D bounds() const override { return m_region; }
        bool contains(const Point3D &p) const override { return m_region.contains(p); }
        bool intersectBoundary(const Ray &ray, float* distToBoundary, bool* enter) const override {
            return m_region.intersectBoundary(ray, distToBoundary, enter);
        }
        SampledSpectrum extinctionCoefficient(const Point3D &p, const WavelengthSamples &wls) const override {
            if (m_region.contains(p))
                return m_sigma_e->evaluate(wls);
                return SampledSpectrum::Zero;
        }
        bool interact(const Ray &ray, const WavelengthSamples &wls, LightPathSampler &pathSampler,
                      MediumInteraction* mi, SampledSpectrum* medThroughput, bool* singleWavelength) const override;
        void getMediumPoint(const MediumInteraction &mi, MediumPoint* medPt) const override;
        void queryCoefficients(const Point3D &p, const WavelengthSamples &wls, SampledSpectrum* sigma_s, SampledSpectrum* sigma_e) const override {
            if (!contains(p)) {
                *sigma_e = SampledSpectrum::Zero;
                *sigma_s = SampledSpectrum::Zero;
            }
            *sigma_s = m_sigma_s->evaluate(wls);
            *sigma_e = m_sigma_e->evaluate(wls);
        }
        float volume() const override { return m_region.volume(); }
        void sample(float u0, float u1, float u2, MediumPoint* medPt, float* volumePDF) const override;
        float evaluateVolumePDF(const MediumPoint& medPt) const override;
    };
}

#endif /* __SLR_HomogeneousMedium__ */

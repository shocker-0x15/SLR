//
//  BidirectionalPathTracingRenderer.h
//
//  Created by 渡部 心 on 2016/01/31.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef __SLR_BidirectionalPathTracingRenderer__
#define __SLR_BidirectionalPathTracingRenderer__

#include "../defines.h"
#include "../references.h"
#include "../Core/Renderer.h"

#include "../Core/geometry.h"
#include "../Core/directional_distribution_functions.h"

namespace SLR {
    class SLR_API BidirectionalPathTracingRenderer : public Renderer {
        struct DDFProxy {
            virtual const void* getDDF() const = 0;
            virtual SampledSpectrum evaluate(const Vector3D &dir_sn, SampledSpectrum* revVal) const = 0;
            virtual float evaluatePDF(const Vector3D &dir_sn, float* revVal = nullptr) const = 0;
        };
        struct EDFProxy : public DDFProxy {
            const EDF* edf;
            EDFQuery query;
            
            EDFProxy(const EDF* _edf, const EDFQuery &_query) : edf(_edf), query(_query) {}
            const void* getDDF() const override { return edf; }
            SampledSpectrum evaluate(const Vector3D &dir_sn, SampledSpectrum* revVal) const override {
                return edf->evaluate(query, dir_sn);
            }
            float evaluatePDF(const Vector3D &dir_sn, float* revVal) const override {
                return edf->evaluatePDF(query, dir_sn);
            }
        };
        struct BSDFProxy : public DDFProxy {
            const BSDF* bsdf;
            BSDFQuery query;
            
            BSDFProxy(const BSDF* _bsdf, const BSDFQuery &_query) : bsdf(_bsdf), query(_query) {}
            const void* getDDF() const override { return bsdf; }
            SampledSpectrum evaluate(const Vector3D &dir_sn, SampledSpectrum* revVal) const override {
                return bsdf->evaluate(query, dir_sn, revVal);
            }
            float evaluatePDF(const Vector3D &dir_sn, float* revVal) const override {
                return bsdf->evaluatePDF(query, dir_sn, revVal);
            }
        };
        struct IDFProxy : public DDFProxy {
            const IDF* idf;
            
            IDFProxy(const IDF* _idf) : idf(_idf) {}
            const void* getDDF() const override { return idf; }
            SampledSpectrum evaluate(const Vector3D &dir_sn, SampledSpectrum* revVal) const override {
                return idf->evaluate(dir_sn);
            }
            float evaluatePDF(const Vector3D &dir_sn, float* revVal) const override {
                return idf->evaluatePDF(dir_sn);
            }
        };
        
        struct BPTVertex {
            SurfacePoint surfPt;
            const DDFProxy* ddf;
            SampledSpectrum alpha;
            float cosIn;
            float areaPDF;
            float RRProb;
            float revAreaPDF;
            float revRRProb;
            DirectionType sampledType;
            bool lambdaSelected;
            BPTVertex(const SurfacePoint &_surfPt, const DDFProxy* _ddf,
                      const SampledSpectrum &_alpha, float _cosIn, float _areaPDF, float _RRProb, DirectionType _sampledType, bool _lambdaSelected) :
            surfPt(_surfPt), ddf(_ddf),
            alpha(_alpha), cosIn(_cosIn), areaPDF(_areaPDF), RRProb(_RRProb), revAreaPDF(NAN), revRRProb(NAN), sampledType(_sampledType), lambdaSelected(_lambdaSelected) {}
        };
        
        struct Job {
            const Scene* scene;
            
            ArenaAllocator* mems;
            IndependentLightPathSampler* pathSamplers;
            
            const Camera* camera;
            ImageSensor* sensor;
            float timeStart;
            float timeEnd;
            uint32_t imageWidth;
            uint32_t imageHeight;
            uint32_t numPixelX;
            uint32_t numPixelY;
            uint32_t basePixelX;
            uint32_t basePixelY;
            
            // working area
            float curPx, curPy;
            int16_t wlHint;
            std::vector<BPTVertex> lightVertices;
            std::vector<BPTVertex> eyeVertices;
            
            ProgressReporter* reporter;
            
            void kernel(uint32_t threadID);
            void generateSubPath(const WavelengthSamples &initWLs, const SampledSpectrum &initAlpha, const Ray &initRay, float dirPDF, DirectionType sampledType,
                                 float cosLast, bool adjoint, IndependentLightPathSampler &pathSampler, ArenaAllocator &mem);
            float calculateMISWeight(float lExtend1stAreaPDF, float lExtend1stRRProb, float lExtend2ndAreaPDF, float lExtend2ndRRProb,
                                     float eExtend1stAreaPDF, float eExtend1stRRProb, float eExtend2ndAreaPDF, float eExtend2ndRRProb,
                                     uint32_t numLVtx, uint32_t numEVtx) const;
        };
        
        uint32_t m_samplesPerPixel;
    public:
        BidirectionalPathTracingRenderer(uint32_t spp);
        void render(const Scene &scene, const RenderSettings &settings) const override;
    };
}

#endif /* __SLR_BidirectionalPathTracingRenderer__ */

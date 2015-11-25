//
//  AMCMCPPMRenderer.h
//
//  Created by 渡部 心 on 2015/08/17.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__AMCMCPPMRenderer__
#define __SLR__AMCMCPPMRenderer__

#include "../defines.h"
#include "../references.h"
#include "../Core/Renderer.h"

#include "../Core/geometry.h"
#include "../Core/RandomNumberGenerator.h"
#include "../Core/directional_distribution_functions.h"

namespace SLR {
    struct Particle {
        Point3D position;
        BoundingBox3D::Axis plane;
        Particle() { };
        Particle(const Point3D &pos) : position(pos) { };
    };
    
    class KDTree {
    protected:
        BoundingBox3D m_BBox;
        uint32_t m_numStored;
        uint32_t m_numHalfStored;
        
        std::vector<Particle*> m_ptrParticles;
        
        void balanceSegment(std::vector<Particle*> &balanced, uint32_t index, uint32_t start, uint32_t end);
        void locateParticlesInternal(uint32_t index, const Point3D &queryPos, float radius2, std::vector<Particle*> &particlesFound) const;
    public:
        void initialize(uint64_t expectedSize);
        uint32_t numStored() const { return m_numStored; };
        void balance();
        void locateParticles(const Point3D &queryPos, float radius2, std::vector<Particle*> &particlesFound) const {
            particlesFound.clear();
            if (m_numStored < 1)
                return;
            locateParticlesInternal(1, queryPos, radius2, particlesFound);
        };
    };
    
    class AMCMCPPMRenderer : public Renderer {
        struct Hitpoint : public Particle {
            float imgX, imgY;
            Normal3D gNormal;
            Vector3D dirOut_sn;
            SampledSpectrum weight;
            const BSDF* bsdf;
            ReferenceFrame shadingFrame;
            
            Hitpoint(float px, float py, const Point3D &pos, const Normal3D &gn, const Vector3D &dirO_sn, const SampledSpectrum &w, const BSDF* f, const ReferenceFrame &frame) :
            imgX(px), imgY(py), Particle(pos), gNormal(gn), dirOut_sn(dirO_sn), weight(w), bsdf(f), shadingFrame(frame) { };
        };
        
        class HitpointMap : public KDTree {
            std::vector<std::vector<Hitpoint>> m_points;
            uint32_t m_numThreads;
        public:
            HitpointMap() {};
            ~HitpointMap() {};
            
            void initialize(uint32_t numThreads, uint32_t numPixels);
            void store(uint32_t thread, float imgX, float imgY,
                       const SampledSpectrum &weight, const Point3D &pos, const Normal3D &gn, const Vector3D &dir_sn,
                       const BSDF* f, const ReferenceFrame &frame);
            void build();
            void queryHitpoints(const Point3D &pos, const Normal3D &gn_sn, float radius, std::vector<Particle*> &pointsFound) const;
        };
        
        struct PrimarySample {
        protected:
            void adaptiveMutateElement(float *value, float u1, float u2, float size) {
                if (u1 < 0.5f) {
                    *value -= powf(u2, 1.0f / size + 1.0f);
                    if (*value < 0.0f) *value += 1.0f;
                }
                else {
                    *value += powf(u2, 1.0f / size + 1.0f);
                    if (*value >= 1.0f) *value -= 1.0f;
                }
            };
        public:
            uint64_t modifiedTime;
            
            PrimarySample() : modifiedTime(0) { };
            
            virtual void init(RandomNumberGenerator &rng) = 0;
            virtual void adaptiveMutate(RandomNumberGenerator &rng, float size) = 0;
        };
        
        struct PathPrimarySample : public PrimarySample {
            float uComponent;
            float uDirection[2];
            float uRR;
            
            void init(RandomNumberGenerator &rng) {
                uComponent = rng.getFloat0cTo1o();
                uDirection[0] = rng.getFloat0cTo1o();
                uDirection[1] = rng.getFloat0cTo1o();
                uRR = rng.getFloat0cTo1o();
            };
            
            void adaptiveMutate(RandomNumberGenerator &rng, float size) {
                adaptiveMutateElement(&uComponent, rng.getFloat0cTo1o(), rng.getFloat0cTo1o(), size);
                adaptiveMutateElement(&uDirection[0], rng.getFloat0cTo1o(), rng.getFloat0cTo1o(), size);
                adaptiveMutateElement(&uDirection[1], rng.getFloat0cTo1o(), rng.getFloat0cTo1o(), size);
                adaptiveMutateElement(&uRR, rng.getFloat0cTo1o(), rng.getFloat0cTo1o(), size);
            };
        };
        
        struct LightPrimarySample : public PrimarySample {
            float uLight;
            float uPosition[2];
            float uComponent;
            float uDirection[2];
            
            void init(RandomNumberGenerator &rng) {
                uLight = rng.getFloat0cTo1o();
                uPosition[0] = rng.getFloat0cTo1o();
                uPosition[1] = rng.getFloat0cTo1o();
                uComponent = rng.getFloat0cTo1o();
                uDirection[0] = rng.getFloat0cTo1o();
                uDirection[1] = rng.getFloat0cTo1o();
            };
            
            void adaptiveMutate(RandomNumberGenerator &rng, float size) {
                adaptiveMutateElement(&uLight, rng.getFloat0cTo1o(), rng.getFloat0cTo1o(), size);
                adaptiveMutateElement(&uPosition[0], rng.getFloat0cTo1o(), rng.getFloat0cTo1o(), size);
                adaptiveMutateElement(&uPosition[1], rng.getFloat0cTo1o(), rng.getFloat0cTo1o(), size);
                adaptiveMutateElement(&uComponent, rng.getFloat0cTo1o(), rng.getFloat0cTo1o(), size);
                adaptiveMutateElement(&uDirection[0], rng.getFloat0cTo1o(), rng.getFloat0cTo1o(), size);
                adaptiveMutateElement(&uDirection[1], rng.getFloat0cTo1o(), rng.getFloat0cTo1o(), size);
            };
        };
        
        class ReplicaExchangeSampler {
            RandomNumberGenerator* m_rng;
            uint64_t m_globalTime;
            uint64_t m_uniformTime;
            bool m_uniform;
            float m_mSize;
            std::vector<float> m_mSizeHistory;
            uint32_t m_histIdx;
            
            LightPrimarySample m_lightSample;
            LightPrimarySample m_lightStack;
            
            uint32_t m_curPathDim;
            std::vector<PathPrimarySample> m_pathSamples;
            std::stack<PathPrimarySample> m_pathStack;
        public:
            enum class Mode {
                Uniform,
                Mutation
            };
            
            ReplicaExchangeSampler() { };
            ReplicaExchangeSampler(RandomNumberGenerator* rng);
            ~ReplicaExchangeSampler() { };
            
            LightPrimarySample getLightPrimarySample();
            PathPrimarySample getPathPrimarySample() ;
            void initDimensions(Mode mode, float size);
            void clearStack();
            void rollBack();
        };
        
        struct DistributedRTJob {
            const Scene* scene;
            WavelengthSamples wls;
            ArenaAllocator* mems;
            RandomNumberGenerator** rngs;
            HitpointMap* hpMap;
            
            const Camera* camera;
            float time;
            
            ImageSensor* sensor;
            uint32_t imageWidth;
            uint32_t imageHeight;
            uint32_t numPixelX;
            uint32_t numPixelY;
            uint32_t basePixelX;
            uint32_t basePixelY;
            
            void kernel(uint32_t threadID);
            SampledSpectrum record(uint32_t threadID, const Scene &scene, const WavelengthSamples &initWLs, float px, float py, const Ray &initRay, const SampledSpectrum &initAlpha,
                                   RandomNumberGenerator &rng, ArenaAllocator &mem) const;
        };
        struct PhotonSplattingJob {
            struct ResultRecord {
                float imgX, imgY;
                SampledSpectrum contribution;
                ResultRecord() { };
                ResultRecord(float px, float py, const SampledSpectrum &c) : imgX(px), imgY(py), contribution(c) { };
            };
            const Scene* scene;
            WavelengthSamples wls;
            float radius;
            ArenaAllocator* mems;
            HitpointMap* hpMap;
            uint32_t numPhotons;
            
            ImageSensor* sensor;
            float time;
            
            ReplicaExchangeSampler sampler;
            uint64_t uniformCount;
            uint64_t mutationCount;
            uint64_t acceptedCount;
            float mutationSize;
            std::vector<ResultRecord> prevResults;
            
            void kernel(uint32_t threadID);
            float photonTracing(const WavelengthSamples &wls, ArenaAllocator &mem, std::vector<ResultRecord> &results);
        };
        
        uint32_t m_numPhotonsPerPass;
        uint32_t m_numPasses;
    public:
        AMCMCPPMRenderer(uint32_t numPhotonsPerPass, uint32_t numPasses);
        void render(const Scene &scene, const RenderSettings &settings) const override;
    };    
}

#endif

//
//  AMCMCPPMRenderer.cpp
//
//  Created by 渡部 心 on 2015/08/17.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "AMCMCPPMRenderer.h"

#include "../Core/RenderSettings.h"
#include "../Helper/ThreadPool.h"
#include "../Core/XORShift.h"
#include "../Core/ImageSensor.h"
#include "../Core/RandomNumberGenerator.h"
#include "../Memory/ArenaAllocator.h"
#include "../Core/cameras.h"
#include "../Core/geometry.h"
#include "../Core/SurfaceObject.h"
#include "../Core/directional_distribution_functions.h"

namespace SLR {
    void KDTree::initialize(uint64_t expectedSize) {
        m_numStored = 0;
        
        m_ptrParticles.clear();
        m_ptrParticles.reserve(expectedSize);
        
        m_ptrParticles.push_back(NULL);
        
        m_BBox = BoundingBox3D();
    }
    
    void KDTree::balance() {
        if (m_numStored > 1) {
            std::vector<Particle*> balancedList(m_numStored + 1);
            balanceSegment(balancedList, 1, 1, m_numStored);
            std::copy(balancedList.begin(), balancedList.end(), m_ptrParticles.begin());
        }
        
        m_numHalfStored = m_numStored / 2 - 1;
    }
    
    void KDTree::balanceSegment(std::vector<Particle*> &balanced, uint32_t index, uint32_t start, uint32_t end) {
        int median = 1;
        // medianを整列する領域の長さの半分に最も近い2のべき乗とする。（この段階では完全にそうはならない）
        // length: 383, median: 128(128)
        // length: 384, median: 128(128, 256)
        // length: 385, median: 128(256)
        // length: 511, median: 128(256)
        // length: 512, median: 256(256)
        // length: 513, median: 256(256)
        // length: 640, median: 256(256)
        // length: 768, median: 256(256, 512)
        // length: 800, median: 256(512)
        while ((4 * median) <= (end - start + 1))
            median += median;
        
        // medianを領域の長さの半分に最も近い2のべき乗にして、適切なオフセットを加える。
        // start: 1, end: 383, median: end - 128 + 1: 255
        // start: 1, end: 384, median: 256 + start - 1: 256
        // start: 1, end: 385, median: 256 + start - 1: 256
        // start: 1, end: 511, median: 256 + start - 1: 256
        // start: 1, end: 512, median: end - 256 + 1: 257
        // start: 1, end: 513, median: end - 256 + 1: 258
        // start: 1, end: 640, median: end - 256 + 1: 485
        // start: 1, end: 768, median: 512 + start - 1: 512
        // start: 1, end: 800, median: 512 + start - 1: 512
        // length : 8 ~ 11: end-, 12 ~ 15: +start,
        //         16 ~ 23: end-, 24 ~ 31: +start,
        //         32 ~ 47: end-, 48 ~ 63: +start,
        //         2^n ~ 2^n + 2^(n-1) - 1: end-
        //         2^n + 2^(n-1) ~ 2^(n+1) - 1: +start
        if ((3 * median) <= (end - start + 1)) {
            median += median;
            median += start - 1;
        }
        else {
            median = end - median + 1;
        }
        
        BoundingBox3D::Axis axis = m_BBox.widestAxis();
        
        // axis軸においてメディアンとなる値より小さい値を前に、大きい値を後ろに並べる。
        auto compareFunc = [&axis](Particle* p0, Particle* p1) { return p0->position[axis] < p1->position[axis]; };
        std::nth_element(m_ptrParticles.begin() + start, m_ptrParticles.begin() + median, m_ptrParticles.begin() + end + 1, compareFunc);
        
        balanced[index] = m_ptrParticles[median];
        balanced[index]->plane = axis;
        
        // balance the left side segment.
        if (start < median) {
            if (start < median - 1) {
                const float tmp = m_BBox.maxP[axis];
                m_BBox.maxP[axis] = balanced[index]->position[axis];
                balanceSegment(balanced, 2 * index, start, median - 1);
                m_BBox.maxP[axis] = tmp;
            }
            else {
                balanced[2 * index] = m_ptrParticles[start];
            }
        }
        
        // balance the right side segment.
        if (median < end) {
            if (median + 1 < end) {
                const float tmp = m_BBox.minP[axis];
                m_BBox.minP[axis] = balanced[index]->position[axis];
                balanceSegment(balanced, 2 * index + 1, median + 1, end);
                m_BBox.minP[axis] = tmp;
            }
            else {
                balanced[2 * index + 1] = m_ptrParticles[end];
            }
        }
    }
    
    void KDTree::locateParticlesInternal(uint32_t index, const Point3D &queryPos, float radius2, std::vector<Particle*> &particlesFound) const {
        Particle* p = m_ptrParticles[index];
        
        // Recursively search.
        if (index < m_numHalfStored) {
            float vDist = queryPos[p->plane] - p->position[p->plane];
            if (vDist > 0.0) {
                locateParticlesInternal(2 * index + 1, queryPos, radius2, particlesFound);
                if (vDist * vDist < radius2)
                    locateParticlesInternal(2 * index, queryPos, radius2, particlesFound);
            }
            else {
                locateParticlesInternal(2 * index, queryPos, radius2, particlesFound);
                if (vDist * vDist < radius2)
                    locateParticlesInternal(2 * index + 1, queryPos, radius2, particlesFound);
            }
        }
        
        // store the particle if it is in the search radius.
        if (sqDistance(queryPos, p->position) < radius2)
            particlesFound.push_back(p);
    }
    
    
    void AMCMCPPMRenderer::HitpointMap::initialize(uint32_t numThreads, uint32_t numPixels) {
        KDTree::initialize(numPixels * 3);
        m_numThreads = numThreads;
        m_points.resize(m_numThreads);
        for (int i = 0; i < m_numThreads; ++i) {
            m_points[i].clear();
            m_points[i].reserve(numPixels * 3 / numThreads);
        }
    }
    
    void AMCMCPPMRenderer::HitpointMap::store(uint32_t thread, float imgX, float imgY,
                                              const SampledSpectrum &weight, const Point3D &pos, const Normal3D &gn, const Vector3D &dir_sn,
                                              const BSDF* f, const ReferenceFrame &frame) {
        m_points[thread].emplace_back(imgX, imgY, pos, gn, dir_sn, weight, f, frame);
    }
    
    void AMCMCPPMRenderer::HitpointMap::build() {
        for (int i = 0; i < m_numThreads; ++i) {
            for (int j = 0; j < m_points[i].size(); ++j) {
                m_ptrParticles.push_back(&(m_points[i][j]));
                m_BBox.unify(m_points[i][j].position);
            }
            m_numStored += m_points[i].size();
        }
        balance();
    }
    
    void AMCMCPPMRenderer::HitpointMap::queryHitpoints(const Point3D &pos, const Normal3D &gn, float radius, std::vector<Particle*> &pointsFound) const {
        locateParticles(pos, radius * radius, pointsFound);
        for (int i = 0; i < pointsFound.size(); ++i) {
            Hitpoint &hp = *(Hitpoint*)pointsFound[i];
            if (dot(Vector3D(gn), Vector3D(hp.gNormal)) < 0.707)
                pointsFound[i] = nullptr;
        }
    }
    
    
    AMCMCPPMRenderer::ReplicaExchangeSampler::ReplicaExchangeSampler(RandomNumberGenerator* rng) : m_rng(rng) {
        m_pathSamples.resize(100);
        m_globalTime = 0;
        m_curPathDim = 0;
        m_mSizeHistory.resize(100);
    };
    
    AMCMCPPMRenderer::LightPrimarySample AMCMCPPMRenderer::ReplicaExchangeSampler::getLightPrimarySample() {
        if (m_lightSample.modifiedTime < m_globalTime) {
            if (m_uniform) {
                m_lightStack = m_lightSample;
                m_lightSample.init(*m_rng);
                m_lightSample.modifiedTime = m_globalTime;
            }
            else {
                if (m_lightSample.modifiedTime < m_uniformTime) {
                    m_lightSample.init(*m_rng);
                    m_lightSample.modifiedTime = m_uniformTime;
                }
                
                int mIdx = 0;
                while (m_lightSample.modifiedTime < m_globalTime - 1) {
                    m_lightSample.adaptiveMutate(*m_rng, m_mSizeHistory[mIdx++]);
                    ++m_lightSample.modifiedTime;
                }
                
                m_lightStack = m_lightSample;
                m_lightSample.adaptiveMutate(*m_rng, m_mSize);
                ++m_lightSample.modifiedTime;
            }
        }
        
        return m_lightSample;
    };
    
    AMCMCPPMRenderer::PathPrimarySample AMCMCPPMRenderer::ReplicaExchangeSampler::getPathPrimarySample() {
        if (m_curPathDim >= m_pathSamples.size())
            m_pathSamples.resize(m_pathSamples.size() * 1.5f);
        
        PathPrimarySample &curSample = m_pathSamples[m_curPathDim];
        if (curSample.modifiedTime < m_globalTime) {
            if (m_uniform) {
                m_pathStack.push(curSample);
                curSample.init(*m_rng);
                curSample.modifiedTime = m_globalTime;
            }
            else {
                if (curSample.modifiedTime < m_uniformTime) {
                    curSample.init(*m_rng);
                    curSample.modifiedTime = m_uniformTime;
                }
                
                int mIdx = 0;
                while (curSample.modifiedTime < m_globalTime - 1) {
                    curSample.adaptiveMutate(*m_rng, m_mSizeHistory[mIdx++]);
                    ++curSample.modifiedTime;
                }
                
                m_pathStack.push(curSample);
                curSample.adaptiveMutate(*m_rng, m_mSize);
                ++curSample.modifiedTime;
            }
        }
        
        return m_pathSamples[m_curPathDim++];
    };
    
    // パスの生成前に呼ばれる．
    void AMCMCPPMRenderer::ReplicaExchangeSampler::initDimensions(Mode mode, float size) {
        ++m_globalTime;
        m_curPathDim = 0;
        m_uniform = mode == Mode::Uniform;
        m_mSize = size;
    };
    
    // 採択時に呼ばれる．
    void AMCMCPPMRenderer::ReplicaExchangeSampler::clearStack() {
        while (m_pathStack.empty() == false)
            m_pathStack.pop();
        
        if (m_uniform) {
            // 一様サンプルとレプリカ交換を行った場合は時間を記録する．
            m_uniformTime = m_globalTime;
            m_histIdx = 0;
        }
        else {
            // 変異させた場合は変異サイズの履歴を保存する．
            m_mSizeHistory[m_histIdx++] = m_mSize;
            if (m_histIdx >= m_mSizeHistory.size())
                m_mSizeHistory.resize(m_mSizeHistory.size() * 1.5f);
        }
    };
    
    // 棄却時に呼ばれる．
    void AMCMCPPMRenderer::ReplicaExchangeSampler::rollBack() {
        m_lightSample = m_lightStack;
        while (m_pathStack.empty() == false) {
            m_pathSamples[--m_curPathDim] = m_pathStack.top();
            m_pathStack.pop();
        }
        --m_globalTime;
    };
    
    
    AMCMCPPMRenderer::AMCMCPPMRenderer() {
    };
    
    void AMCMCPPMRenderer::render(const Scene &scene, const RenderSettings &settings) const {        
        uint32_t numSamples = settings.getInt(RenderSettingItem::NumSamples);
        uint32_t imageWidth = settings.getInt(RenderSettingItem::ImageWidth);
        uint32_t imageHeight = settings.getInt(RenderSettingItem::ImageHeight);
        ImageSensor sensor(imageWidth, imageHeight);
        
#ifdef DEBUG
        uint32_t numThreads = 1;
#else
        uint32_t numThreads = std::thread::hardware_concurrency();
#endif
        XORShift topRand(settings.getInt(RenderSettingItem::RNGSeed));
        auto mems = std::unique_ptr<ArenaAllocator[]>(new ArenaAllocator[numThreads]);
        auto memPTs = std::unique_ptr<ArenaAllocator[]>(new ArenaAllocator[numThreads]);
        auto rngs = std::unique_ptr<XORShift[]>(new XORShift[numThreads]);
        auto rngRefs = std::unique_ptr<RandomNumberGenerator*[]>(new RandomNumberGenerator*[numThreads]);
        auto scales = std::unique_ptr<float[]>(new float[numThreads]);
        for (int i = 0; i < numThreads; ++i) {
            new (mems.get() + i) ArenaAllocator();
            new (memPTs.get() + i) ArenaAllocator();
            rngRefs[i] = new (rngs.get() + i) XORShift(topRand.getUInt());
        }
        sensor.addSeparatedBuffers(numThreads);
        
        HitpointMap hitpointMap;
        
        DistributedRTJob jobDRT;
        jobDRT.scene = &scene;
        jobDRT.sensor = &sensor;
        jobDRT.camera = scene.getCamera();
        jobDRT.imageWidth = imageWidth;
        jobDRT.imageHeight = imageHeight;
        jobDRT.numPixelX = sensor.tileWidth();
        jobDRT.numPixelY = sensor.tileHeight();
        jobDRT.mems = mems.get();
        jobDRT.rngs = rngRefs.get();
        jobDRT.hpMap = &hitpointMap;
        
        float radius = 0.0025f;
        const uint32_t numMCMCs = 16;
        const uint32_t numPhotonsPerPass = 100000;
        PhotonSplattingJob jobPSs[numMCMCs];
        std::unique_ptr<XORShift[]> rngPSs = std::unique_ptr<XORShift[]>(new XORShift[numMCMCs]);
        for (int i = 0; i < numMCMCs; ++i) {
            PhotonSplattingJob &jobPS = jobPSs[i];
            jobPS.scene = &scene;
            jobPS.sensor = &sensor;
            jobPS.mems = memPTs.get();
            jobPS.hpMap = &hitpointMap;
            jobPS.numPhotons = numPhotonsPerPass;
            
            jobPS.sampler = ReplicaExchangeSampler(rngPSs.get());
            jobPS.uniformCount = 0;
            jobPS.mutationCount = 0;
            jobPS.acceptedCount = 0;
            jobPS.mutationSize = 1;
        }
        
        uint32_t imgIdx = 0;
        uint32_t exportIdx = 1;
        uint32_t endIdx = 16;
        
        float timeStart = settings.getFloat(RenderSettingItem::TimeStart);
        float timeEnd = settings.getFloat(RenderSettingItem::TimeEnd);
        float sensorReponce = settings.getFloat(RenderSettingItem::SensorResponse);
        
        for (int s = 0; s < numSamples; ++s) {
            float time = timeStart + (timeEnd - timeStart) * topRand.getFloat0cTo1o();
            jobDRT.time = time;
            
            // FIXME: wavelength selection should be done in each path.
            float selectWLPDF;
            WavelengthSamples wls = WavelengthSamples::createWithEqualOffsets(topRand.getFloat0cTo1o(), topRand.getFloat0cTo1o(), &selectWLPDF);
            
            jobDRT.wls = wls;
            
            wls.flags |= WavelengthSamples::LambdaSelected;
            for (int i = 0; i < numMCMCs; ++i)
                jobPSs[i].wls = wls;
            
            hitpointMap.initialize(numThreads, imageWidth * imageHeight);
            
            // Distributed Ray Tracing Pass: record hitpoints in a k-d tree.
            ThreadPool DRTPool(numThreads);
            for (int ty = 0; ty < sensor.numTileY(); ++ty) {
                for (int tx = 0; tx < sensor.numTileX(); ++tx) {
                    jobDRT.basePixelX = tx * sensor.tileWidth();
                    jobDRT.basePixelY = ty * sensor.tileHeight();
                    DRTPool.enqueue(std::bind(&DistributedRTJob::kernel, jobDRT, std::placeholders::_1));
                }
            }
            DRTPool.wait();
            
            // Build a balanced k-d tree.
            hitpointMap.build();
            
            // Photon Tracing Pass: splatting photon's contribution to hitpoints near the photon.
            ThreadPool PTPool(numThreads);
            for (int i = 0; i < numMCMCs; ++i) {
                jobPSs[i].time = time;
                jobPSs[i].radius = radius;
                PTPool.enqueue(std::bind(&PhotonSplattingJob::kernel, std::ref(jobPSs[i]), std::placeholders::_1));
            }
            PTPool.wait();
            
            for (int i = 0; i < numThreads; ++i)
                mems[i].reset();
            
            uint64_t numPasses = s + 1;
            const float alpha = 2.0f / 3.0f;
            radius = std::sqrt((numPasses + alpha) / (numPasses + 1)) * radius;
            
            uint64_t numPhotons = numPasses * numPhotonsPerPass * numMCMCs;
            uint64_t uniformCount = 0;
            for (int i = 0; i < numMCMCs; ++i) {
                printf("MCMC: %2u: #uniform: %llu, #mutation: %llu, #accept: %llu, mutationSize: %g\n",
                       i, jobPSs[i].uniformCount, jobPSs[i].mutationCount, jobPSs[i].acceptedCount, jobPSs[i].mutationSize);
                uniformCount += jobPSs[i].uniformCount;
            }
            float scale = float(uniformCount) / numPhotons * (sensorReponce / (numPasses * numMCMCs));
            for (int i = 0; i < numThreads; ++i)
                scales[i] = scale;
            printf("\n");
            
            ++imgIdx;
            if (imgIdx == exportIdx) {
                exportIdx += exportIdx;
                char filename[256];
                sprintf(filename, "%03u.bmp", imgIdx);
                sensor.saveImage(filename, settings.getFloat(RenderSettingItem::SensorResponse) / (s + 1));
                printf("%u samples: %s\n", s, filename);
                ++imgIdx;
                if (imgIdx == endIdx)
                    break;
            }
        }
        
        //    sensor.saveImage("output.png", settings.getFloat(RenderSettingItem::SensorResponse) / numSamples);
    }
    
    
    void AMCMCPPMRenderer::DistributedRTJob::kernel(uint32_t threadID) {
        ArenaAllocator &mem = mems[threadID];
        RandomNumberGenerator &rng = *rngs[threadID];
        for (int ly = 0; ly < numPixelY; ++ly) {
            for (int lx = 0; lx < numPixelX; ++lx) {
                float px = basePixelX + lx + rng.getFloat0cTo1o();
                float py = basePixelY + ly + rng.getFloat0cTo1o();
                
                LensPosQuery lensQuery(time, wls);
                LensPosSample lensSample(rng.getFloat0cTo1o(), rng.getFloat0cTo1o());
                LensPosQueryResult lensResult;
                SampledSpectrum We0 = camera->sample(lensQuery, lensSample, &lensResult);
                
                IDFSample WeSample(px / imageWidth, py / imageHeight);
                IDFQueryResult WeResult;
                IDF* idf = camera->createIDF(lensResult.surfPt, wls, mem);
                SampledSpectrum We1 = idf->sample(WeSample, &WeResult);
                
                Ray ray(lensResult.surfPt.p, lensResult.surfPt.shadingFrame.fromLocal(WeResult.dirLocal), time);
                SampledSpectrum weight = (We0 * We1) / (lensResult.areaPDF * WeResult.dirPDF);
                SampledSpectrum C = record(threadID, *scene, wls, px, py, ray, weight, rng, mem);
                SLRAssert(weight.hasNaN() == false && weight.hasInf() == false,
                          "Unexpected value detected: %s\n"
                          "pix: (%f, %f)", weight.toString().c_str(), px, py);
                SLRAssert(C.hasNaN() == false && C.hasInf() == false,
                          "Unexpected value detected: %s\n"
                          "pix: (%f, %f)", C.toString().c_str(), px, py);
                
                sensor->add(px, py, wls, C);
            }
        }
    }
    
    SampledSpectrum AMCMCPPMRenderer::DistributedRTJob::record(uint32_t threadID, const Scene &scene, const WavelengthSamples &initWLs, float px, float py, const Ray &initRay, const SampledSpectrum &initAlpha,
                                                               RandomNumberGenerator &rng, ArenaAllocator &mem) const {
        WavelengthSamples wls = initWLs;
        Ray ray = initRay;
        SurfacePoint surfPt;
        SampledSpectrum alpha = initAlpha;
        float initY = alpha.luminance();
        SampledSpectrumSum sp(SampledSpectrum::Zero);
        uint32_t pathLength = 0;
        
        Intersection isect;
        if (!scene.intersect(ray, &isect))
            return SampledSpectrum::Zero;
        isect.getSurfacePoint(&surfPt);
        
        Vector3D dirOut_sn = surfPt.shadingFrame.toLocal(-ray.dir);
        if (surfPt.isEmitting()) {
            EDF* edf = surfPt.createEDF(wls, mem);
            SampledSpectrum Le = surfPt.emittance(wls) * edf->evaluate(EDFQuery(), dirOut_sn);
            sp += alpha * Le;
        }
        if (surfPt.atInfinity)
            return sp;
        
        while (true) {
            ++pathLength;
            Normal3D gNorm_sn = surfPt.shadingFrame.toLocal(surfPt.gNormal);
            BSDF* bsdf = surfPt.createBSDF(wls, mem);
            
            DirectionType receiveType = DirectionType(DirectionType::WholeSphere | DirectionType::NonDelta);
            if (bsdf->matches(receiveType))
                hpMap->store(threadID, px, py, alpha, surfPt.p, surfPt.gNormal, dirOut_sn, bsdf, surfPt.shadingFrame);
            
            DirectionType scatterType = DirectionType(DirectionType::WholeSphere | DirectionType::Delta);
            if (!bsdf->matches(scatterType))
                break;
            
            BSDFQuery fsQuery(dirOut_sn, gNorm_sn, wls.selectedLambda);
            BSDFSample fsSample(rng.getFloat0cTo1o(), rng.getFloat0cTo1o(), rng.getFloat0cTo1o());
            fsQuery.flags = scatterType;
            BSDFQueryResult fsResult;
            SampledSpectrum fs = bsdf->sample(fsQuery, fsSample, &fsResult);
            if (fs == SampledSpectrum::Zero || fsResult.dirPDF == 0.0f)
                break;
            if (fsResult.dirType.isDispersive()) {
                fsResult.dirPDF /= WavelengthSamples::NumComponents;
                wls.flags |= WavelengthSamples::LambdaSelected;
            }
            alpha *= fs * (std::fabs(fsResult.dir_sn.z) / fsResult.dirPDF);
            SLRAssert(!alpha.hasInf() && !alpha.hasNaN(), "alpha: unexpected value detected: %s", alpha.toString().c_str());
            
            Vector3D dirIn = surfPt.shadingFrame.fromLocal(fsResult.dir_sn);
            ray = Ray(surfPt.p + Ray::Epsilon * dirIn, dirIn, ray.time);
            
            isect = Intersection();
            if (!scene.intersect(ray, &isect))
                break;
            isect.getSurfacePoint(&surfPt);
            
            dirOut_sn = surfPt.shadingFrame.toLocal(-ray.dir);
            
            if (surfPt.isEmitting()) {
                EDF* edf = surfPt.createEDF(wls, mem);
                SampledSpectrum Le = surfPt.emittance(wls) * edf->evaluate(EDFQuery(), dirOut_sn);
                SLRAssert(!Le.hasNaN() && !Le.hasInf(), "Le: unexpected value detected: %s", Le.toString().c_str());
                
                sp += alpha * Le;
            }
            if (surfPt.atInfinity)
                break;
            
            // Russian roulette
            float continueProb = std::min(alpha.luminance() / initY, 1.0f);
            if (rng.getFloat0cTo1o() < continueProb)
                alpha /= continueProb;
            else
                break;
        }
        
        return sp;
    }
    
    
    void AMCMCPPMRenderer::PhotonSplattingJob::kernel(uint32_t threadID) {
        ArenaAllocator &mem = mems[threadID];
        std::vector<ResultRecord> results;
        
        // 一様サンプリングでフォトンパスが可視になるまで続ける．
        if (uniformCount == 0) {
            while (true) {
                results.clear();
                sampler.initDimensions(ReplicaExchangeSampler::Mode::Uniform, 0);
                if (photonTracing(wls, mem, results) != 0.0f) {
                    sampler.clearStack();
                    break;
                }
                sampler.rollBack();
            }
            uniformCount = 1;
            acceptedCount = 1;
            
            prevResults.resize(results.size());
            std::copy(results.begin(), results.end(), prevResults.begin());
        }
        
        FloatSum cMutationSize = mutationSize;
        const float targetAcceptanceRate = 0.234f;
        for (int i = 0; i < numPhotons; ++i) {
            results.clear();
            sampler.initDimensions(ReplicaExchangeSampler::Mode::Uniform, 0);
            bool success = true;
            float V = photonTracing(wls, mem, results);
            if (V != 0.0f) {
                // 一様サンプリングパスが可視だった場合．
                ++uniformCount;
            }
            else {
                // 一様サンプリングパスが不可視だった場合は，前回のパスに戻して変異を行う．
                sampler.rollBack();
                ++mutationCount;
                
                sampler.initDimensions(ReplicaExchangeSampler::Mode::Mutation, mutationSize < 0.0f ? 0.0f : mutationSize);
                V = photonTracing(wls, mem, results);
                if (V != 0.0f)
                    ++acceptedCount;
                else
                    success = false;
                
                // 変異サイズを採択率に応じて適応的に変更する．
                float acceptanceRate = float(acceptedCount) / mutationCount;
                cMutationSize += (acceptanceRate - targetAcceptanceRate) / mutationCount;
            }
            
            if (success) {
                // 一様サンプリングもしくは変異によってパスを変化させることに成功した場合．
                // 新たなパスの寄与の蓄積と今後のための寄与の保存を行う．
                for (int j = 0; j < results.size(); ++j)
                    sensor->add(threadID, results[j].imgX, results[j].imgY, wls, results[j].contribution);
                
                prevResults.resize(results.size());
                std::copy(results.begin(), results.end(), prevResults.begin());
                
                sampler.clearStack();
            }
            else {
                // パスの変化に失敗した場合は前回の寄与をそのまま蓄積する．
                for (int j = 0; j < prevResults.size(); ++j)
                    sensor->add(threadID, prevResults[j].imgX, prevResults[j].imgY, wls, prevResults[j].contribution);
                
                sampler.rollBack();
            }
            
            mem.reset();
        }
        mutationSize = cMutationSize;
    }
    
    float AMCMCPPMRenderer::PhotonSplattingJob::photonTracing(const WavelengthSamples &wls, ArenaAllocator &mem, std::vector<ResultRecord> &results) {
        float I = 0.0f;
        
        LightPrimarySample psLight = sampler.getLightPrimarySample();
        Light light;
        float lightProb;
        scene->selectLight(psLight.uLight, &light, &lightProb);
        
        LightPosQuery xpQuery(time, wls);
        LightPosSample xpSample(psLight.uPosition[0], psLight.uPosition[1]);
        LightPosQueryResult xpResult;
        EDFQuery LeQuery;
        EDFSample LeSample(psLight.uComponent, psLight.uDirection[0], psLight.uDirection[1]);
        EDFQueryResult LeResult;
        SampledSpectrum M = light.sample(xpQuery, xpSample, &xpResult);
        EDF* edf = xpResult.surfPt.createEDF(wls, mem);
        SampledSpectrum Le = M * edf->sample(LeQuery, LeSample, &LeResult);
        
        SampledSpectrum alpha = Le / (lightProb * xpResult.areaPDF * LeResult.dirPDF) * std::fabs(LeResult.dir_sn.z) / numPhotons;
        Vector3D dirOut = xpResult.surfPt.shadingFrame.fromLocal(LeResult.dir_sn);
        Ray ray(xpResult.surfPt.p + Ray::Epsilon * dirOut, dirOut, time);
        float initY = alpha.luminance();
        uint32_t pathLength = 0;
        
        std::vector<Particle*> hitpoints;
        while (true) {
            ++pathLength;
            Intersection isect;
            SurfacePoint surfPt;
            if (!scene->intersect(ray, &isect) || std::isinf(ray.distMax))
                break;
            isect.getSurfacePoint(&surfPt);
            
            // FIXME: should I consider the delta check using hitpoint's BSDF?
            BSDF* bsdf = surfPt.createBSDF(wls, mem);
            if (bsdf->hasNonDelta()) {
                hpMap->queryHitpoints(surfPt.p, surfPt.gNormal, radius, hitpoints);
                float kernelWeight = 1.0f / (M_PI * radius * radius);
                for (int i = 0; i < hitpoints.size(); ++i) {
                    if (!hitpoints[i])
                        continue;
                    Hitpoint &hp = *(Hitpoint*)hitpoints[i];
                    Vector3D dirIn_sn = hp.shadingFrame.toLocal(-ray.dir);
                    Normal3D gNorm_sn = hp.shadingFrame.toLocal(hp.gNormal);
                    BSDFQuery queryBSDF(dirIn_sn, gNorm_sn, wls.selectedLambda);
                    SampledSpectrum fs = hp.bsdf->evaluate(queryBSDF, hp.dirOut_sn);
                    SampledSpectrum contribution = hp.weight * fs * alpha * kernelWeight * std::fabs(dirIn_sn.z / dot(Vector3D(surfPt.gNormal), ray.dir));
                    SLRAssert(!contribution.hasInf() && !contribution.hasNaN(), "contribution: unexpected value detected: %s",
                              contribution.toString().c_str());
                    results.emplace_back(hp.imgX, hp.imgY, contribution);
                    I = 1.0;
                }
            }
            
            Vector3D dirIn_sn = surfPt.shadingFrame.toLocal(-ray.dir);
            Normal3D gNorm_sn = surfPt.shadingFrame.toLocal(surfPt.gNormal);
            
            // TODO: consider adjoint correction.
            PathPrimarySample psPath = sampler.getPathPrimarySample();
            BSDFQuery fsQuery(dirIn_sn, gNorm_sn, wls.selectedLambda);
            BSDFSample fsSample(psPath.uComponent, psPath.uDirection[0], psPath.uDirection[1]);
            BSDFQueryResult fsResult;
            SampledSpectrum fs = bsdf->sample(fsQuery, fsSample, &fsResult);
            if (fs == SampledSpectrum::Zero || fsResult.dirPDF == 0.0f)
                break;
            alpha *= fs * (std::fabs(fsResult.dir_sn.z) / fsResult.dirPDF);
            SLRAssert(!alpha.hasInf() && !alpha.hasNaN(), "alpha: unexpected value detected: %s", alpha.toString().c_str());
            
            dirOut = surfPt.shadingFrame.fromLocal(fsResult.dir_sn);
            ray = Ray(surfPt.p + Ray::Epsilon * dirOut, dirOut, ray.time);
            
            // Russian roulette
            float continueProb = std::min(alpha.luminance() / initY, 1.0f);
            if (psPath.uRR < continueProb)
                alpha /= continueProb;
            else
                break;
        }
        
        return I;
    }    
}

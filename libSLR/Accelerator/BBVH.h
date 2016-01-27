//
//  BBVH.h
//
//  Created by 渡部 心 on 2014/05/04.
//  Copyright (c) 2014年 渡部 心. All rights reserved.
//

#ifndef SLR_BBVH_h
#define SLR_BBVH_h

#include "../defines.h"
#include "../references.h"
#include "../Core/geometry.h"
#include "../Core/SurfaceObject.h"

namespace SLR {
    class SLR_API BBVH {
    public:
        enum class Partitioning {
            Median = 0,
            Midpoint = 1,
            BinnedSAH = 2
        };
    private:
        struct Node {
            BoundingBox3D bbox;
            uint32_t c0, c1;
            uint32_t offsetFirstLeaf;
            uint32_t numLeaves;
            BoundingBox3D::Axis axis;
            
            Node() : c0(0), c1(0), offsetFirstLeaf(0), numLeaves(0) { };
            
            void initAsLeaf(const BoundingBox3D &bb, uint32_t offset, uint32_t numLs) {
                bbox = bb;
                c0 = c1 = 0;
                offsetFirstLeaf = offset;
                numLeaves = numLs;
            };
            void initAsInternal(const BoundingBox3D &bb, uint32_t cIdx0, uint32_t cIdx1, BoundingBox3D::Axis ax) {
                bbox = bb;
                c0 = cIdx0;
                c1 = cIdx1;
                axis = ax;
                offsetFirstLeaf = numLeaves = 0;
            };
        };
        
        struct ObjInfos {
            std::vector<SurfaceObject*>* objs;
            std::vector<BoundingBox3D> bboxes;
            std::vector<Point3D> centroids;
            std::vector<uint32_t> indices;
        };
        
        Partitioning m_method;
        uint32_t m_depth;
        BoundingBox3D m_bounds;
        std::vector<Node> m_nodes;
        std::vector<SurfaceObject*> m_objLists;
        
        uint32_t buildRecursive(ObjInfos &infos, uint32_t start, uint32_t end, uint32_t depth) {
            auto &indices = infos.indices;
            auto &centroids = infos.centroids;
            
            uint32_t nodeIdx = (uint32_t)m_nodes.size();
            m_nodes.emplace_back();
            
            if (++depth > m_depth)
                m_depth = depth;
            
            BoundingBox3D bbox;
            BoundingBox3D centroidBB;
            for (uint32_t i = start; i < end; ++i) {
                uint32_t idx = indices[i];
                bbox.unify(infos.bboxes[idx]);
                centroidBB.unify(centroids[idx]);
            }
            
            uint32_t numObjs = end - start;
            SLRAssert(numObjs >= 1, "Number of objects is zero.");
            BoundingBox3D::Axis widestAxis = centroidBB.widestAxis();
            
            if (numObjs == 1) {
                m_nodes[nodeIdx].initAsLeaf(bbox, (uint32_t)m_objLists.size(), 1);
                m_objLists.push_back(infos.objs->at(indices[start]));
                return nodeIdx;
            }
            
            uint32_t splitIdx;
            switch (m_method) {
                case Partitioning::Median: { // 子の数が同じになるように分割する。
                    splitIdx = (start + end) / 2;
                    std::nth_element(indices.begin() + start, indices.begin() + splitIdx, indices.begin() + end, [&centroids, &widestAxis](uint32_t idx0, uint32_t idx1) {
                        return centroids[idx0][widestAxis] < centroids[idx1][widestAxis];
                    });
                    break;
                }
                case Partitioning::Midpoint: { // 重心点の範囲中で分布が広い次元の中点で分割する。
                    float pivot = centroidBB.centerOfAxis(widestAxis);
                    auto firstOf2ndGroup = std::partition(indices.begin() + start, indices.begin() + end, [&centroids, &widestAxis, &pivot](uint32_t idx) {
                        return centroids[idx][widestAxis] < pivot;
                    });
                    splitIdx = std::max((uint32_t)std::distance(indices.begin() + start, firstOf2ndGroup), 1u) + start;
                    break;
                }
                case Partitioning::BinnedSAH: { // SAH(の近似)
                    struct BinInfo {
                        BoundingBox3D bbox;
                        uint32_t numObjs;
                        BinInfo() : numObjs(0) { };
                    };
                    
                    const float travCost = 0.125f;
                    const float isectCost = 1.0f;
                    const uint32_t numBins = 16;
                    BinInfo binInfos[numBins];
                    
                    for (uint32_t i = start; i < end; ++i) {
                        uint32_t idx = indices[i];
                        uint32_t bin = numBins * ((centroids[idx][widestAxis] - centroidBB.minP[widestAxis]) / (centroidBB.maxP[widestAxis] - centroidBB.minP[widestAxis]));
                        bin = std::min(bin, numBins - 1);
                        ++binInfos[bin].numObjs;
                        binInfos[bin].bbox.unify(infos.bboxes[idx]);
                    }
                    
                    float surfaceAreaParent = bbox.surfaceArea();
                    float costs[numBins - 1];
                    for (uint32_t i = 0; i < numBins - 1; ++i) {
                        BoundingBox3D b0, b1;
                        uint32_t count0 = 0, count1 = 0;
                        for (int j = 0; j <= i; ++j) {
                            b0.unify(binInfos[j].bbox);
                            count0 += binInfos[j].numObjs;
                        }
                        for (int j = i + 1; j < numBins; ++j) {
                            b1.unify(binInfos[j].bbox);
                            count1 += binInfos[j].numObjs;
                        }
                        costs[i] = travCost + (b0.surfaceArea() * isectCost * count0 + b1.surfaceArea() * isectCost * count1) / surfaceAreaParent;
                    }
                    
                    float minCost = costs[0];
                    uint32_t splitPlane = 0;
                    for (int i = 1; i < numBins - 1; ++i) {
                        if (costs[i] < minCost) {
                            minCost = costs[i];
                            splitPlane = i;
                        }
                    }
                    
                    if (minCost < numObjs * isectCost) {
                        float pivot = centroidBB.minP[widestAxis] + (centroidBB.maxP[widestAxis] - centroidBB.minP[widestAxis]) / numBins * (splitPlane + 1);
                        auto firstOf2ndGroup = std::partition(indices.begin() + start, indices.begin() + end, [&centroids, &widestAxis, &pivot](uint32_t idx) {
                            return centroids[idx][widestAxis] < pivot;
                        });
                        splitIdx = std::max((uint32_t)std::distance(indices.begin() + start, firstOf2ndGroup), 1u) + start;
                    }
                    else {
                        m_nodes[nodeIdx].initAsLeaf(bbox, (uint32_t)m_objLists.size(), numObjs);
                        for (uint32_t i = start; i < end; ++i)
                            m_objLists.push_back(infos.objs->at(indices[i]));
                        return nodeIdx;
                    }
                    break;
                }
                default:
                    break;
            }
            
            uint32_t c0 = buildRecursive(infos, start, splitIdx, depth);
            uint32_t c1 = buildRecursive(infos, splitIdx, end, depth);
            m_nodes[nodeIdx].initAsInternal(bbox, c0, c1, widestAxis);
            return nodeIdx;
        }
        
        float calcSAHCost() const {
            const float Ci = 0.125f;
            const float Cl = 0.0f;
            const float Ct = 1.0f;
            float costInt = 0.0f;
            float costLeaf = 0.0f;
            float costObj = 0.0f;
            for (int i = 0; i < m_nodes.size(); ++i) {
                const Node &node = m_nodes[i];
                float surfaceArea = node.bbox.surfaceArea();
                if (node.numLeaves == 0) {
                    costInt += surfaceArea;
                }
                else {
                    costLeaf += surfaceArea;
                    costObj += surfaceArea * node.numLeaves;
                }
            }
            float rootSA = m_nodes[0].bbox.surfaceArea();
            costInt *= Ci / rootSA;
            costLeaf *= Cl / rootSA;
            costObj *= Ct / rootSA;
            
            return costInt + costLeaf + costObj;
        };
        
    public:
        BBVH(std::vector<SurfaceObject*> &objs, Partitioning method = Partitioning::BinnedSAH) {
            m_method = method;
            
            ObjInfos infos;
            infos.objs = &objs;
            for (int i = 0; i < objs.size(); ++i) {
                BoundingBox3D bb = objs[i]->bounds();
                m_bounds.unify(bb);
                infos.bboxes.push_back(bb);
                infos.centroids.push_back(bb.centroid());
                infos.indices.push_back(i);
            }
            
            m_depth = 0;
            buildRecursive(infos, 0, (uint32_t)objs.size(), 0);
#ifdef DEBUG
            printf("depth: %u, cost: %g\n", m_depth, calcSAHCost());
#endif
        };
        
        BoundingBox3D bounds() const {
            return m_bounds;
        };
        
        bool intersect(Ray &ray, Intersection* isect) const {
            uint32_t objDepth = (uint32_t)isect->obj.size();
            Vector3D invRayDir = ray.dir.reciprocal();
            bool dirSigns[] = {ray.dir.x > 0, ray.dir.y > 0, ray.dir.z > 0};
            auto intersectAABB = [&ray, &invRayDir](const BoundingBox3D bb) {
                float dist0 = ray.distMin, dist1 = ray.distMax;
                Vector3D tNear = (bb.minP - ray.org) * invRayDir;
                Vector3D tFar = (bb.maxP - ray.org) * invRayDir;
                for (int i = 0; i < 3; ++i) {
                    if (tNear[i] > tFar[i])
                        std::swap(tNear[i], tFar[i]);
                    dist0 = tNear[i] > dist0 ? tNear[i] : dist0;
                    dist1 = tFar[i] < dist1 ? tFar[i]  : dist1;
                    if (dist0 > dist1)
                        return false;
                }
                return true;
            };
            
            uint32_t idxStack[64];
            uint32_t depth = 0;
            idxStack[depth++] = 0;
            while (depth > 0) {
                const Node &node = m_nodes[idxStack[--depth]];
                if (!intersectAABB(node.bbox))
                    continue;
                if (node.numLeaves == 0) {
                    bool positiveDir = dirSigns[node.axis];
                    idxStack[depth++] = positiveDir ? node.c1 : node.c0;
                    idxStack[depth++] = positiveDir ? node.c0 : node.c1;
                }
                else {
                    for (uint32_t i = 0; i < node.numLeaves; ++i)
                        if (m_objLists[node.offsetFirstLeaf + i]->intersect(ray, isect))
                            ray.distMax = isect->dist;
                }
            }
            return isect->obj.size() > objDepth;
        }
        
        bool intersect(Ray &ray, SurfacePoint* surfPt) const {
            Intersection isect;
            if (!intersect(ray, &isect))
                return false;
            isect.obj.top()->getSurfacePoint(isect, surfPt);
            return true;
        }
    };
    
}

#endif

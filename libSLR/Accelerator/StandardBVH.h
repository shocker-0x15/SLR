//
//  StandardBVH.h
//
//  Created by 渡部 心 on 2014/05/04.
//  Copyright (c) 2014年 渡部 心. All rights reserved.
//

#ifndef __SLR_StandardBVH__
#define __SLR_StandardBVH__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/accelerator.h"

namespace SLR {
    class SLR_API StandardBVH : public Accelerator {
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
            const std::vector<SurfaceObject*>* objs;
            std::vector<BoundingBox3D> bboxes;
            std::vector<Point3D> centroids;
            std::vector<uint32_t> indices;
        };
        
        Partitioning m_method;
        uint32_t m_depth;
        float m_cost;
        BoundingBox3D m_bounds;
        std::vector<Node> m_nodes;
        std::vector<const SurfaceObject*> m_objLists;
        
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
            BoundingBox3D::Axis widestAxis = centroidBB.widestAxis();
            const float pcBBMin = centroidBB.minP[widestAxis];
            const float pcBBMax = centroidBB.maxP[widestAxis];
            
            uint32_t numObjs = end - start;
            SLRAssert(numObjs >= 1, "Number of objects is zero.");
            
            if (numObjs == 1) {
                m_nodes[nodeIdx].initAsLeaf(bbox, (uint32_t)m_objLists.size(), 1);
                m_objLists.push_back(infos.objs->at(indices[start]));
                return nodeIdx;
            }
            
            uint32_t splitIdx;
            switch (m_method) {
                case Partitioning::Median: {
                    // partitions so that the numbers of children of both side become the same.
                    splitIdx = (start + end) / 2;
                    std::nth_element(indices.begin() + start, indices.begin() + splitIdx, indices.begin() + end, [&centroids, &widestAxis](uint32_t idx0, uint32_t idx1) {
                        return centroids[idx0][widestAxis] < centroids[idx1][widestAxis];
                    });
                    break;
                }
                case Partitioning::Midpoint: {
                    // partitions at the midpoint of a dimension in which primitives' centroid distribution is the widest.
                    float pivot = centroidBB.centerOfAxis(widestAxis);
                    auto firstOf2ndGroup = std::partition(indices.begin() + start, indices.begin() + end, [&centroids, &widestAxis, &pivot](uint32_t idx) {
                        return centroids[idx][widestAxis] < pivot;
                    });
                    splitIdx = std::max((uint32_t)std::distance(indices.begin() + start, firstOf2ndGroup), 1u) + start;
                    break;
                }
                case Partitioning::BinnedSAH: {
                    if ((pcBBMax - pcBBMin) <= 0) {
                        // partitions so that the numbers of children of both side become the same.
                        splitIdx = (start + end) / 2;
                        std::nth_element(indices.begin() + start, indices.begin() + splitIdx, indices.begin() + end, [&centroids, &widestAxis](uint32_t idx0, uint32_t idx1) {
                            return centroids[idx0][widestAxis] < centroids[idx1][widestAxis];
                        });
                        break;
                    }
                    
                    struct BinInfo {
                        BoundingBox3D bbox;
                        uint32_t numObjs;
                        float sumCost;
                        BinInfo() : numObjs(0), sumCost(0.0f) { };
                    };
                    
                    const float travCost = 1.2f;
                    const uint32_t numBins = 16;
                    BinInfo binInfos[numBins];
                    
                    // Binning and calculate cost of leaf node from all the primitives.
                    float leafNodeCost = 0.0f;
                    for (uint32_t i = start; i < end; ++i) {
                        uint32_t idx = indices[i];
                        float isectCost = infos.objs->at(idx)->costForIntersect();
                        leafNodeCost += isectCost;
                        
                        uint32_t bin = numBins * ((centroids[idx][widestAxis] - pcBBMin) / (pcBBMax - pcBBMin));
                        bin = std::min(bin, numBins - 1);
                        ++binInfos[bin].numObjs;
                        binInfos[bin].sumCost += isectCost;
                        binInfos[bin].bbox.unify(infos.bboxes[idx]);
                    }
                    
                    // evaluate SAH cost for every pair of child partitions and determine a plane with the minimum cost.
                    uint32_t splitPlane = 0;
                    float minCost = INFINITY;
                    float surfaceAreaParent = bbox.surfaceArea();
                    for (uint32_t i = 0; i < numBins - 1; ++i) {
                        BoundingBox3D b0, b1;
                        float cost0 = 0.0f, cost1 = 0.0f;
                        for (int j = 0; j <= i; ++j) {
                            b0.unify(binInfos[j].bbox);
                            cost0 += binInfos[j].sumCost;
                        }
                        for (int j = i + 1; j < numBins; ++j) {
                            b1.unify(binInfos[j].bbox);
                            cost1 += binInfos[j].sumCost;
                        }
                        float cost = travCost + (b0.surfaceArea() * cost0 + b1.surfaceArea() * cost1) / surfaceAreaParent;
                        if (cost < minCost) {
                            minCost = cost;
                            splitPlane = i;
                        }
                    }
                    SLRAssert(std::isfinite(minCost), "invalid cost value: %g", minCost);
                    
                    // perform object partitioning if the cost is less than the leaf node cost. 
                    if (minCost < leafNodeCost) {
                        float pivot = centroidBB.minP[widestAxis] + (pcBBMax - pcBBMin) / numBins * (splitPlane + 1);
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
            const float Ci = 1.2f;
            const float Cl = 0.0f;
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
                    float costPrims = 0.0f;
                    for (uint32_t j = 0; j < node.numLeaves; ++j)
                        costPrims += m_objLists[node.offsetFirstLeaf + j]->costForIntersect();
                    costObj += surfaceArea * costPrims;
                }
            }
            float rootSA = m_nodes[0].bbox.surfaceArea();
            costInt *= Ci / rootSA;
            costLeaf *= Cl / rootSA;
            costObj /= rootSA;
            
            return costInt + costLeaf + costObj;
        }
        
    public:
        StandardBVH(const std::vector<SurfaceObject*> &objs, Partitioning method = Partitioning::BinnedSAH) {
            m_method = method;
            
            ObjInfos infos;
            infos.objs = &objs;
            infos.bboxes.resize(objs.size());
            infos.centroids.resize(objs.size());
            infos.indices.resize(objs.size());
            for (int i = 0; i < objs.size(); ++i) {
                BoundingBox3D bb = objs[i]->bounds();
                m_bounds.unify(bb);
                infos.bboxes[i] = bb;
                infos.centroids[i] = bb.centroid();
                infos.indices[i] = i;
            }
            
            m_depth = 0;
            buildRecursive(infos, 0, (uint32_t)objs.size(), 0);
            m_cost = calcSAHCost();
#ifdef DEBUG
            printf("depth: %u, cost: %g\n", m_depth, m_cost);
#endif
        }
        
        float costForIntersect() const override {
            return m_cost;
        }
        
        BoundingBox3D bounds() const override {
            return m_bounds;
        }
        
        bool intersect(const Ray &ray, const RaySegment &segment, SurfaceInteraction* si, const SurfaceObject** closestObject) const override {
            *closestObject = nullptr;
            bool dirIsPositive[] = {ray.dir.x >= 0, ray.dir.y >= 0, ray.dir.z >= 0};
            RaySegment isectRange = segment;
            
            const uint32_t StackSize = 64;
            uint32_t idxStack[StackSize];
            uint32_t depth = 0;
            idxStack[depth++] = 0;
            while (depth > 0) {
                const Node &node = m_nodes[idxStack[--depth]];
                if (!node.bbox.intersect(ray, isectRange))
                    continue;
                if (node.numLeaves == 0) {
                    SLRAssert(depth < StackSize, "StandardBVH::intersect: stack overflow");
                    bool positiveDir = dirIsPositive[node.axis];
                    idxStack[depth++] = positiveDir ? node.c1 : node.c0;
                    idxStack[depth++] = positiveDir ? node.c0 : node.c1;
                }
                else {
                    for (uint32_t i = 0; i < node.numLeaves; ++i) {
                        if (m_objLists[node.offsetFirstLeaf + i]->intersect(ray, isectRange, si)) {
                            *closestObject = m_objLists[node.offsetFirstLeaf + i];
                            isectRange.distMax = si->getDistance();
                        }
                    }
                }
            }
            return *closestObject != nullptr;
        }
    };
}

#endif /* __SLR_StandardBVH__ */

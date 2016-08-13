//
//  SBVH.h
//
//  Created by 渡部 心 on 2016/07/03.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef SBVH_h
#define SBVH_h

#include "../defines.h"
#include "../references.h"
#include "../Core/Accelerator.h"

#include <nmmintrin.h>

namespace SLR {
    class SLR_API SBVH : public Accelerator {
        friend class QBVH;
        
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
        
        struct Fragment {
            const SurfaceObject* obj;
            BoundingBox3D bbox;
            float costForIntersect;
        };
        
        uint32_t m_depth;
        float m_cost;
        BoundingBox3D m_bounds;
        std::vector<Node> m_nodes;
        std::vector<const SurfaceObject*> m_objLists;
        
        uint32_t buildRecursive(Fragment* fragments, uint32_t currentSize, uint32_t maximumBudget, uint32_t start, uint32_t end, uint32_t depth, uint32_t* numAdded) {
//#define PRINT_PROCESSING_TIME
            *numAdded = 0;
            uint32_t nodeIdx = (uint32_t)m_nodes.size();
            m_nodes.emplace_back();
            
            if (++depth > m_depth)
                m_depth = depth;
            
            std::chrono::system_clock::time_point tpStart, tpEnd;
            double elapsed;
            
            tpStart = std::chrono::system_clock::now();
            
            // calculate AABBs and a cost of leaf node from all the primitives.
            BoundingBox3D parentBB;
            BoundingBox3D parentCentroidBB;
            float leafNodeCost = 0.0f;
            for (uint32_t i = start; i < end; ++i) {
                parentBB.unify(fragments[i].bbox);
                parentCentroidBB.unify(fragments[i].bbox.centroid());
                float isectCost = fragments[i].costForIntersect;
                leafNodeCost += isectCost;
            }
            const BoundingBox3D::Axis widestAxisOP = parentCentroidBB.widestAxis();
            const BoundingBox3D::Axis widestAxisSP = parentBB.widestAxis();
            const float surfaceAreaParent = parentBB.surfaceArea();
            const float pcBBMin = parentCentroidBB.minP[widestAxisOP];
            const float pcBBMax = parentCentroidBB.maxP[widestAxisOP];
            const float pBBMin = parentBB.minP[widestAxisSP];
            const float pBBMax = parentBB.maxP[widestAxisSP];
            
            uint32_t numObjs = end - start;
            SLRAssert(numObjs >= 1, "Number of objects is zero.");
            
            tpEnd = std::chrono::system_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::microseconds>(tpEnd - tpStart).count();
#ifdef PRINT_PROCESSING_TIME
            if (depth == 1 || depth == 2)
                printf("calculate parent BBox (%u-%u): %g[ms]\n", start, end, elapsed * 0.001f);
#endif
            
            if (numObjs == 1) {
                m_nodes[nodeIdx].initAsLeaf(parentBB, (uint32_t)m_objLists.size(), 1);
                m_objLists.push_back(fragments[start].obj);
                return nodeIdx;
            }
            
            const float travCost = 1.2f;
            
            struct ObjectBinInfo {
                BoundingBox3D bbox;
                uint32_t numObjs;
                float sumCost;
                ObjectBinInfo() : numObjs(0), sumCost(0.0f) {}
            };
            
            struct SpatialBinInfo {
                BoundingBox3D bbox;
                uint32_t numEntries;
                uint32_t numExits;
                float sumCostEntries;
                float sumCostExits;
                SpatialBinInfo() : numEntries(0), numExits(0), sumCostEntries(0.0f), sumCostExits(0.0f) {}
            };
            
            const uint32_t numObjBins = 32;
            ObjectBinInfo objBinInfos[numObjBins];
            uint32_t splitPlaneOP = 0;
            float minCostByOP = INFINITY;
            
            if (parentCentroidBB.surfaceArea() > 0) {
                tpStart = std::chrono::system_clock::now();
                
                // Object Binning
                for (uint32_t i = start; i < end; ++i) {
                    const Fragment &fragment = fragments[i];
                    
                    uint32_t binIdx = numObjBins * ((fragment.bbox.centerOfAxis(widestAxisOP) - pcBBMin) / (pcBBMax - pcBBMin));
                    binIdx = std::min(binIdx, numObjBins - 1);
                    
                    ++objBinInfos[binIdx].numObjs;
                    objBinInfos[binIdx].sumCost += fragment.costForIntersect;
                    objBinInfos[binIdx].bbox.unify(fragment.bbox);
                }
                
                // evaluate SAH cost for every pair of child partitions and determine a plane with the minimum cost.
                for (uint32_t i = 0; i < numObjBins - 1; ++i) {
                    BoundingBox3D b0, b1;
                    float cost0 = 0.0f, cost1 = 0.0f;
                    for (int j = 0; j <= i; ++j) {
                        b0.unify(objBinInfos[j].bbox);
                        cost0 += objBinInfos[j].sumCost;
                    }
                    for (int j = i + 1; j < numObjBins; ++j) {
                        b1.unify(objBinInfos[j].bbox);
                        cost1 += objBinInfos[j].sumCost;
                    }
                    float cost = travCost + (b0.surfaceArea() * cost0 + b1.surfaceArea() * cost1) / surfaceAreaParent;
                    if (cost < minCostByOP) {
                        minCostByOP = cost;
                        splitPlaneOP = i;
                    }
                }
                
                tpEnd = std::chrono::system_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::microseconds>(tpEnd - tpStart).count();
#ifdef PRINT_PROCESSING_TIME
                if (depth == 1 || depth == 2)
                    printf("Object Binning: %g[ms]\n", elapsed * 0.001f);
#endif
            }
            
            // calculate surface area of intersection of two bounding boxes resulted from object partitioning.
            // This becomes a criteria for determining whether it performs spatial splitting.
            BoundingBox3D bbLeftOP, bbRightOP;
            for (int j = 0; j <= splitPlaneOP; ++j)
                bbLeftOP.unify(objBinInfos[j].bbox);
            for (int j = splitPlaneOP + 1; j < numObjBins; ++j)
                bbRightOP.unify(objBinInfos[j].bbox);
            BoundingBox3D overlappedBB = intersection(bbLeftOP, bbRightOP);
            float overlappedSA = 0;
            if (overlappedBB.isValid())
                overlappedSA = overlappedBB.surfaceArea();
            
            const uint32_t numSBins = 16;
            const float spatialBinWidth = parentBB.width(widestAxisSP) / numSBins;
            SpatialBinInfo sBinInfos[numSBins];
            uint32_t splitPlaneSP = 0;
            float minCostBySP = INFINITY;
            
            float alpha = 1e-5;
            if (overlappedSA / m_bounds.surfaceArea() > alpha) {
                tpStart = std::chrono::system_clock::now();
                
                // Spatial Binning
                for (int i = start; i < end; ++i) {
                    const Fragment &fragment = fragments[i];
                    
                    const BoundingBox3D &bbox = fragment.bbox;
                    uint32_t entryBin = numSBins * ((bbox.minP[widestAxisSP] - pBBMin) / (pBBMax - pBBMin));
                    uint32_t exitBin = numSBins * ((bbox.maxP[widestAxisSP] - pBBMin) / (pBBMax - pBBMin));
                    entryBin = std::min(entryBin, numSBins - 1);
                    exitBin = std::min(exitBin, numSBins - 1);
                    
                    ++sBinInfos[entryBin].numEntries;
                    ++sBinInfos[exitBin].numExits;
                    
                    float isectCost = fragment.costForIntersect;
                    sBinInfos[entryBin].sumCostEntries += isectCost;
                    sBinInfos[exitBin].sumCostExits += isectCost;
                    
                    for (int binIdx = entryBin; binIdx <= exitBin; ++binIdx) {
                        float splitPosMin = binIdx * spatialBinWidth + pBBMin;
                        BoundingBox3D choppedBB = fragment.obj->choppedBounds(widestAxisSP, splitPosMin, splitPosMin + spatialBinWidth);
                        sBinInfos[binIdx].bbox.unify(intersection(choppedBB, bbox));
                    }
                }
                
                // evaluate SAH cost for every pair of child partitions and determine a plane with the minimum cost.
                for (uint32_t i = 0; i < numSBins - 1; ++i) {
                    BoundingBox3D b0, b1;
                    float cost0 = 0.0f, cost1 = 0.0f;
                    for (int j = 0; j <= i; ++j) {
                        b0.unify(sBinInfos[j].bbox);
                        cost0 += sBinInfos[j].sumCostEntries;
                    }
                    for (int j = i + 1; j < numSBins; ++j) {
                        b1.unify(sBinInfos[j].bbox);
                        cost1 += sBinInfos[j].sumCostExits;
                    }
                    float cost = travCost + (b0.surfaceArea() * cost0 + b1.surfaceArea() * cost1) / surfaceAreaParent;
                    if (cost < minCostBySP) {
                        minCostBySP = cost;
                        splitPlaneSP = i;
                    }
                }
                
                tpEnd = std::chrono::system_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::microseconds>(tpEnd - tpStart).count();
#ifdef PRINT_PROCESSING_TIME
                if (depth == 1 || depth == 2)
                    printf("Spatial Binning: %g[ms]\n", elapsed * 0.001f);
#endif
            }
            
            if (leafNodeCost < minCostByOP && leafNodeCost < minCostBySP) {
                m_nodes[nodeIdx].initAsLeaf(parentBB, (uint32_t)m_objLists.size(), numObjs);
                for (uint32_t i = start; i < end; ++i)
                    m_objLists.push_back(fragments[i].obj);
                return nodeIdx;
            }
            else if (minCostByOP < minCostBySP) {
                tpStart = std::chrono::system_clock::now();
                
                float pivot = pcBBMin + (pcBBMax - pcBBMin) / numObjBins * (splitPlaneOP + 1);
                auto firstOf2ndGroup = std::partition(fragments + start, fragments + end, [&widestAxisOP, &pivot](const Fragment &fragment) {
                    return fragment.bbox.centerOfAxis(widestAxisOP) < pivot;
                });
                uint32_t splitIdx = std::max((uint32_t)std::distance(fragments + start, firstOf2ndGroup), 1u) + start;
                SLRAssert(splitIdx > start && splitIdx < end, "Invalid partitioning.");
                
                tpEnd = std::chrono::system_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::microseconds>(tpEnd - tpStart).count();
#ifdef PRINT_PROCESSING_TIME
                if (depth == 1 || depth == 2)
                    printf("Object Partitioning: %g[ms]\n", elapsed * 0.001f);
#endif
                
                uint32_t numLeftAdded, numRightAdded;
                uint32_t c0 = buildRecursive(fragments, currentSize, maximumBudget, start, splitIdx, depth, &numLeftAdded);
                uint32_t c1 = buildRecursive(fragments, currentSize + numLeftAdded, maximumBudget, splitIdx + numLeftAdded, end + numLeftAdded, depth, &numRightAdded);
                m_nodes[nodeIdx].initAsInternal(parentBB, c0, c1, widestAxisOP);
                *numAdded += numLeftAdded + numRightAdded;
                return nodeIdx;
            }
            else {
                tpStart = std::chrono::system_clock::now();
                
                uint32_t numLeftsBySplit = 0;
                uint32_t numRightsBySplit = 0;
                for (int j = 0; j <= splitPlaneSP; ++j)
                    numLeftsBySplit += sBinInfos[j].numEntries;
                for (int j = splitPlaneSP + 1; j < numSBins; ++j)
                    numRightsBySplit += sBinInfos[j].numExits;
                
                uint32_t numLeftIndices = 0, numRightIndices = 0;
                Fragment* newFragments = new Fragment[numLeftsBySplit + numRightsBySplit];
                Fragment* leftFragments = newFragments + 0;
                Fragment* rightFragments = newFragments + numLeftsBySplit;
                for (int i = start; i < end; ++i) {
                    BoundingBox3D &bbox = fragments[i].bbox;
                    uint32_t entryBin = numSBins * ((bbox.minP[widestAxisSP] - pBBMin) / (pBBMax - pBBMin));
                    uint32_t exitBin = numSBins * ((bbox.maxP[widestAxisSP] - pBBMin) / (pBBMax - pBBMin));
                    entryBin = std::min(entryBin, numSBins - 1);
                    exitBin = std::min(exitBin, numSBins - 1);
                    
//                    // consider unsplitting
//                    float isectCost = fragments[i].obj->costForIntersect();
//                    float splitCost = splitPlane.leftBBox.surfaceArea() * splitPlane.leftCost + splitPlane.rightBBox.surfaceArea() * splitPlane.rightCost;
//                    float leftAlignCost = calcUnion(splitPlane.leftBBox, bbox).surfaceArea() * splitPlane.leftCost + splitPlane.rightBBox.surfaceArea() * (splitPlane.rightCost - isectCost);
//                    float rightAlignCost = splitPlane.leftBBox.surfaceArea() * (splitPlane.leftCost - isectCost) + calcUnion(splitPlane.rightBBox, bbox).surfaceArea() * splitPlane.rightCost;
//                    int32_t cheapest = cheapest = splitCost < leftAlignCost ? (splitCost < rightAlignCost ? 0 : 1) : (leftAlignCost < rightAlignCost ? -1 : 1);
                    
                    if (exitBin <= splitPlaneSP) {
                        leftFragments[numLeftIndices++] = fragments[i];
                    }
                    else if (entryBin > splitPlaneSP) {
                        rightFragments[numRightIndices++] = fragments[i];
                    }
                    else {
                        Fragment &dstLeft = leftFragments[numLeftIndices++];
                        Fragment &dstRight = rightFragments[numRightIndices++];
                        
                        float splitPos = (splitPlaneSP + 1) * spatialBinWidth + pBBMin;
                        BoundingBox3D splitLeftBBox, splitRightBBox;
                        fragments[i].obj->splitBounds(widestAxisSP, splitPos, &splitLeftBBox, &splitRightBBox);
                        
                        dstLeft.obj = dstRight.obj = fragments[i].obj;
                        dstLeft.costForIntersect = dstRight.costForIntersect = fragments[i].costForIntersect;
                        dstLeft.bbox = intersection(splitLeftBBox, bbox);
                        dstRight.bbox = intersection(splitRightBBox, bbox);
                    }
                }
                *numAdded = (numLeftIndices + numRightIndices) - numObjs;
                std::copy_backward(fragments + end, fragments + currentSize, fragments + currentSize + *numAdded);
                std::copy(leftFragments, leftFragments + numLeftIndices, fragments + start);
                std::copy(rightFragments, rightFragments + numRightIndices, fragments + start + numLeftIndices);
                delete[] newFragments;
                uint32_t splitIdx = start + numLeftIndices;
                SLRAssert(splitIdx > start && splitIdx < end + *numAdded, "Invalid partitioning.");
                
                tpEnd = std::chrono::system_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::microseconds>(tpEnd - tpStart).count();
#ifdef PRINT_PROCESSING_TIME
                if (depth == 1 || depth == 2)
                    printf("Spatial Partitioning: %g[ms]\n", elapsed * 0.001f);
#endif
                
                uint32_t numLeftAdded, numRightAdded;
                uint32_t c0 = buildRecursive(fragments, currentSize + *numAdded, maximumBudget, start, splitIdx, depth, &numLeftAdded);
                uint32_t c1 = buildRecursive(fragments, currentSize + *numAdded + numLeftAdded, maximumBudget, splitIdx + numLeftAdded, end + *numAdded + numLeftAdded, depth, &numRightAdded);
                m_nodes[nodeIdx].initAsInternal(parentBB, c0, c1, widestAxisSP);
                *numAdded += numLeftAdded + numRightAdded;
                return nodeIdx;
            }
#undef PRINT_PROCESSING_TIME
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
        SBVH(const std::vector<SurfaceObject*> &objs) {
            std::chrono::system_clock::time_point tpStart, tpEnd;
            double elapsed;
            
            tpStart = std::chrono::system_clock::now();
            
            const uint32_t MemoryBudget = 5;
            Fragment* fragments = new Fragment[MemoryBudget * objs.size()];
            for (int i = 0; i < objs.size(); ++i) {
                BoundingBox3D bb = objs[i]->bounds();
                m_bounds.unify(bb);
                fragments[i].obj = objs[i];
                fragments[i].bbox = bb;
                fragments[i].costForIntersect = objs[i]->costForIntersect();
            }
            
            m_depth = 0;
            uint32_t numAdded;
            buildRecursive(fragments, (uint32_t)objs.size(), MemoryBudget * (uint32_t)objs.size(), 0, (uint32_t)objs.size(), 0, &numAdded);
            delete[] fragments;
            
            tpEnd = std::chrono::system_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(tpEnd - tpStart).count();
            
            m_cost = calcSAHCost();
//#ifdef DEBUG
            printf("num fragments: %u => %u, depth: %u, cost: %g, time: %g[s]\n", (uint32_t)objs.size(), (uint32_t)(objs.size() + numAdded), m_depth, m_cost, elapsed * 0.001f);
//#endif
        }
        
        float costForIntersect() const override {
            return m_cost;
        }
        
        BoundingBox3D bounds() const override {
            return m_bounds;
        }
        
        bool intersect(Ray &ray, Intersection* isect) const override {
            uint32_t objDepth = (uint32_t)isect->obj.size();
            bool dirIsPositive[] = {ray.dir.x >= 0, ray.dir.y >= 0, ray.dir.z >= 0};
            
            const uint32_t StackSize = 64;
            uint32_t idxStack[StackSize];
            uint32_t depth = 0;
            idxStack[depth++] = 0;
            while (depth > 0) {
                const Node &node = m_nodes[idxStack[--depth]];
                if (!node.bbox.intersect(ray))
                    continue;
                if (node.numLeaves == 0) {
                    SLRAssert(depth < StackSize, "SBVH::intersect: stack overflow");
                    bool positiveDir = dirIsPositive[node.axis];
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
    };
    
    inline __m128 _mm_sel_ps(const __m128 &mask, const __m128 &t, const __m128 &f) {
        // (((t ^ f) & mask)^f)
        return _mm_xor_ps(f, _mm_and_ps(mask, _mm_xor_ps(t, f)));
    }
    
    class QBVH : public Accelerator {
        struct Children {
            union {
                uint32_t asUInt;
                struct {
                    uint idx : 27;
                    uint numLeaves : 4;
                    bool isLeafNode : 1;
                };
            };
            
            bool isValid() const {
                return asUInt != UINT32_MAX;
            }
        };
        
        struct Node {
            __m128 min_x;
            __m128 min_y;
            __m128 min_z;
            __m128 max_x;
            __m128 max_y;
            __m128 max_z;
            Children children[4];
            BoundingBox3D::Axis topAxis;
            BoundingBox3D::Axis leftAxis;
            BoundingBox3D::Axis rightAxis;
            uint32_t padding;
            
            uint32_t intersect(const Ray &ray) const {
                const Vector3D invRayDir = ray.dir.reciprocal();
                
                __m128 rOrg_x = _mm_set_ps1(ray.org.x);
                __m128 rOrg_y = _mm_set_ps1(ray.org.y);
                __m128 rOrg_z = _mm_set_ps1(ray.org.z);
                __m128 invRayDir_x = _mm_set_ps1(invRayDir.x);
                __m128 invRayDir_y = _mm_set_ps1(invRayDir.y);
                __m128 invRayDir_z = _mm_set_ps1(invRayDir.z);
                
                __m128 tNear = _mm_set_ps1(ray.distMin);
                __m128 tFar = _mm_set_ps1(ray.distMax);
                
                tNear = _mm_max_ps(_mm_mul_ps(_mm_sub_ps(invRayDir.x > 0.0f ? min_x : max_x, rOrg_x), invRayDir_x), tNear);
                tNear = _mm_max_ps(_mm_mul_ps(_mm_sub_ps(invRayDir.y > 0.0f ? min_y : max_y, rOrg_y), invRayDir_y), tNear);
                tNear = _mm_max_ps(_mm_mul_ps(_mm_sub_ps(invRayDir.z > 0.0f ? min_z : max_z, rOrg_z), invRayDir_z), tNear);
                tFar = _mm_min_ps(_mm_mul_ps(_mm_sub_ps(invRayDir.x > 0.0f ? max_x : min_x, rOrg_x), invRayDir_x), tFar);
                tFar = _mm_min_ps(_mm_mul_ps(_mm_sub_ps(invRayDir.y > 0.0f ? max_y : min_y, rOrg_y), invRayDir_y), tFar);
                tFar = _mm_min_ps(_mm_mul_ps(_mm_sub_ps(invRayDir.z > 0.0f ? max_z : min_z, rOrg_z), invRayDir_z), tFar);
                
                return _mm_movemask_ps(_mm_cmple_ps(tNear, tFar));
            }
        };
        
        uint32_t m_depth;
        float m_cost;
        BoundingBox3D m_bounds;
        std::vector<Node> m_nodes;
        std::vector<const SurfaceObject*> m_objLists;
        
        Children collapseBBVH(const SBVH &baseBBVH, uint32_t grandparent, uint32_t depth) {
            Children ret;
            Children invalidChild = {UINT32_MAX};
            
            const SBVH::Node* root = &baseBBVH.m_nodes[grandparent];
            if (root->numLeaves > 0) {
                uint32_t baseIdx = (uint32_t)m_objLists.size();
                uint32_t numLeaves = root->numLeaves;
                SLRAssert(numLeaves <= 16, "The number of leaves for QBVH node is currently limited to 16.");
                for (int i = 0; i < numLeaves; ++i)
                    m_objLists.push_back(baseBBVH.m_objLists[root->offsetFirstLeaf + i]);
                ret.isLeafNode = true;
                ret.numLeaves = numLeaves;
                ret.idx = baseIdx;
                return ret;
            }
            
            if (++depth > m_depth)
                m_depth = depth;
            
            uint32_t childIdx[2] = {root->c0, root->c1};
            const SBVH::Node* child[2] = {&baseBBVH.m_nodes[childIdx[0]], &baseBBVH.m_nodes[childIdx[1]]};
            uint32_t grandchildIdx[4] = {UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX};
            const SBVH::Node* grandchild[4] = {nullptr, nullptr, nullptr, nullptr};
            if (child[0]->numLeaves == 0) {
                grandchildIdx[0] = child[0]->c0;
                grandchildIdx[1] = child[0]->c1;
                grandchild[0] = &baseBBVH.m_nodes[grandchildIdx[0]];
                grandchild[1] = &baseBBVH.m_nodes[grandchildIdx[1]];
            }
            if (child[1]->numLeaves == 0) {
                grandchildIdx[2] = child[1]->c0;
                grandchildIdx[3] = child[1]->c1;
                grandchild[2] = &baseBBVH.m_nodes[grandchildIdx[2]];
                grandchild[3] = &baseBBVH.m_nodes[grandchildIdx[3]];
            }
            
            uint32_t nodeIdx = (uint32_t)m_nodes.size();
            m_nodes.emplace_back();
            Node &node = m_nodes.back();
            node.leftAxis = node.rightAxis = BoundingBox3D::Axis_X;
            
            Children children[4] = {invalidChild, invalidChild, invalidChild, invalidChild};
            if (child[0]->numLeaves == 0) {
                if (child[1]->numLeaves == 0) {
                    node.min_x = _mm_setr_ps(grandchild[0]->bbox.minP.x, grandchild[1]->bbox.minP.x, grandchild[2]->bbox.minP.x, grandchild[3]->bbox.minP.x);
                    node.min_y = _mm_setr_ps(grandchild[0]->bbox.minP.y, grandchild[1]->bbox.minP.y, grandchild[2]->bbox.minP.y, grandchild[3]->bbox.minP.y);
                    node.min_z = _mm_setr_ps(grandchild[0]->bbox.minP.z, grandchild[1]->bbox.minP.z, grandchild[2]->bbox.minP.z, grandchild[3]->bbox.minP.z);
                    node.max_x = _mm_setr_ps(grandchild[0]->bbox.maxP.x, grandchild[1]->bbox.maxP.x, grandchild[2]->bbox.maxP.x, grandchild[3]->bbox.maxP.x);
                    node.max_y = _mm_setr_ps(grandchild[0]->bbox.maxP.y, grandchild[1]->bbox.maxP.y, grandchild[2]->bbox.maxP.y, grandchild[3]->bbox.maxP.y);
                    node.max_z = _mm_setr_ps(grandchild[0]->bbox.maxP.z, grandchild[1]->bbox.maxP.z, grandchild[2]->bbox.maxP.z, grandchild[3]->bbox.maxP.z);
                    
                    node.topAxis = root->axis;
                    node.leftAxis = child[0]->axis;
                    node.rightAxis = child[1]->axis;

                    children[0] = collapseBBVH(baseBBVH, grandchildIdx[0], depth);
                    children[1] = collapseBBVH(baseBBVH, grandchildIdx[1], depth);
                    children[2] = collapseBBVH(baseBBVH, grandchildIdx[2], depth);
                    children[3] = collapseBBVH(baseBBVH, grandchildIdx[3], depth);
                }
                else {
                    node.min_x = _mm_setr_ps(grandchild[0]->bbox.minP.x, grandchild[1]->bbox.minP.x, child[1]->bbox.minP.x, INFINITY);
                    node.min_y = _mm_setr_ps(grandchild[0]->bbox.minP.y, grandchild[1]->bbox.minP.y, child[1]->bbox.minP.y, INFINITY);
                    node.min_z = _mm_setr_ps(grandchild[0]->bbox.minP.z, grandchild[1]->bbox.minP.z, child[1]->bbox.minP.z, INFINITY);
                    node.max_x = _mm_setr_ps(grandchild[0]->bbox.maxP.x, grandchild[1]->bbox.maxP.x, child[1]->bbox.maxP.x, -INFINITY);
                    node.max_y = _mm_setr_ps(grandchild[0]->bbox.maxP.y, grandchild[1]->bbox.maxP.y, child[1]->bbox.maxP.y, -INFINITY);
                    node.max_z = _mm_setr_ps(grandchild[0]->bbox.maxP.z, grandchild[1]->bbox.maxP.z, child[1]->bbox.maxP.z, -INFINITY);
                    
                    node.topAxis = root->axis;
                    node.leftAxis = child[0]->axis;
                    
                    children[0] = collapseBBVH(baseBBVH, grandchildIdx[0], depth);
                    children[1] = collapseBBVH(baseBBVH, grandchildIdx[1], depth);
                    children[2] = collapseBBVH(baseBBVH, childIdx[1], depth);
                }
            }
            else {
                if (child[1]->numLeaves == 0) {
                    node.min_x = _mm_setr_ps(child[0]->bbox.minP.x, INFINITY, grandchild[2]->bbox.minP.x, grandchild[3]->bbox.minP.x);
                    node.min_y = _mm_setr_ps(child[0]->bbox.minP.y, INFINITY, grandchild[2]->bbox.minP.y, grandchild[3]->bbox.minP.y);
                    node.min_z = _mm_setr_ps(child[0]->bbox.minP.z, INFINITY, grandchild[2]->bbox.minP.z, grandchild[3]->bbox.minP.z);
                    node.max_x = _mm_setr_ps(child[0]->bbox.maxP.x, -INFINITY, grandchild[2]->bbox.maxP.x, grandchild[3]->bbox.maxP.x);
                    node.max_y = _mm_setr_ps(child[0]->bbox.maxP.y, -INFINITY, grandchild[2]->bbox.maxP.y, grandchild[3]->bbox.maxP.y);
                    node.max_z = _mm_setr_ps(child[0]->bbox.maxP.z, -INFINITY, grandchild[2]->bbox.maxP.z, grandchild[3]->bbox.maxP.z);
                    
                    node.topAxis = root->axis;
                    node.rightAxis = child[1]->axis;
                    
                    children[0] = collapseBBVH(baseBBVH, childIdx[0], depth);
                    children[2] = collapseBBVH(baseBBVH, grandchildIdx[2], depth);
                    children[3] = collapseBBVH(baseBBVH, grandchildIdx[3], depth);
                }
                else {
                    node.min_x = _mm_setr_ps(child[0]->bbox.minP.x, INFINITY, child[1]->bbox.minP.x, INFINITY);
                    node.min_y = _mm_setr_ps(child[0]->bbox.minP.y, INFINITY, child[1]->bbox.minP.y, INFINITY);
                    node.min_z = _mm_setr_ps(child[0]->bbox.minP.z, INFINITY, child[1]->bbox.minP.z, INFINITY);
                    node.max_x = _mm_setr_ps(child[0]->bbox.maxP.x, -INFINITY, child[1]->bbox.maxP.x, -INFINITY);
                    node.max_y = _mm_setr_ps(child[0]->bbox.maxP.y, -INFINITY, child[1]->bbox.maxP.y, -INFINITY);
                    node.max_z = _mm_setr_ps(child[0]->bbox.maxP.z, -INFINITY, child[1]->bbox.maxP.z, -INFINITY);
                    
                    node.topAxis = root->axis;
                    
                    children[0] = collapseBBVH(baseBBVH, childIdx[0], depth);
                    children[2] = collapseBBVH(baseBBVH, childIdx[1], depth);
                }
            }
            // do NOT use the variable "node" instead of "m_nodes[nodeIdx]".
            m_nodes[nodeIdx].children[0] = children[0];
            m_nodes[nodeIdx].children[1] = children[1];
            m_nodes[nodeIdx].children[2] = children[2];
            m_nodes[nodeIdx].children[3] = children[3];
            
            ret.isLeafNode = false;
            ret.numLeaves = 0;
            ret.idx = nodeIdx;
            return ret;
        }
        
        float calcSAHCost() const {
            const float Ci = 1.2f;
            float costInt = 0.0f;
            float costObj = 0.0f;
            for (int i = 0; i < m_nodes.size(); ++i) {
                const Node &node = m_nodes[i];
                BoundingBox3D nodeBB;
                BoundingBox3D cBBs[4];
                for (int c = 0; c < 4; ++c) {
                    if (!node.children[c].isValid())
                        continue;
                    const float* xMin = (const float*)&node.min_x;
                    const float* yMin = (const float*)&node.min_y;
                    const float* zMin = (const float*)&node.min_z;
                    const float* xMax = (const float*)&node.max_x;
                    const float* yMax = (const float*)&node.max_y;
                    const float* zMax = (const float*)&node.max_z;
                    cBBs[c].minP = Point3D(xMin[c], yMin[c], zMin[c]);
                    cBBs[c].maxP = Point3D(xMax[c], yMax[c], zMax[c]);
                    
                    nodeBB.unify(cBBs[c]);
                }
                
                float surfaceArea = nodeBB.surfaceArea();
                costInt += surfaceArea;
                SLRAssert(std::isfinite(surfaceArea), "invalid surface area value.");
                
                for (int c = 0; c < 4; ++c) {
                    Children child = node.children[c];
                    if (!child.isValid())
                        continue;
                    
                    float cSurfaceArea = cBBs[c].surfaceArea();
                    if (child.isLeafNode) {
                        float costPrims = 0.0f;
                        for (uint32_t j = 0; j < child.numLeaves; ++j)
                            costPrims += m_objLists[child.idx + j]->costForIntersect();
                        costObj += cSurfaceArea * costPrims;
                    }
                }
            }
            float rootSA = m_bounds.surfaceArea();
            costInt *= Ci / rootSA;
            costObj /= rootSA;
            
            return costInt + costObj;
        }
        
    public:
        QBVH(const SBVH &baseBBVH) {
            std::chrono::system_clock::time_point tpStart, tpEnd;
            double elapsed;
            
            tpStart = std::chrono::system_clock::now();
            
            m_depth = 0;
            m_bounds = baseBBVH.bounds();
            collapseBBVH(baseBBVH, 0, 0);
            
            tpEnd = std::chrono::system_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(tpEnd - tpStart).count();
            
            m_cost = calcSAHCost();
            printf("depth: %u, cost: %g, time: %g[s]\n", m_depth, m_cost, elapsed * 0.001f);
        }
        
        float costForIntersect() const override {
            return m_cost;
        }
        
        BoundingBox3D bounds() const override {
            return m_bounds;
        }
        
        bool intersect(Ray &ray, Intersection* isect) const override {
            uint32_t objDepth = (uint32_t)isect->obj.size();
            bool dirIsPositive[] = {ray.dir.x >= 0, ray.dir.y >= 0, ray.dir.z >= 0};
            
            const uint32_t StackSize = 64;
            uint32_t idxStack[StackSize];
            uint32_t depth = 0;
            idxStack[depth++] = 0;
            while (depth > 0) {
                const Node &node = m_nodes[idxStack[--depth]];
                uint32_t hitFlags = node.intersect(ray);
                if (hitFlags == 0)
                    continue;
                
                const uint32_t OrderTable[] = {
                    0x0123, 0x0132, 0x1023, 0x1032,
                    0x2301, 0x3201, 0x2310, 0x3210
                };
                uint32_t encodedOrder = OrderTable[4 * dirIsPositive[node.topAxis] + 2 * dirIsPositive[node.leftAxis] + 1 * dirIsPositive[node.rightAxis]];
                uint32_t order[4] = {(encodedOrder >> 0) & 0xF, (encodedOrder >> 4) & 0xF, (encodedOrder >> 8) & 0xF, (encodedOrder >> 12) & 0xF};
                Children children[] = {node.children[order[0]], node.children[order[1]], node.children[order[2]], node.children[order[3]]};
                for (int i = 0; i < 4; ++i) {
                    bool hit = ((hitFlags >> order[i]) & 0x1) == 0x1;
                    if (!hit)
                        children[i].asUInt = UINT32_MAX;
                }
                
                for (int i = 3; i >= 0; --i) {
                    const Children &child = children[i];
                    if (!child.isValid() || child.isLeafNode)
                        continue;
                    SLRAssert(depth < StackSize, "QBVH::intersect: stack overflow");
                    idxStack[depth++] = child.idx;
                }
                for (int i = 0; i < 4; ++i) {
                    const Children &child = children[i];
                    if (!child.isValid() || !child.isLeafNode)
                        continue;
                    for (uint32_t j = 0; j < child.numLeaves; ++j)
                        if (m_objLists[child.idx + j]->intersect(ray, isect))
                            ray.distMax = isect->dist;
                }
            }
            return isect->obj.size() > objDepth;
        }
    };
}

#endif /* SBVH_h */

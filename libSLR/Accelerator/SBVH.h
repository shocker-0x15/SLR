//
//  SBVH.h
//
//  Created by 渡部 心 on 2016/07/03.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef __SLR_SBVH__
#define __SLR_SBVH__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/accelerator.h"

namespace SLR {
    // References
    // Spatial Splits in Bounding Volume Hierarchies
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
            
            if ((pcBBMax - pcBBMin) > 0) {
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
                        if (!splitLeftBBox.isValid())
                            --numLeftIndices;
                        if (!splitRightBBox.isValid())
                            --numRightIndices;
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
                    SLRAssert(depth < StackSize, "SBVH::intersect: stack overflow");
                    bool positiveDir = dirIsPositive[node.axis];
                    idxStack[depth++] = positiveDir ? node.c1 : node.c0;
                    idxStack[depth++] = positiveDir ? node.c0 : node.c1;
                }
                else {
                    for (uint32_t i = 0; i < node.numLeaves; ++i) {
#ifdef DEBUG
                        if (Accelerator::traceTraverse) {
                            debugPrintf("%s%u\n",
                                        Accelerator::traceTraversePrefix.c_str(), node.offsetFirstLeaf + i);
                        }
#endif
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

#endif /* __SLR_SBVH__ */

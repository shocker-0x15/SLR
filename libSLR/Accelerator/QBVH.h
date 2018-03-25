//
//  QBVH.h
//
//  Created by 渡部 心 on 2016/08/14.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef __SLR_QBVH__
#define __SLR_QBVH__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/accelerator.h"
#include "../Accelerator/SBVH.h"
#include <nmmintrin.h>

namespace SLR {
    inline __m128 _mm_sel_ps(const __m128 &mask, const __m128 &t, const __m128 &f) {
        // (((t ^ f) & mask)^f)
        return _mm_xor_ps(f, _mm_and_ps(mask, _mm_xor_ps(t, f)));
    }
    
    // References
    // Shallow Bounding Volume Hierarchies for Fast SIMD Ray Tracing of Incoherent Rays
    class QBVH : public Accelerator {
        struct Children {
            union {
                uint32_t asUInt;
                struct {
                    unsigned int idx : 27;
                    unsigned int numLeaves : 4;
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
            
            uint32_t intersect(const Ray &ray, const RaySegment &segment) const {
                const Vector3D invRayDir = ray.dir.reciprocal();
                
                __m128 rOrg_x = _mm_set_ps1(ray.org.x);
                __m128 rOrg_y = _mm_set_ps1(ray.org.y);
                __m128 rOrg_z = _mm_set_ps1(ray.org.z);
                __m128 invRayDir_x = _mm_set_ps1(invRayDir.x);
                __m128 invRayDir_y = _mm_set_ps1(invRayDir.y);
                __m128 invRayDir_z = _mm_set_ps1(invRayDir.z);
                
                __m128 tNear = _mm_set_ps1(segment.distMin);
                __m128 tFar = _mm_set_ps1(segment.distMax);
                
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
            const Children invalidChild = {{UINT32_MAX}};
            
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
            Children rootResult = collapseBBVH(baseBBVH, 0, 0);
            if (rootResult.isLeafNode) {
                const Children invalidChild = {UINT32_MAX};
                
                m_nodes.emplace_back();
                Node &node = m_nodes.back();
                
                const SBVH::Node* orgRoot = &baseBBVH.m_nodes[0];
                node.min_x = _mm_setr_ps(orgRoot->bbox.minP.x, INFINITY, INFINITY, INFINITY);
                node.min_y = _mm_setr_ps(orgRoot->bbox.minP.y, INFINITY, INFINITY, INFINITY);
                node.min_z = _mm_setr_ps(orgRoot->bbox.minP.z, INFINITY, INFINITY, INFINITY);
                node.max_x = _mm_setr_ps(orgRoot->bbox.maxP.x, -INFINITY, -INFINITY, -INFINITY);
                node.max_y = _mm_setr_ps(orgRoot->bbox.maxP.y, -INFINITY, -INFINITY, -INFINITY);
                node.max_z = _mm_setr_ps(orgRoot->bbox.maxP.z, -INFINITY, -INFINITY, -INFINITY);
                
                node.children[1] = node.children[2] = node.children[3] = invalidChild;
                node.children[0] = rootResult;
            }
            
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
        
        bool intersectWithoutAlpha(const Ray &ray, const RaySegment &segment, SurfaceInteraction* si) const override {
            const SurfaceObject* closestObject = nullptr;
            bool dirIsPositive[] = {ray.dir.x >= 0, ray.dir.y >= 0, ray.dir.z >= 0};
            RaySegment isectRange = segment;
            
            const uint32_t StackSize = 64;
            uint32_t idxStack[StackSize];
            uint32_t depth = 0;
            idxStack[depth++] = 0;
            while (depth > 0) {
                const Node &node = m_nodes[idxStack[--depth]];
                uint32_t hitFlags = node.intersect(ray, isectRange);
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
                    for (uint32_t j = 0; j < child.numLeaves; ++j) {
                        if (m_objLists[child.idx + j]->intersectWithoutAlpha(ray, isectRange, si)) {
                            closestObject = m_objLists[child.idx + j];
                            isectRange.distMax = si->getDistance();
                        }
                    }
                }
            }
            return closestObject != nullptr;
        }
        
        bool intersect(const Ray &ray, const RaySegment &segment, LightPathSampler &pathSampler, SurfaceInteraction* si, const SurfaceObject** closestObject) const override {
            *closestObject = nullptr;
            bool dirIsPositive[] = {ray.dir.x >= 0, ray.dir.y >= 0, ray.dir.z >= 0};
            RaySegment isectRange = segment;
            
            const uint32_t StackSize = 64;
            uint32_t idxStack[StackSize];
            uint32_t depth = 0;
            idxStack[depth++] = 0;
            while (depth > 0) {
                const Node &node = m_nodes[idxStack[--depth]];
                uint32_t hitFlags = node.intersect(ray, isectRange);
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
                    for (uint32_t j = 0; j < child.numLeaves; ++j) {
                        if (m_objLists[child.idx + j]->intersect(ray, isectRange, pathSampler, si)) {
                            *closestObject = m_objLists[child.idx + j];
                            isectRange.distMax = si->getDistance();
                        }
                    }
                }
            }
            return *closestObject != nullptr;
        }
        
        float testVisibility(const Ray &ray, const RaySegment &segment) const override {
            float fractionalVisibility = 1.0f;
            bool dirIsPositive[] = {ray.dir.x >= 0, ray.dir.y >= 0, ray.dir.z >= 0};
            
            const uint32_t StackSize = 64;
            uint32_t idxStack[StackSize];
            uint32_t depth = 0;
            idxStack[depth++] = 0;
            while (depth > 0) {
                const Node &node = m_nodes[idxStack[--depth]];
                uint32_t hitFlags = node.intersect(ray, segment);
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
                    for (uint32_t j = 0; j < child.numLeaves; ++j) {
                        fractionalVisibility *= m_objLists[child.idx + j]->testVisibility(ray, segment);
                        if (fractionalVisibility == 0.0f)
                            return 0.0f;
                    }
                }
            }
            return fractionalVisibility;
        }
    };
}

#endif /* __SLR_QBVH__ */

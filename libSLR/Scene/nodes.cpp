//
//  nodes.cpp
//
//  Created by 渡部 心 on 2017/01/05.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "nodes.h"
#include "../Memory/ArenaAllocator.h"
#include "../Core/Transform.h"
#include "../Core/SurfaceObject.h"
#include "../SurfaceMaterials/IBLEmission.h"

namespace SLR {
    void InternalNode::addChildNode(Node *node) {
        m_childNodes.push_back(node);
    }
    
    void InternalNode::createRenderingData(Allocator* mem, const Transform* subTF, RenderingData *data) {
        if (subTF)
            m_appliedTransform = ChainedTransform(subTF, m_localToWorld).reduce(mem);// &m_localToWorld or m_localToWorld.copy(tfMem)?
        else
            m_appliedTransform = m_localToWorld->copy(mem);
        
        if (m_localToWorld->isStatic()) {
            for (int i = 0; i < m_childNodes.size(); ++i) {
                Node* child = m_childNodes[i];
                if (child->isDirectlyTransformable()) {
                    child->createRenderingData(mem, m_appliedTransform, data);
                }
                else {
                    RenderingData subData(nullptr);
                    child->createRenderingData(mem, nullptr, &subData);
                    m_TFObjs.push_back(nullptr);
                    if (subData.surfObjs.size() > 0) {
                        m_TFObjs.back() = mem->create<TransformedSurfaceObject>(subData.surfObjs[0], m_appliedTransform);
                        data->surfObjs.push_back(m_TFObjs.back());
                    }
                    
                    // TODO: need to check the program not to do more than once.
                    if (subData.camera) {
                        data->camera = subData.camera;
                        data->camTransform = m_appliedTransform;
                    }
                }
            }
        }
        else {
            RenderingData subData(nullptr);
            for (int i = 0; i < m_childNodes.size(); ++i)
                m_childNodes[i]->createRenderingData(mem, nullptr, &subData);
            
            if (subData.surfObjs.size() > 0) {
                SurfaceObject* child;
                if (subData.surfObjs.size() > 1) {
                    m_subObj = mem->create<SurfaceObjectAggregate>(subData.surfObjs);
                    child = m_subObj;
                }
                else if (subData.surfObjs.size() == 1) {
                    child = subData.surfObjs[0];
                }
                
                m_TFObjs.push_back(mem->create<TransformedSurfaceObject>(child, m_appliedTransform));
                data->surfObjs.push_back(m_TFObjs.back());
            }
            
            if (subData.camera) {
                data->camera = subData.camera;
                data->camTransform = m_appliedTransform;
            }
        }
    }
    
    void InternalNode::destroyRenderingData(Allocator* mem) {
        if (m_localToWorld->isStatic()) {
            for (int i = (int)m_childNodes.size() - 1; i >= 0; --i) {
                Node* child = m_childNodes[i];
                if (child->isDirectlyTransformable()) {
                    child->destroyRenderingData(mem);
                }
                else {
                    TransformedSurfaceObject* tfObj = m_TFObjs.back();
                    if (tfObj)
                        mem->destroy(tfObj);
                    m_TFObjs.pop_back();
                }
            }
            SLRAssert(m_TFObjs.size() == 0, "Destroying transformed objects is inconsitent.");
        }
        else {
            if (m_TFObjs.size() > 0) {
                mem->destroy(m_TFObjs.back());
                m_TFObjs.clear();
                
                if (m_subObj)
                    mem->destroy(m_subObj);
                m_subObj = nullptr;
            }
            
            for (int i = (int)m_childNodes.size() - 1; i >= 0; --i)
                m_childNodes[i]->destroyRenderingData(mem);
        }
        
        mem->destroy(m_appliedTransform);
        m_appliedTransform = nullptr;
    }
    
    
    
    void ReferenceNode::createRenderingData(Allocator* mem, const Transform *subTF, RenderingData *data) {
        if (!m_ready) {
            RenderingData subData(nullptr);
            m_node->createRenderingData(mem, nullptr, &subData);
            if (subData.surfObjs.size() > 1) {
                m_obj = mem->create<SurfaceObjectAggregate>(subData.surfObjs);
                m_isAggregate = true;
            }
            else {
                m_obj = subData.surfObjs[0];
                m_isAggregate = false;
            }
            m_ready = true;
        }
        
        data->surfObjs.push_back(m_obj);
    }
    
    void ReferenceNode::destroyRenderingData(Allocator* mem) {
        if (m_ready) {
            if (m_isAggregate)
                mem->destroy(m_obj);
            m_node->destroyRenderingData(mem);
            m_ready = false;
        }
    }
    
    
    
    void InfiniteSphereNode::createRenderingData(Allocator* mem, const Transform* subTF, RenderingData* data) {
        SLRAssert(subTF == nullptr, "Transformation to InfiniteSphereNode is currently not supported.");
        m_emission = mem->create<IBLEmission>(m_scene, m_IBLTex, m_scale);
        m_obj = mem->create<InfiniteSphereSurfaceObject>(data->scene, m_emission);
        data->envObj = m_obj;
    }
    
    void InfiniteSphereNode::destroyRenderingData(Allocator* mem) {
        if (m_obj)
            mem->destroy(m_obj);
        m_obj = nullptr;
        if (m_emission)
            mem->destroy(m_emission);
        m_emission = nullptr;
    }
}

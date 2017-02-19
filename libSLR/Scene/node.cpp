//
//  nodes.cpp
//
//  Created by 渡部 心 on 2017/01/05.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "node.h"
#include "../MemoryAllocators/ArenaAllocator.h"
#include "../Core/transform.h"
#include "../Core/surface_object.h"
#include "../Core/medium_object.h"
#include "../SurfaceMaterial/IBLEmitterSurfaceProperty.h"
#include "../SurfaceShape/InfinitesimalPointSurfaceShape.h"

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
                    m_TFSurfObjs.push_back(nullptr);
                    if (subData.surfObjs.size() > 0) {
                        m_TFSurfObjs.back() = mem->create<TransformedSurfaceObject>(subData.surfObjs[0], m_appliedTransform);
                        data->surfObjs.push_back(m_TFSurfObjs.back());
                    }
                    m_TFMedObjs.push_back(nullptr);
                    if (subData.medObjs.size() > 0) {
                        m_TFMedObjs.back() = mem->create<TransformedMediumObject>(subData.medObjs[0], m_appliedTransform);
                        data->medObjs.push_back(m_TFMedObjs.back());
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
                    m_subSurfObj = mem->create<SurfaceObjectAggregate>(subData.surfObjs);
                    child = m_subSurfObj;
                }
                else if (subData.surfObjs.size() == 1) {
                    child = subData.surfObjs[0];
                }
                
                m_TFSurfObjs.push_back(mem->create<TransformedSurfaceObject>(child, m_appliedTransform));
                data->surfObjs.push_back(m_TFSurfObjs.back());
            }
            
            if (subData.medObjs.size() > 0) {
                MediumObject* child;
                if (subData.surfObjs.size() > 1) {
                    m_subMedObj = mem->create<MediumObjectAggregate>(subData.medObjs);
                    child = m_subMedObj;
                }
                else if (subData.surfObjs.size() == 1) {
                    child = subData.medObjs[0];
                }
                
                m_TFMedObjs.push_back(mem->create<TransformedMediumObject>(child, m_appliedTransform));
                data->medObjs.push_back(m_TFMedObjs.back());
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
                    TransformedMediumObject* tfMedObj = m_TFMedObjs.back();
                    if (tfMedObj)
                        mem->destroy(tfMedObj);
                    m_TFMedObjs.pop_back();
                    
                    TransformedSurfaceObject* tfSurfObj = m_TFSurfObjs.back();
                    if (tfSurfObj)
                        mem->destroy(tfSurfObj);
                    m_TFSurfObjs.pop_back();
                }
            }
            SLRAssert(m_TFSurfObjs.size() == 0, "Destroying transformed objects is inconsitent.");
        }
        else {
            if (m_TFMedObjs.size() > 0) {
                mem->destroy(m_TFMedObjs.back());
                m_TFMedObjs.clear();
                
                if (m_subMedObj)
                    mem->destroy(m_subMedObj);
                m_subMedObj = nullptr;
            }
            
            if (m_TFSurfObjs.size() > 0) {
                mem->destroy(m_TFSurfObjs.back());
                m_TFSurfObjs.clear();
                
                if (m_subSurfObj)
                    mem->destroy(m_subSurfObj);
                m_subSurfObj = nullptr;
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
    
    
    
    void InfinitesimalPointNode::createRenderingData(Allocator* mem, const Transform* subTF, RenderingData *data) {
        // apply transform
        StaticTransform transform;
        if (subTF) {
            SLRAssert(subTF->isStatic(), "Transformation given to InfinitesimalPointNode must be static.");
            subTF->sample(0.0f, &transform);
        }
        Point3D p = transform * m_position;
        Vector3D d = transform * m_direction;
        
        m_surface = mem->create<InfinitesimalPointSurfaceShape>(p, d);
        m_obj = mem->create<SingleSurfaceObject>(m_surface, m_material);
        data->surfObjs.push_back(m_obj);
    }
    
    void InfinitesimalPointNode::destroyRenderingData(Allocator* mem) {
        if (m_obj)
            mem->destroy(m_obj);
        m_obj = nullptr;
        if (m_surface)
            mem->destroy(m_surface);
        m_surface = nullptr;
    }
    
    
    
    void InfiniteSphereNode::createRenderingData(Allocator* mem, const Transform* subTF, RenderingData* data) {
        SLRAssert(subTF == nullptr, "Transformation to InfiniteSphereNode is currently not supported.");
        m_emission = mem->create<IBLEmitterSurfaceProperty>(m_scene, m_IBLTex, m_scale);
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

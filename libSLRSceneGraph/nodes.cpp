//
//  nodes.cpp
//
//  Created by 渡部 心 on 2015/03/24.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "nodes.h"
#include <libSLR/Core/SurfaceObject.h>
#include <libSLR/Memory/ArenaAllocator.h>
#include <libSLR/Core/cameras.h>
#include <libSLR/Core/Transform.h>

namespace SLRSceneGraph {
    bool InternalNode::addChildNode(const NodeRef &node) {
        // Create a shared_ptr for passing it to contains() with a No-Op deleter so that "this" will not be deleted.
        NodeRef thisRef(this, [](void* ptr){});
        if (node->contains(thisRef) || this == node.get()) {
            printf("This causes a circular reference.\n");
            return false;
        }
        if (node->isInstanced()) {
            auto it = std::find(m_childNodes.begin(), m_childNodes.end(), node);
            if (it != m_childNodes.end()) {
                printf("Another instanced node already exists in this node.\n");
                return false;
            }
        }
        else {
            if (contains(node)) {
                printf("This node already has the given node.\n");
                return false;
            }
        }
        m_childNodes.push_back(node);
        return true;
    }
    
    const NodeRef &InternalNode::childNodeAt(int i) const {
        return m_childNodes[i];
    }
    
    NodeRef &InternalNode::childNodeAt(int i) {
        return m_childNodes[i];
    }
    
    
    void InternalNode::setTransform(const TransformRef &tf) {
        m_localToWorld = tf;
    }
    
    const TransformRef InternalNode::getTransform() const {
        return m_localToWorld;
    }
    
    
    bool InternalNode::contains(const NodeRef &obj) const {
        for (int i = 0; i < m_childNodes.size(); ++i) {
            if (m_childNodes[i] == obj || m_childNodes[i]->contains(obj))
                return true;
        }
        return false;
    }
    
    bool InternalNode::hasChildren() const {
        return m_childNodes.size() > 0;
    }
    
    
    void InternalNode::applyTransform() {
        TransformRef tf = getTransform();
        if (tf->isStatic()) {
            for (int i = 0; i < m_childNodes.size(); ++i)
                m_childNodes[i]->applyTransform(*(SLR::StaticTransform*)tf.get());
            m_localToWorld = createShared<SLR::StaticTransform>(SLR::Matrix4x4::Identity);
        }
        else {
            SLRAssert(false, "Non static transform cannot be applied.");
        }
    }
    
    void InternalNode::applyTransform(const SLR::StaticTransform &tf) {
        SLRAssert_NotImplemented();
        //    m_localToWorld = tf * m_localToWorld;
    }
    
    void InternalNode::applyTransformToLeaf(const SLR::StaticTransform &tf) {
        SLRAssert_NotImplemented();
        //    StaticTransform cumT = tf * m_localToWorld;
        //    for (int i = 0; i < m_childNodes.size(); ++i) {
        //        NodeRef &child = m_childNodes[i];
        //        if (child->hasChildren())
        //            child->applyTransformToLeaf(cumT);
        //        else
        //            child->applyTransform(cumT);
        //    }
        //    m_localToWorld = Matrix4x4::Identity;
    }
    
    void InternalNode::getRenderingData(SLR::ArenaAllocator &mem, const SLR::Transform* subTF, RenderingData *data) {
        const SLR::Transform* reduced;
        if (subTF)
            reduced = SLR::ChainedTransform(subTF, m_localToWorld.get()).reduce(mem);// &m_localToWorld or m_localToWorld.copy(tfMem)?
        else
            reduced = m_localToWorld.get();
        
        if (m_localToWorld->isStatic()) {
            for (int i = 0; i < m_childNodes.size(); ++i)
                m_childNodes[i]->getRenderingData(mem, reduced, data);
        }
        else {
            RenderingData subData;
            for (int i = 0; i < m_childNodes.size(); ++i)
                m_childNodes[i]->getRenderingData(mem, nullptr, &subData);
            if (subData.surfObjs.size() > 1) {
                SLR::SurfaceObjectAggregate* aggr = mem.create<SLR::SurfaceObjectAggregate>(subData.surfObjs);
                SLR::TransformedSurfaceObject* tfobj = mem.create<SLR::TransformedSurfaceObject>(aggr, reduced);
                data->surfObjs.push_back(tfobj);
            }
            else if (subData.surfObjs.size() == 1) {
                SLR::TransformedSurfaceObject* tfobj = mem.create<SLR::TransformedSurfaceObject>(subData.surfObjs[0], reduced);
                data->surfObjs.push_back(tfobj);
            }
            
            if (subData.camera) {
                data->camera = subData.camera;
                data->camTransform = SLR::ChainedTransform(reduced, subData.camTransform).reduce(mem);
            }
        }
    }
    
    
    SurfaceObjectNode::~SurfaceObjectNode() {
        if (!m_ready)
            return;
        for (int i = 0; i < m_numRefinedObjs; ++i)
            delete m_refinedObjs[i];
        delete[] m_refinedObjs;
    }
    
    void SurfaceObjectNode::getRenderingData(SLR::ArenaAllocator &mem, const SLR::Transform* subTF, RenderingData *data) {
        if (!m_ready) {
            createSurfaceObjects();
            m_ready = true;
        }
        size_t curSize = data->surfObjs.size();
        data->surfObjs.resize(curSize + m_numRefinedObjs);
        
        SLR::StaticTransform transform;
        if (subTF) {
            SLRAssert(subTF->isStatic(), "Transformation given to SurfaceObjectNode must be static.");
            subTF->sample(0.0f, &transform);
        }
        if (transform.isIdentity()) {
            for (int i = 0; i < m_numRefinedObjs; ++i) {
                SLR::SingleSurfaceObject* obj = m_refinedObjs[i];
                data->surfObjs[curSize + i] = obj;
            }
        }
        else {
            for (int i = 0; i < m_numRefinedObjs; ++i) {
                SLR::SingleSurfaceObject* obj = m_refinedObjs[i];
                SLR::TransformedSurfaceObject* tfObj = mem.create<SLR::TransformedSurfaceObject>(obj, subTF);
                data->surfObjs[curSize + i] = tfObj;
            }
        }
    }
    
    
    void ReferenceNode::getRenderingData(SLR::ArenaAllocator &mem, const SLR::Transform* subTF, RenderingData *data) {
        if (!m_ready) {
            m_node->getRenderingData(mem, nullptr, &m_subData);
            if (m_subData.surfObjs.size() > 1)
                m_surfObj = mem.create<SLR::SurfaceObjectAggregate>(m_subData.surfObjs);
            else
                m_surfObj = m_subData.surfObjs[0];
            m_ready = true;
        }
        data->surfObjs.push_back(mem.create<SLR::TransformedSurfaceObject>(m_surfObj, subTF));
    }
    
    
    CameraNode::~CameraNode() {
        if (!m_ready)
            return;
        delete m_camera;
    }
    
    void CameraNode::getRenderingData(SLR::ArenaAllocator &mem, const SLR::Transform* subTF, RenderingData *data) {
        if (!m_ready) {
            createCamera();
            m_ready = true;
        }
        data->camera = m_camera;
        data->camTransform = subTF;
    }
}

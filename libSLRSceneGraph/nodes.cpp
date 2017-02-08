//
//  nodes.cpp
//
//  Created by 渡部 心 on 2015/03/24.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "nodes.h"
#include <libSLR/Memory/ArenaAllocator.h>
#include <libSLR/Core/SurfaceObject.h>
#include <libSLR/Core/MediumObject.h>
#include <libSLR/Core/Transform.h>
#include <libSLR/Scene/nodes.h>

namespace SLRSceneGraph {
    Node::~Node() {
        if (m_rawData)
            free(m_rawData);
    }
    
    
    
    void InternalNode::allocateRawData() {
        m_rawData = (SLR::Node*)malloc(sizeof(SLR::InternalNode));
    }
    
    void InternalNode::setupRawData() {
        new (m_rawData) SLR::InternalNode(m_localToWorld.get());
        m_setup = true;
    }
    
    void InternalNode::terminateRawData() {
        SLR::InternalNode &raw = *(SLR::InternalNode*)m_rawData;
        if (m_setup)
            raw.~InternalNode();
        m_setup = false;
    }
    
    InternalNode::InternalNode(const TransformRef &localToWorld) :
    m_localToWorld(localToWorld) {
        allocateRawData();
    }
    
    bool InternalNode::addChildNode(const NodeRef &node) {
        // Create a shared_ptr for passing it to contains() with a No-Op deleter so that "this" will not be deleted.
        NodeRef thisRef(this, [](void* ptr){});
        if (node->contains(thisRef) || this == node.get()) {
            printf("This causes a circular reference.\n");
            return false;
        }
        if (node->isUniqueInTree()) {
            if (contains(node)) {
                printf("This node already has the given node.\n");
                return false;
            }
        }
        else {
            // multiple instances in the same internal node make no sense.
            auto it = std::find(m_childNodes.begin(), m_childNodes.end(), node);
            if (it != m_childNodes.end()) {
                printf("Another instanced node already exists in this node.\n");
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
    
    NodeRef InternalNode::copy() const {
        InternalNodeRef ret = createShared<InternalNode>(m_localToWorld);
        ret->m_name = m_name;
        for (int i = 0; i < m_childNodes.size(); ++i) {
            NodeRef c = m_childNodes[i]->copy();
            ret->m_childNodes.push_back(c);
        }
        ret->setupRawData();
        return ret;
    }
    
    void InternalNode::prepareForRendering() {
        terminateRawData();
        setupRawData();
        
        SLR::InternalNode &raw = *(SLR::InternalNode*)getRaw();
        for (int i = 0; i < m_childNodes.size(); ++i) {
            m_childNodes[i]->prepareForRendering();
            raw.addChildNode(m_childNodes[i]->getRaw());
        }
    }
    
    void InternalNode::propagateTransform() {
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
    
    
    
    void ReferenceNode::allocateRawData() {
        m_rawData = (SLR::Node*)malloc(sizeof(SLR::ReferenceNode));
    }
    
    void ReferenceNode::setupRawData() {
        new (m_rawData) SLR::ReferenceNode(m_node->getRaw());
        m_setup = true;
    }
    
    void ReferenceNode::terminateRawData() {
        SLR::ReferenceNode &raw = *(SLR::ReferenceNode*)m_rawData;
        if (m_setup)
            raw.~ReferenceNode();
        m_setup = false;
    }
    
    ReferenceNode::ReferenceNode(const NodeRef &node) : m_node(node) {
        allocateRawData();
    }
    
    NodeRef ReferenceNode::copy() const {
        NodeRef ret = createShared<ReferenceNode>(m_node->copy());
        return ret;
    }
    
    void ReferenceNode::prepareForRendering() {
        terminateRawData();
        setupRawData();
        m_node->prepareForRendering();
    }
}

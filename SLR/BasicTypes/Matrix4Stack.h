//
//  Matrix4x4Stack.h
//
//  Created by 渡部 心 on 11/11/25.
//  Copyright (c) 2014年 渡部 心. All rights reserved.
//

#ifndef SLR_Matrix4fStack_h
#define SLR_Matrix4fStack_h

#include "../defines.h"
#include "../references.h"
#include "Matrix4f.h"
#include <stack>

class Matrix4x4Stack {
    std::stack<Matrix4x4> m_stack;
    bool m_mulFromLeft;
    
public:
    Matrix4x4Stack(bool mulFromLeft = false) {
        m_stack.push(Matrix4x4::Identity);
        m_mulFromLeft = mulFromLeft;
    };
    
    Matrix4x4Stack(const Matrix4x4& mat) {
        m_stack.push(mat);
    };
    
    void push() {
        m_stack.push(m_stack.top());
    };
    
    void pop() {
        if (m_stack.size() > 1)
            m_stack.pop();
        else
            m_stack.top() = Matrix4x4::Identity;
    };
    
    Matrix4x4Stack& operator=(const Matrix4x4& mat) {
        m_stack.top() = mat;
        return *this;
    };
    
    Matrix4x4Stack& operator*=(const Matrix4x4& mat) {
        if (m_mulFromLeft)
            m_stack.top() = mat * m_stack.top();
        else
            m_stack.top() *= mat;
        return *this;
    };
    
    const Matrix4x4& top() const {
        return m_stack.top();
    };
    
    Matrix4x4& top() {
        return m_stack.top();
    };
    
    operator Matrix4x4() const {
        return m_stack.top();
    };
};

#endif

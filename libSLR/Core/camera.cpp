//
//  Camera.cpp
//
//  Created by 渡部 心 on 2015/05/30.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "camera.h"

namespace SLR {
    void Camera::setTransform(const Transform *t) {
        m_transform = t;
    }
}

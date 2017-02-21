//
//  InfinitesimalPointSurfaceShape.cpp
//
//  Created by 渡部 心 on 2017/02/17.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "InfinitesimalPointSurfaceShape.h"

#include "../BasicTypes/CompensatedSum.h"
#include "../Core/distributions.h"
#include "../Core/surface_object.h"
#include "../Core/textures.h"

namespace SLR {
    void InfinitesimalPointSurfaceShape::calculateSurfacePoint(const SurfaceInteraction &si, SurfacePoint* surfPt) const {
        ReferenceFrame shadingFrame;
        shadingFrame.z = m_direction;
        shadingFrame.z.makeCoordinateSystem(&shadingFrame.x, &shadingFrame.y);
        
        *surfPt = SurfacePoint(si, false, shadingFrame, shadingFrame.x);
    }
    
    void InfinitesimalPointSurfaceShape::sample(float u0, float u1, SurfacePoint* surfPt, float* areaPDF, DirectionType* posType) const {
        ReferenceFrame shadingFrame;
        shadingFrame.z = m_direction;
        shadingFrame.z.makeCoordinateSystem(&shadingFrame.x, &shadingFrame.y);
        
        *surfPt = SurfacePoint(m_position,
                               false,
                               shadingFrame,
                               m_direction,
                               0.0f, 0.0f,
                               TexCoord2D(0.0f, 0.0f),
                               shadingFrame.x
                               );
        *areaPDF = 1.0f;
        *posType = DirectionType::Delta0D;
    }
}

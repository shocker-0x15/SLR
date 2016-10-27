//
//  MediumObject.h
//
//  Created by 渡部 心 on 2016/09/19.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef __SLR_MediumObject_h__
#define __SLR_MediumObject_h__

#include "../defines.h"
#include "../references.h"
#include "Object.h"

namespace SLR {
    class SLR_API MediumObject : public Object {
        
    };
    
    class SLR_API HomogeneousMediumObject : public MediumObject {
        
    };
    
    class SLR_API SpectrallyUniformExtinctionGridMediumObject : public MediumObject {
        
    };
    
    class SLR_API GridMediumObject : public MediumObject {
        
    };
}

#endif /* MediumObject_h */

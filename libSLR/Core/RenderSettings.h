//
//  RenderSettings.h
//
//  Created by 渡部 心 on 2015/08/01.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_RenderSettings__
#define __SLR_RenderSettings__

#include "../defines.h"
#include "../references.h"

namespace SLR {
    enum class RenderSettingItem {
        NumThreads,
        ImageWidth,
        ImageHeight,
        TimeStart,
        TimeEnd,
        Brightness,
        RNGSeed,
    };
    
    class SLR_API RenderSettings {
        std::map<RenderSettingItem, bool> m_boolValues;
        std::map<RenderSettingItem, int32_t> m_int32Values;
        std::map<RenderSettingItem, float> m_floatValues;
        std::map<RenderSettingItem, std::string> m_stringValues;
        std::map<RenderSettingItem, void*> m_pointerValue;
    public:
        void addItem(RenderSettingItem item, bool value) { m_boolValues[item] = value; };
        void addItem(RenderSettingItem item, int32_t value) { m_int32Values[item] = value; };
        void addItem(RenderSettingItem item, float value) { m_floatValues[item] = value; };
        void addItem(RenderSettingItem item, const std::string &value) { m_stringValues[item] = value; };
        void addItem(RenderSettingItem item, void* value) { m_pointerValue[item] = value; };
        
        bool getBool(RenderSettingItem item) const { return m_boolValues.at(item); };
        int32_t getInt(RenderSettingItem item) const { return m_int32Values.at(item); };
        float getFloat(RenderSettingItem item) const { return m_floatValues.at(item); };
        std::string getString(RenderSettingItem item) const { return m_stringValues.at(item); };
        void* getPointer(RenderSettingItem item) const { return m_pointerValue.at(item); };
    };    
}

#endif /* __SLR_RenderSettings__ */

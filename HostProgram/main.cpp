//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#define _USE_MATH_DEFINES
#include <cstdio>
#include <thread>

#include <libSLR/defines.h>
#include <libSLR/MemoryAllocators/ArenaAllocator.h>
#include <libSLR/BasicTypes/spectrum_base.h>
#include <libSLR/Core/renderer.h>
#include <libSLR/Core/RenderSettings.h>
#include <libSLR/Scene/Scene.h>
#include <libSLRSceneGraph/declarations.h>
#include <libSLRSceneGraph/Scene/Scene.h>
#include <libSLRSceneGraph/API.h>

#include "StopWatch.h"

int main(int argc, const char * argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Too few command line arguments.\n");
        return -1;
    }
    
    // print launching time
    using namespace std::chrono;
    std::time_t ctimeLaunch = system_clock::to_time_t(system_clock::now());
    slrprintf("%s\n", std::ctime(&ctimeLaunch));
    
    SLR::initializeColorSystem();
    
    StopWatch stopwatch;
    
    // read a scene
    stopwatch.start();
    SLRSceneGraph::SceneRef scene = createShared<SLRSceneGraph::Scene>();
    SLRSceneGraph::RenderingContext context;
    bool readSceneSuccess = SLRSceneGraph::readScene(argv[1], scene, &context);
    if (!readSceneSuccess) {
        slrprintf("Failed to read a scene file.\n");
        exit(-1);
    }
    slrprintf("read scene: %g [s]\n", stopwatch.stop() * 1e-3f);
    
    // setup render settings
    SLR::RenderSettings settings;
#ifdef DEBUG
    settings.addItem(SLR::RenderSettingItem::NumThreads, 1);
#else
    settings.addItem(SLR::RenderSettingItem::NumThreads, (int32_t)std::thread::hardware_concurrency());
#endif
    settings.addItem(SLR::RenderSettingItem::ImageWidth, context.width);
    settings.addItem(SLR::RenderSettingItem::ImageHeight, context.height);
    settings.addItem(SLR::RenderSettingItem::TimeStart, context.timeStart);
    settings.addItem(SLR::RenderSettingItem::TimeEnd, context.timeEnd);
    settings.addItem(SLR::RenderSettingItem::Brightness, context.brightness);
    settings.addItem(SLR::RenderSettingItem::RNGSeed, context.rngSeed);
    
    scene->prepareForRendering();
    SLR::Scene* rawScene = scene->getRaw();
    SLR::ArenaAllocator sceneMem;
    rawScene->build(&sceneMem);
    context.renderer->render(*rawScene, settings);
    rawScene->destory();
    
    return 0;
}

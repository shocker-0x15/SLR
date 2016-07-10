//
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include <cstdio>

#include <libSLR/defines.h>
#include <libSLRSceneGraph/references.h>
#include <libSLR/BasicTypes/Spectrum.h>
#include <libSLR/Memory/ArenaAllocator.h>
#include <libSLR/Core/Renderer.h>
#include <libSLR/Core/RenderSettings.h>
#include <libSLRSceneGraph/Scene.h>
#include <libSLRSceneGraph/API.hpp>

#include "StopWatch.h"

int main(int argc, const char * argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Too few command line arguments.\n");
        return -1;
    }
    
    using namespace std::chrono;
    SLR::initSpectrum();
    
    StopWatch stopwatch;
    StopWatchHiRes stopwatchHiRes;
    
    std::time_t ctimeLaunch = system_clock::to_time_t(stopwatch.start());
    printf("%s\n", std::ctime(&ctimeLaunch));
    
    stopwatch.start();
    SLRSceneGraph::SceneRef scene = createShared<SLRSceneGraph::Scene>();
    SLRSceneGraph::RenderingContext context;
    bool readSceneSuccess = SLRSceneGraph::readScene(argv[1], scene, &context);
    if (!readSceneSuccess) {
        printf("Failed to read a scene file.\n");
        exit(-1);
    }
    printf("read scene: %g [s]\n", stopwatch.stop() * 1e-3f);
    
    stopwatch.start();
    const SLR::Scene* rawScene;
    SLR::ArenaAllocator mem;
    scene->build(&rawScene, mem);
    printf("build scene: %g [s]\n", stopwatch.stop() * 1e-3f);
    
    SLR::RenderSettings settings;
    settings.addItem(SLR::RenderSettingItem::ImageWidth, context.width);
    settings.addItem(SLR::RenderSettingItem::ImageHeight, context.height);
    settings.addItem(SLR::RenderSettingItem::TimeStart, context.timeStart);
    settings.addItem(SLR::RenderSettingItem::TimeEnd, context.timeEnd);
    settings.addItem(SLR::RenderSettingItem::Brightness, context.brightness);
    settings.addItem(SLR::RenderSettingItem::RNGSeed, context.rngSeed);
    
    context.renderer->render(*rawScene, settings);
    
    return 0;
}

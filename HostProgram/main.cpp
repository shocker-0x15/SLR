//
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include <cstdio>

#include <libSLR/defines.h>
#include <libSLRSceneGraph/references.h>
#include <libSLR/BasicTypes/Spectrum.h>
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
    SLRSceneGraph::Scene scene;
    bool readSceneSuccess = SLRSceneGraph::readScene(argv[1], &scene);
    SLRAssert(readSceneSuccess, "Failed to read a scene file.");
    printf("build: %g [s]\n", stopwatch.stop() * 1e-3f);
    
    return 0;
}

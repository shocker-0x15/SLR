//
//  main.cpp
//
//  Created by 渡部 心 on 2017/05/22.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include <gtest/gtest.h>
#include <cstdio>
#include <iostream>

#include <libSLR/BasicTypes/spectrum_base.h>

class CustomPrinter : public ::testing::EmptyTestEventListener {
    
};

int main(int argc, char * argv[]) {
    SLR::initializeColorSystem();
    ::testing::InitGoogleTest(&argc, argv);
    
//    // Gets hold of the event listener list.
//    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
//    delete listeners.Release(listeners.default_result_printer());
//    // Adds a listener to the end.  Google Test takes the ownership.
//    listeners.Append(new CustomPrinter());
    
    return RUN_ALL_TESTS();
}

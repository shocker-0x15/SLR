//
//  ProgressReporter.cpp
//
//  Created by 渡部 心 on 2016/09/09.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "ProgressReporter.h"

namespace SLR {
    ProgressReporter::ProgressReporter() : m_sleep(false), m_numLastLines(0), m_finishable(false) {
        m_printThread = std::thread([this]() {
            const std::chrono::milliseconds sleepDuration(250);
            while (!m_finishable) {
                std::this_thread::sleep_for(sleepDuration);
                
                {
                    std::unique_lock<std::mutex> lock(m_jobMutex);
                    
                    if (m_sleep) {
                        m_printCondVar.wait(lock, [this]() {
                            return !m_sleep;
                        });
                    }
                
                    const size_t entireLength = 64;
                    char buffer[2 + entireLength];
                    
                    if (m_numLastLines > 0)
                        slrprintf("\033[%uF", m_numLastLines);
                    
                    for (int i = (int32_t)m_jobStask.size() - 1; i >= 0; --i) {
                        Job &job = m_jobStask[i];
                        
                        size_t headLength = job.title.length() + 1;
                        const size_t tailLength = 1 + 6;
                        size_t barLength = entireLength - (headLength + tailLength);
                        
                        snprintf(buffer, headLength + 1, "%s|", job.title.c_str());
                        char* curPtr = buffer + headLength;
                        
                        float percentage = (float)job.numWorksDone / job.totalWork;
                        uint32_t numFills = uint32_t(barLength * percentage);
                        for (int i = 0; i < numFills; ++i)
                            *(curPtr++) = '+';
                        for (int i = 0; i < barLength - numFills; ++i)
                            *(curPtr++) = '-';
                        snprintf(curPtr, tailLength + 1, "|%5.1f%%", percentage * 100);
                        curPtr += tailLength;
                        *curPtr = '\0';
                        
                        bool printFinished = job.numWorksDone / job.totalWork;
                        if (job.printFinished == false && printFinished) {
                            job.printFinished = true;
                            m_printCondVar.notify_one();
                        }
                        
                        slrprintf("%s\n", buffer);
                    }
                    m_numLastLines = (uint32_t)m_jobStask.size();
                    
                    fflush(stdout);
                }
            }
            slrprintf("\n");
        });
    }
    
    void ProgressReporter::pushJob(const std::string &title, uint64_t totalWork, const std::chrono::system_clock::time_point &startTime) {
        std::lock_guard<std::mutex> lock(m_jobMutex);
        m_jobStask.emplace_back(title, totalWork, startTime);
    }
    
    void ProgressReporter::popJob() {
        //std::lock_guard<std::mutex> lock(m_jobMutex);
        //m_jobStask.pop_back();
        std::unique_lock<std::mutex> lock(m_jobMutex);
        m_printCondVar.wait(lock, [this]() {
            return m_jobStask.back().printFinished;
        });
        m_jobStask.pop_back();
    }
    
    void ProgressReporter::beginOtherThreadPrint() {
        std::lock_guard<std::mutex> lock(m_jobMutex);
        slrprintf("\033[%uF", m_numLastLines);
        for (int i = 0; i < m_numLastLines; ++i)
            slrprintf("                                "
                      "                                \n"); // 64 spaces
        slrprintf("\033[%uF", m_numLastLines);
        m_sleep = true;
    }
    
    void ProgressReporter::endOtherThreadPrint() {
        std::lock_guard<std::mutex> lock(m_jobMutex);
        m_numLastLines = 0;
        m_sleep = false;
        m_printCondVar.notify_one();
    }
    
    void ProgressReporter::update(uint64_t numWorks) {
        for (Job &job : m_jobStask)
            job.numWorksDone += numWorks;
    }
    
    void ProgressReporter::finish() {
        m_finishable = true;
        m_printThread.join();
    }
}

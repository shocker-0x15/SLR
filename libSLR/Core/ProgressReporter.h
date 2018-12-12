//
//  ProgressReporter.h
//
//  Created by 渡部 心 on 2016/09/09.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef __SLR_ProgressReporter__
#define __SLR_ProgressReporter__

#include "../defines.h"
#include "../declarations.h"
#include <atomic>
#include <thread>
#include <condition_variable>
#include <mutex>

namespace SLR {
    class SLR_API ProgressReporter {
        struct Job {
            std::string title;
            uint64_t totalWork;
            std::chrono::system_clock::time_point startTime;
            std::atomic<uint64_t> numWorksDone;
            bool printFinished;
            
            Job(const Job& job) : title(job.title), totalWork(job.totalWork), startTime(job.startTime), numWorksDone(0), printFinished(false) {
                numWorksDone += job.numWorksDone;
            }
            Job(const std::string &title_, uint64_t totalWork_, const std::chrono::system_clock::time_point &startTime_ = std::chrono::system_clock::now()) :
            title(title_), totalWork(totalWork_), startTime(startTime_), numWorksDone(0), printFinished(false) { }
        };
        
        std::mutex m_jobMutex;
        std::vector<Job> m_jobStask;
        std::condition_variable m_printCondVar;
        std::thread m_printThread;
        bool m_sleep;
        uint32_t m_numLastLines;
        bool m_finishable;
    public:
        ProgressReporter();
        ~ProgressReporter() {}

        void pushJob(const std::string &title, uint64_t totalWork, const std::chrono::system_clock::time_point &startTime = std::chrono::system_clock::now());
        void popJob();
        void beginOtherThreadPrint();
        void endOtherThreadPrint();
        void update(uint64_t numWorks = 1);
        void finish();
        std::chrono::system_clock::duration elapsed() const {
            return std::chrono::system_clock::now() - m_jobStask.front().startTime;
        }
    };
}

#endif /* __SLR_ProgressReporter__ */

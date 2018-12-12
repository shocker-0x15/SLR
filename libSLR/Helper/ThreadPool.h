//
//  ThreadPool.h
//
//  Created by 渡部 心 on 2015/07/30.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_ThreadPool__
#define __SLR_ThreadPool__

#include "../defines.h"
#include "../declarations.h"
#include <thread>
#include <condition_variable>
#include <mutex>

class ThreadPool {
    typedef std::function<void(uint32_t threadID)> JobFunctionObject;
    
    class Worker {
        uint32_t m_threadID;
        ThreadPool& m_pool;
    public:
        Worker(uint32_t threadID, ThreadPool &pool) : m_threadID(threadID), m_pool(pool) { };
        
        void job() {
            JobFunctionObject task;
            while (true) {
                {
                    std::unique_lock<std::mutex> lock(m_pool.m_mutex);
                    while (!m_pool.m_finishable && m_pool.m_taskQueue.empty())
                        m_pool.m_condVar.wait(lock);
                    
                    if (m_pool.m_finishable && m_pool.m_taskQueue.empty())
                        return;
                    
                    task = m_pool.m_taskQueue.front();
                    m_pool.m_taskQueue.pop_front();
                }
                
                task(m_threadID);
            }
        }
    };
    
    friend class Worker;
    std::vector<std::thread> m_workers;
    std::deque<JobFunctionObject> m_taskQueue;
    std::mutex m_mutex;
    std::condition_variable m_condVar;
    bool m_finishable;
public:
    ThreadPool(uint32_t numThreads = std::thread::hardware_concurrency()) : m_finishable(false) {
        for (int i = 0; i < numThreads; ++i) {
            m_workers.push_back(std::thread(std::bind(&Worker::job, Worker(i, *this))));
        }
    };
    ~ThreadPool() {
        m_finishable = true;
        m_condVar.notify_all();
        for (int i = 0; i < m_workers.size(); ++i)
            if (m_workers[i].joinable())
                m_workers[i].join();
    };
    
    void enqueue(const JobFunctionObject &task) {
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_taskQueue.push_back(task);
        }
        m_condVar.notify_all();
    };
    
    void wait() {
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_finishable = true;
        }
        m_condVar.notify_all();
        for (int i = 0; i < m_workers.size(); ++i)
            m_workers[i].join();
    };
    
    uint32_t numThreads() const { return (uint32_t)m_workers.size(); };
};

#endif /* __SLR_ThreadPool__ */

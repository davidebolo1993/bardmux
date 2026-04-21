#pragma once

#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <future>
#include <atomic>

class ThreadPool {
public:
    explicit ThreadPool(int n_threads);
    ~ThreadPool();

    void submit(std::function<void()> task);
    void wait_all();   // block until all submitted tasks have completed

private:
    std::vector<std::thread>          workers_;
    std::queue<std::function<void()>> queue_;
    std::mutex                        mtx_;
    std::condition_variable           cv_task_;
    std::condition_variable           cv_done_;
    std::atomic<int>                  active_{0};
    bool                              stop_{false};
};

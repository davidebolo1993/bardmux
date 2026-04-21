#pragma once
#include <functional>
#include <thread>
#include <vector>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <atomic>

class ThreadPool {
public:
    // max_queue: max number of pending tasks before submit() blocks the producer.
    // Set to 2*num_threads so workers are always fed but RAM is bounded.
    explicit ThreadPool(int num_threads, int max_queue = 0);
    ~ThreadPool();

    // Blocks if the queue is full (backpressure).
    void submit(std::function<void()> task);

    // Wait for all submitted tasks to finish.
    void wait_all();

private:
    std::vector<std::thread>          workers_;
    std::queue<std::function<void()>> tasks_;

    std::mutex              mtx_;
    std::condition_variable cv_task_;    // workers wait here
    std::condition_variable cv_space_;   // producer waits here when queue full
    std::condition_variable cv_done_;    // wait_all() waits here

    std::atomic<int> active_{0};   // tasks currently running
    int  max_queue_;
    bool stop_ = false;
};

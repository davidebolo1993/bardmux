#include "thread_pool.h"

ThreadPool::ThreadPool(int n) {
    for (int i = 0; i < n; ++i) {
        workers_.emplace_back([this] {
            for (;;) {
                std::function<void()> task;
                {
                    std::unique_lock<std::mutex> lock(mtx_);
                    cv_task_.wait(lock, [this]{ return stop_ || !queue_.empty(); });
                    if (stop_ && queue_.empty()) return;
                    task = std::move(queue_.front());
                    queue_.pop();
                    ++active_;
                }
                task();
                {
                    std::unique_lock<std::mutex> lock(mtx_);
                    --active_;
                }
                cv_done_.notify_all();
            }
        });
    }
}

ThreadPool::~ThreadPool() {
    { std::unique_lock<std::mutex> lock(mtx_); stop_ = true; }
    cv_task_.notify_all();
    for (auto& w : workers_) w.join();
}

void ThreadPool::submit(std::function<void()> task) {
    { std::unique_lock<std::mutex> lock(mtx_); queue_.push(std::move(task)); }
    cv_task_.notify_one();
}

void ThreadPool::wait_all() {
    std::unique_lock<std::mutex> lock(mtx_);
    cv_done_.wait(lock, [this]{ return queue_.empty() && active_.load() == 0; });
}

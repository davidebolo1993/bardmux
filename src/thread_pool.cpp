#include "thread_pool.h"

ThreadPool::ThreadPool(int num_threads, int max_queue)
    : max_queue_(max_queue > 0 ? max_queue : num_threads * 2)
{
    for (int i = 0; i < num_threads; ++i) {
        workers_.emplace_back([this] {
            for (;;) {
                std::function<void()> task;
                {
                    std::unique_lock<std::mutex> lock(mtx_);
                    cv_task_.wait(lock, [this] {
                        return stop_ || !tasks_.empty();
                    });
                    if (stop_ && tasks_.empty()) return;
                    task = std::move(tasks_.front());
                    tasks_.pop();
                    ++active_;
                }
                // Notify producer: a slot freed up in the queue
                cv_space_.notify_one();

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
    {
        std::unique_lock<std::mutex> lock(mtx_);
        stop_ = true;
    }
    cv_task_.notify_all();
    for (auto& w : workers_) w.join();
}

void ThreadPool::submit(std::function<void()> task) {
    std::unique_lock<std::mutex> lock(mtx_);
    // Block producer until there is space in the queue
    cv_space_.wait(lock, [this] {
        return stop_ || static_cast<int>(tasks_.size()) < max_queue_;
    });
    if (stop_) return;
    tasks_.push(std::move(task));
    cv_task_.notify_one();
}

void ThreadPool::wait_all() {
    std::unique_lock<std::mutex> lock(mtx_);
    cv_done_.wait(lock, [this] {
        return tasks_.empty() && active_.load() == 0;
    });
}

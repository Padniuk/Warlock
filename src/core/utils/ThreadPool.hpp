/**
 * @file ThreadPool.hpp
 * @brief A persistent worker-thread pool shared by the whole pipeline.
 */

#ifndef WARLOCK_THREADPOOL_HPP
#define WARLOCK_THREADPOOL_HPP

#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <future>
#include <type_traits>

namespace framework {

    /**
     * @brief Fixed-size pool of worker threads processing a shared task
     * queue, constructed once by FrameworkManager and handed to every
     * Module. A module submits chunks of its own per-batch work via
     * submit(); tasks from different modules (including concurrently
     * running ones, see Module::runsConcurrently()) share the same queue
     * and worker threads rather than each getting a dedicated set.
     */
    class ThreadPool {
    public:
        /// @param threads Number of worker threads to start.
        explicit ThreadPool(size_t threads);
        /// Signals all workers to stop and joins them.
        ~ThreadPool();

        /**
         * @brief Enqueues `f(args...)` to run on a worker thread.
         * @return A future for the call's result.
         * @throws std::runtime_error if called after the pool has been
         * (or is being) destroyed.
         */
        template<class F, class... Args>
        auto submit(F&& f, Args&&... args) -> std::future<std::invoke_result_t<F, Args...>> {
            using return_type = std::invoke_result_t<F, Args...>;

            auto task = std::make_shared<std::packaged_task<return_type()>>(
                std::bind(std::forward<F>(f), std::forward<Args>(args)...)
            );

            std::future<return_type> res = task->get_future();
            {
                std::unique_lock<std::mutex> lock(queue_mutex_);
                if(stop_) throw std::runtime_error("enqueue on stopped ThreadPool");
                tasks_.emplace([task](){ (*task)(); });
            }
            condition_.notify_one();
            return res;
        }

        /// @return The number of worker threads in this pool.
        size_t getThreadCount() const { return workers_.size(); }

    private:
        std::vector<std::thread> workers_;
        std::queue<std::function<void()>> tasks_;
        std::mutex queue_mutex_;
        std::condition_variable condition_;
        bool stop_;
    };
}
#endif

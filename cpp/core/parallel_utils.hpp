#pragma once

#include <algorithm>
#include <cstddef>
#include <future>
#include <thread>
#include <utility>
#include <vector>

namespace xpongecpp {

inline std::size_t automatic_thread_count(std::size_t task_count, std::size_t min_tasks_per_thread = 1) {
    if (task_count <= 1 || min_tasks_per_thread == 0) {
        return 1;
    }
    const std::size_t hardware = std::max<std::size_t>(1, std::thread::hardware_concurrency());
    const std::size_t bounded = std::max<std::size_t>(1, task_count / min_tasks_per_thread);
    return std::max<std::size_t>(1, std::min(hardware, bounded));
}

template <class Func>
inline void parallel_for_chunks(std::size_t task_count, std::size_t min_tasks_per_thread, Func&& func) {
    const std::size_t thread_count = automatic_thread_count(task_count, min_tasks_per_thread);
    if (thread_count <= 1 || task_count <= 1) {
        if (task_count != 0) {
            func(0, task_count);
        }
        return;
    }

    const std::size_t chunk_size = (task_count + thread_count - 1) / thread_count;
    std::vector<std::future<void>> workers;
    workers.reserve(thread_count);
    for (std::size_t begin = 0; begin < task_count; begin += chunk_size) {
        const std::size_t end = std::min(task_count, begin + chunk_size);
        workers.push_back(std::async(std::launch::async, [begin, end, &func]() {
            func(begin, end);
        }));
    }
    for (auto& worker : workers) {
        worker.get();
    }
}

template <class Func>
inline auto parallel_collect_chunks(std::size_t task_count, std::size_t min_tasks_per_thread, Func&& func)
    -> std::vector<decltype(func(std::size_t{}, std::size_t{}))> {
    using Result = decltype(func(std::size_t{}, std::size_t{}));
    const std::size_t thread_count = automatic_thread_count(task_count, min_tasks_per_thread);
    if (thread_count <= 1 || task_count <= 1) {
        std::vector<Result> results;
        if (task_count != 0) {
            results.push_back(func(0, task_count));
        }
        return results;
    }

    const std::size_t chunk_size = (task_count + thread_count - 1) / thread_count;
    std::vector<std::future<Result>> workers;
    workers.reserve(thread_count);
    for (std::size_t begin = 0; begin < task_count; begin += chunk_size) {
        const std::size_t end = std::min(task_count, begin + chunk_size);
        workers.push_back(std::async(std::launch::async, [begin, end, &func]() {
            return func(begin, end);
        }));
    }

    std::vector<Result> results;
    results.reserve(workers.size());
    for (auto& worker : workers) {
        results.push_back(worker.get());
    }
    return results;
}

}  // namespace xpongecpp

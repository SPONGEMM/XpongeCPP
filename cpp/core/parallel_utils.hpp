#pragma once

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <future>
#include <optional>
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

inline std::size_t compute_chunk_size(std::size_t task_count, std::size_t chunk_hint, std::size_t thread_count) {
    const std::size_t hint = std::max<std::size_t>(1, chunk_hint);
    const std::size_t target_chunk_count = std::max<std::size_t>(thread_count * 4, 1);
    const std::size_t balanced = std::max<std::size_t>(1, (task_count + target_chunk_count - 1) / target_chunk_count);
    return std::max(hint, balanced);
}

template <class Func>
inline void parallel_for_chunks(std::size_t task_count, std::size_t chunk_hint, Func&& func) {
    const std::size_t thread_count = automatic_thread_count(task_count, chunk_hint);
    if (thread_count <= 1 || task_count <= 1) {
        if (task_count != 0) {
            func(0, task_count);
        }
        return;
    }

    const std::size_t chunk_size = compute_chunk_size(task_count, chunk_hint, thread_count);
    const std::size_t chunk_count = (task_count + chunk_size - 1) / chunk_size;
    std::atomic<std::size_t> next_chunk{0};
    std::vector<std::future<void>> workers;
    workers.reserve(thread_count);
    for (std::size_t worker_id = 0; worker_id < thread_count; ++worker_id) {
        workers.push_back(std::async(std::launch::async, [&func, &next_chunk, chunk_count, chunk_size, task_count]() {
            while (true) {
                const std::size_t chunk_index = next_chunk.fetch_add(1, std::memory_order_relaxed);
                if (chunk_index >= chunk_count) {
                    return;
                }
                const std::size_t begin = chunk_index * chunk_size;
                const std::size_t end = std::min(task_count, begin + chunk_size);
                func(begin, end);
            }
        }));
    }
    for (auto& worker : workers) {
        worker.get();
    }
}

template <class Func>
inline auto parallel_collect_chunks(std::size_t task_count, std::size_t chunk_hint, Func&& func)
    -> std::vector<decltype(func(std::size_t{}, std::size_t{}))> {
    using Result = decltype(func(std::size_t{}, std::size_t{}));
    const std::size_t thread_count = automatic_thread_count(task_count, chunk_hint);
    if (thread_count <= 1 || task_count <= 1) {
        std::vector<Result> results;
        if (task_count != 0) {
            results.push_back(func(0, task_count));
        }
        return results;
    }

    const std::size_t chunk_size = compute_chunk_size(task_count, chunk_hint, thread_count);
    const std::size_t chunk_count = (task_count + chunk_size - 1) / chunk_size;
    std::atomic<std::size_t> next_chunk{0};
    std::vector<std::optional<Result>> chunk_results(chunk_count);
    std::vector<std::future<void>> workers;
    workers.reserve(thread_count);
    for (std::size_t worker_id = 0; worker_id < thread_count; ++worker_id) {
        workers.push_back(std::async(std::launch::async,
                                     [&func, &next_chunk, &chunk_results, chunk_count, chunk_size, task_count]() {
            while (true) {
                const std::size_t chunk_index = next_chunk.fetch_add(1, std::memory_order_relaxed);
                if (chunk_index >= chunk_count) {
                    return;
                }
                const std::size_t begin = chunk_index * chunk_size;
                const std::size_t end = std::min(task_count, begin + chunk_size);
                chunk_results[chunk_index] = func(begin, end);
            }
        }));
    }
    for (auto& worker : workers) {
        worker.get();
    }

    std::vector<Result> results;
    results.reserve(chunk_count);
    for (auto& chunk_result : chunk_results) {
        if (chunk_result.has_value()) {
            results.push_back(std::move(*chunk_result));
        }
    }
    return results;
}

}  // namespace xpongecpp

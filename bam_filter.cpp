#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/thread_pool.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <mutex>
#include <thread>
#include <atomic>
#include <condition_variable>

const int CHUNK_SIZE = 10000;

struct ReadPair {
    bam1_t *read1;
    bam1_t *read2;
};

struct ThreadData {
    std::vector<ReadPair> readPairs;
    std::mutex mutex;
    std::condition_variable cv;
    bool ready = false;
    bool done = false;
    std::atomic<bool> error{false};
};

void processReads(ThreadData &data, samFile *out, bam_hdr_t *header, int mapq_threshold) {
    std::vector<ReadPair> localPairs;

    while (true) {
        {
            std::unique_lock<std::mutex> lock(data.mutex);
            data.cv.wait(lock, [&]{ return data.ready || data.done; });

            if (data.done && data.readPairs.empty()) {
                return;
            }

            localPairs.swap(data.readPairs);
            data.ready = false;
        }

        for (auto &pair : localPairs) {
            if (pair.read1->core.qual >= mapq_threshold && pair.read2->core.qual >= mapq_threshold) {
                if (sam_write1(out, header, pair.read1) < 0 || sam_write1(out, header, pair.read2) < 0) {
                    fprintf(stderr, "Error writing to output BAM file\n");
                    data.error = true;
                    return;
                }
            }
            bam_destroy1(pair.read1);
            bam_destroy1(pair.read2);
        }

        localPairs.clear();
    }
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s input.bam output.bam MAPQ_threshold [num_threads]\n", argv[0]);
        return 1;
    }
    
    // Parse command line arguments
    char *input_bam = argv[1];
    char *output_bam = argv[2];
    int mapq_threshold = atoi(argv[3]);
    int num_threads = 1;
    if (argc > 4) {
        num_threads = atoi(argv[4]);
    }

    // Open input and output BAM files
    samFile *in = sam_open(input_bam, "rb");
    if (!in) {
        fprintf(stderr, "Failed to open input BAM file %s\n", input_bam);
        return 1;
    }

    // Read header from input BAM file
    bam_hdr_t *header = sam_hdr_read(in);
    if (!header) {
        fprintf(stderr, "Failed to read header from %s\n", input_bam);
        sam_close(in);
        return 1;
    }

    // Open output BAM file
    samFile *out = sam_open(output_bam, "wb");
    if (!out) {
        fprintf(stderr, "Failed to open output BAM file %s\n", output_bam);
        bam_hdr_destroy(header);
        sam_close(in);
        return 1;
    }

    // Set up multithreading
    htsThreadPool thread_pool = {NULL, 0};
    if (num_threads > 1) {
        thread_pool.pool = hts_tpool_init(num_threads);
        if (!thread_pool.pool) {
            fprintf(stderr, "Failed to create thread pool\n");
            bam_hdr_destroy(header);
            sam_close(in);
            sam_close(out);
            return 1;
        }
        // Set thread pool for input and output files
        hts_set_thread_pool(in, &thread_pool);
        hts_set_thread_pool(out, &thread_pool);
    }

    // Write header to output BAM file
    if (sam_hdr_write(out, header) < 0) {
        fprintf(stderr, "Failed to write header to %s\n", output_bam);
        bam_hdr_destroy(header);
        sam_close(in);
        sam_close(out);
        return 1;
    }

    std::vector<ThreadData> threadData(num_threads);
    std::vector<std::thread> threads;

    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back(processReads, std::ref(threadData[i]), out, header, mapq_threshold);
    }

    int current_thread = 0;
    int ret;
    bam1_t *b = bam_init1();

    // Read and process reads in pairs
    while ((ret = sam_read1(in, header, b)) >= 0) {
        ThreadData &data = threadData[current_thread];
        
        {
            std::lock_guard<std::mutex> lock(data.mutex);
            data.readPairs.push_back({bam_dup1(b), bam_init1()});
        }

        // Read second read of pair
        if ((ret = sam_read1(in, header, data.readPairs.back().read2)) < 0) {
            fprintf(stderr, "Unexpected end of file or error reading second read of pair\n");
            break;
        }

        if (data.readPairs.size() >= CHUNK_SIZE) {
            {
                std::lock_guard<std::mutex> lock(data.mutex);
                data.ready = true;
            }
            data.cv.notify_one();
            current_thread = (current_thread + 1) % num_threads;
        }

        // Check for errors in worker threads
        bool error_occurred = false;
        for (const auto &data : threadData) {
            if (data.error) {
                error_occurred = true;
                break;
            }
        }
        if (error_occurred) {
            fprintf(stderr, "Error occurred in worker thread. Stopping processing.\n");
            break;
        }
    }

    // Signal all threads to finish processing
    for (auto &data : threadData) {
        {
            std::lock_guard<std::mutex> lock(data.mutex);
            data.ready = true;
            data.done = true;
        }
        data.cv.notify_one();
    }

    // Wait for all threads to finish
    for (auto &thread : threads) {
        thread.join();
    }

    // Check for any errors in worker threads
    bool error_occurred = false;
    for (const auto &data : threadData) {
        if (data.error) {
            error_occurred = true;
            break;
        }
    }

    // Clean up
    bam_destroy1(b);
    bam_hdr_destroy(header);
    sam_close(in);
    sam_close(out);

    // Clean up thread pool
    if (thread_pool.pool) {
        hts_tpool_destroy(thread_pool.pool);
    }

    return error_occurred ? 1 : 0;
}
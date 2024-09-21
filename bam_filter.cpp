#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/thread_pool.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>  // for uint64_t

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

    // Open input file (SAM format from stdin if "-")
    samFile *in = sam_open(strcmp(input_bam, "-") == 0 ? "-" : input_bam, "r"); // "r" for SAM input
    if (!in) {
        fprintf(stderr, "Failed to open input SAM file %s\n", input_bam);
        return 1;
    }

    // Read header from input SAM file
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

    // Use uint64_t for large counter values
    uint64_t single_end_count = 0;
    uint64_t unmapped_mate_count = 0;
    uint64_t low_mapq_count = 0;
    uint64_t written_count = 0;

    // Read and write reads in pairs
    bam1_t *b1 = bam_init1();
    bam1_t *b2 = bam_init1();
    int ret1, ret2;

    // Read first read of pair
    while ((ret1 = sam_read1(in, header, b1)) >= 0) {
        // Read second read of pair
        ret2 = sam_read1(in, header, b2);

        // Check for unexpected end of file or error reading second read
        if (ret2 < 0) {
            fprintf(stderr, "Unexpected end of file or error reading second read of pair\n");
            break;
        }

        // Check if both reads are paired (flag 0x1 set) and mate is not unmapped (flag 0x8 not set)
        bool b1_paired = (b1->core.flag & 0x1); // Check if first read is paired
        bool b2_paired = (b2->core.flag & 0x1); // Check if second read is paired
        bool b1_mate_unmapped = (b1->core.flag & 0x8); // Check if mate of first read is unmapped
        bool b2_mate_unmapped = (b2->core.flag & 0x8); // Check if mate of second read is unmapped

        // Filter out single-end reads or pairs where one mate is unmapped
        if (b1_paired && b2_paired) {
            if (!b1_mate_unmapped && !b2_mate_unmapped) {
                // Check if both reads have MAPQ >= threshold
                int mapq1 = b1->core.qual;
                int mapq2 = b2->core.qual;

                if (mapq1 >= mapq_threshold && mapq2 >= mapq_threshold) {
                    // Write reads to output if both have MAPQ >= threshold
                    if (sam_write1(out, header, b1) < 0 || sam_write1(out, header, b2) < 0) {
                        fprintf(stderr, "Failed to write reads to output\n");
                        break;
                    }
                    written_count += 2;  // Increment count for both reads in the pair
                } else {
                    low_mapq_count += 2;  // Both reads are discarded for low MAPQ
                }
            } else {
                unmapped_mate_count += 2;  // Discarded because one or both mates are unmapped
            }
        } else {
            single_end_count += 2;  // Discarded as single-end reads (not paired)
        }
    }

    // Clean up
    bam_destroy1(b1);
    bam_destroy1(b2);
    bam_hdr_destroy(header);
    sam_close(in);
    sam_close(out);

    // Clean up thread pool
    if (thread_pool.pool) {
        hts_tpool_destroy(thread_pool.pool);
    }

    // Report results
    printf("Report:\n");
    printf("Single-end reads discarded: %" PRIu64 "\n", single_end_count);
    printf("Reads discarded due to unmapped mate: %" PRIu64 "\n", unmapped_mate_count);
    printf("Reads discarded due to low MAPQ: %" PRIu64 "\n", low_mapq_count);
    printf("Reads written to output BAM: %" PRIu64 "\n", written_count);

    return 0;
}

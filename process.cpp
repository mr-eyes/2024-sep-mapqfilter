#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/thread_pool.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

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

        // Check if both reads have MAPQ >= threshold
        int mapq1 = b1->core.qual;
        int mapq2 = b2->core.qual;

        // Write reads to output if both have MAPQ >= threshold
        if (mapq1 >= mapq_threshold && mapq2 >= mapq_threshold) {
            if (sam_write1(out, header, b1) < 0 || sam_write1(out, header, b2) < 0) {
                fprintf(stderr, "Failed to write reads to output\n");
                break;
            }
        }
        // Else discard both reads
    }

    // Check for errors reading first read
    bam_destroy1(b1);
    bam_destroy1(b2);
    bam_hdr_destroy(header);
    sam_close(in);
    sam_close(out);

    // Clean up thread pool
    if (thread_pool.pool) {
        hts_tpool_destroy(thread_pool.pool);
    }

    return 0;
}
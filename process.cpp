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
    char *input_bam = argv[1];
    char *output_bam = argv[2];
    int mapq_threshold = atoi(argv[3]);
    int num_threads = 1;
    if (argc > 4) {
        num_threads = atoi(argv[4]);
    }

    samFile *in = sam_open(input_bam, "rb");
    if (!in) {
        fprintf(stderr, "Failed to open input BAM file %s\n", input_bam);
        return 1;
    }
    bam_hdr_t *header = sam_hdr_read(in);
    if (!header) {
        fprintf(stderr, "Failed to read header from %s\n", input_bam);
        sam_close(in);
        return 1;
    }

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
        hts_set_thread_pool(in, &thread_pool);
        hts_set_thread_pool(out, &thread_pool);
    }

    if (sam_hdr_write(out, header) < 0) {
        fprintf(stderr, "Failed to write header to %s\n", output_bam);
        bam_hdr_destroy(header);
        sam_close(in);
        sam_close(out);
        return 1;
    }

    bam1_t *b1 = bam_init1();
    bam1_t *b2 = bam_init1();
    int ret1, ret2;

    while ((ret1 = sam_read1(in, header, b1)) >= 0) {
        ret2 = sam_read1(in, header, b2);
        if (ret2 < 0) {
            fprintf(stderr, "Unexpected end of file or error reading second read of pair\n");
            break;
        }

        // Check that the two reads are paired
        const char *qname1 = bam_get_qname(b1);
        const char *qname2 = bam_get_qname(b2);
        if (strcmp(qname1, qname2) != 0) {
            fprintf(stderr, "Read names do not match: %s vs %s\n", qname1, qname2);
            continue; // Skip to the next pair
        }

        int mapq1 = b1->core.qual;
        int mapq2 = b2->core.qual;

        if (mapq1 >= mapq_threshold && mapq2 >= mapq_threshold) {
            if (sam_write1(out, header, b1) < 0 || sam_write1(out, header, b2) < 0) {
                fprintf(stderr, "Failed to write reads to output\n");
                break;
            }
        }
        // Else discard both reads
    }

    bam_destroy1(b1);
    bam_destroy1(b2);
    bam_hdr_destroy(header);
    sam_close(in);
    sam_close(out);

    if (thread_pool.pool) {
        hts_tpool_destroy(thread_pool.pool);
    }

    return 0;
}

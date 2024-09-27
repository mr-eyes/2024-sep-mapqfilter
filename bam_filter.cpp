// filter_bam.cpp

#include <iostream>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <string>
#include <cstdlib>
#include <cstdio>

int main(int argc, char *argv[])
{
    // Check and parse command-line arguments
    if (argc < 4)
    {
        std::cerr << "Usage: " << argv[0] << " input.bam output.bam mapq_cutoff\n";
        return 1;
    }
    const char *input_bam = argv[1];
    const char *output_bam = argv[2];
    int mapq_cutoff = atoi(argv[3]);

    int n_threads = 128; // Number of threads to use

    // Open input BAM file
    htsFile *in = hts_open(input_bam, "r");
    if (in == NULL)
    {
        std::cerr << "Error opening input BAM file\n";
        return 1;
    }
    // Set threads for reading
    hts_set_threads(in, n_threads);

    // Read header from input BAM file
    bam_hdr_t *header = sam_hdr_read(in);
    if (header == NULL)
    {
        std::cerr << "Error reading header from input BAM file\n";
        return 1;
    }

    // Open output BAM file
    htsFile *out = hts_open(output_bam, "wb");
    if (out == NULL)
    {
        std::cerr << "Error opening output BAM file\n";
        return 1;
    }
    // Set threads for writing
    hts_set_threads(out, n_threads);

    // Write header to output BAM file
    if (sam_hdr_write(out, header) < 0)
    {
        std::cerr << "Error writing header to output BAM file\n";
        return 1;
    }

    // Variables for filtration statistics
    int64_t total_pairs = 0;
    int64_t passed_pairs = 0;
    int64_t failed_unmapped = 0;
    int64_t failed_not_paired = 0;
    int64_t failed_mapq = 0;
    int64_t failed_pos_pnext = 0;

    bam1_t *b1 = bam_init1();
    bam1_t *b2 = bam_init1();

    // Process the BAM file two reads at a time (assuming mates are adjacent)
    while (true)
    {
        // Read first read
        if (sam_read1(in, header, b1) < 0)
            break; // End of file
        // Read second read
        if (sam_read1(in, header, b2) < 0)
            break; // End of file

        total_pairs++;

        // Retrieve flags
        uint32_t flag1 = b1->core.flag;
        uint32_t flag2 = b2->core.flag;

        // Check if reads are unmapped
        bool unmapped1 = (flag1 & BAM_FUNMAP) != 0;
        bool unmapped2 = (flag2 & BAM_FUNMAP) != 0;

        // Check if reads are paired
        bool paired1 = (flag1 & BAM_FPAIRED) != 0;
        bool paired2 = (flag2 & BAM_FPAIRED) != 0;

        // Get MAPQ scores
        uint8_t mapq1 = b1->core.qual;
        uint8_t mapq2 = b2->core.qual;

        // Get POS and PNEXT
        int32_t pos1 = b1->core.pos;    // 0-based leftmost coordinate
        int32_t pnext1 = b1->core.mpos; // Position of the mate
        int32_t pos2 = b2->core.pos;
        int32_t pnext2 = b2->core.mpos;

        // Apply filters
        if (unmapped1 || unmapped2)
        {
            failed_unmapped++;
            continue;
        }

        if (!paired1 || !paired2)
        {
            failed_not_paired++;
            continue;
        }

        if (mapq1 <= mapq_cutoff || mapq2 <= mapq_cutoff)
        {
            failed_mapq++;
            continue;
        }

        if (pos1 != pnext2 || pos2 != pnext1)
        {
            failed_pos_pnext++;
            continue;
        }

        // Passed all filters, write reads to output BAM file
        if (sam_write1(out, header, b1) < 0)
        {
            std::cerr << "Error writing to output BAM file\n";
            return 1;
        }
        if (sam_write1(out, header, b2) < 0)
        {
            std::cerr << "Error writing to output BAM file\n";
            return 1;
        }

        passed_pairs++;
    }

    // Clean up
    bam_destroy1(b1);
    bam_destroy1(b2);
    bam_hdr_destroy(header);
    hts_close(in);
    hts_close(out);

    // Output filtration statistics
    std::cout << "Filtration Statistics:\n";
    std::cout << "Total pairs processed: " << total_pairs << "\n";
    std::cout << "Passed pairs: " << passed_pairs << "\n";
    std::cout << "Failed due to unmapped: " << failed_unmapped << "\n";
    std::cout << "Failed due to not paired: " << failed_not_paired << "\n";
    std::cout << "Failed due to low MAPQ: " << failed_mapq << "\n";
    std::cout << "Failed due to POS != PNEXT: " << failed_pos_pnext << "\n";

    return 0;
}

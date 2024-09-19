# BAM Filter: Setup, Compilation, and Usage Guide

## Environment Setup

```
mamba create -n bamfilter bioconda::samtools bioconda::htslib zlib gcc 
```

## Compilation

```
g++ -O2 -o filter_bam bam_filter.cpp -I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib -lhts -lz -lpthread -Wl,-rpath,$CONDA_PREFIX/lib
```

## Usage

Run the program with the following command:
```
./filter_bam input.bam output.bam MAPQ_threshold [num_threads]
```

- `input.bam`: Path to the input BAM file
- `output.bam`: Path for the output BAM file
- `MAPQ_threshold`: Minimum MAPQ value to keep a read pair
- `num_threads` (optional): Number of threads to use (default is 1)

Example:
```
./filter_bam test_reads.bam filtered_reads.bam 20 4
```

Example: This command filters `all_reads.bam`, keeping only read pairs where both reads have a MAPQ score of at least 20, using 4 threads, and saves the result to `filtered_reads.bam`.

## Troubleshooting

If you encounter library loading errors, ensure that the Conda environment is activated and try setting the `LD_LIBRARY_PATH` (Linux) or `DYLD_LIBRARY_PATH` (macOS):

```
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH  # Linux
export DYLD_LIBRARY_PATH=$CONDA_PREFIX/lib:$DYLD_LIBRARY_PATH  # macOS
```

Then run the program again.
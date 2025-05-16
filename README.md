# levx
Computing ultralow-res pairwise Levehnstein distances across large chromosomes

Resolution is as follows:
* 10bp when the distance between two loci is <100kb
* 100bp when the distance between two loci is <1Mb
* 1kb when the distance between two loci is >1Mb

This can be changed in the source before building.

Build: `cargo build --release`
Usage: `./chromosome_distance_calculator <fasta_file> <output_ipc_file>`

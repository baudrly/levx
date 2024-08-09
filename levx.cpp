#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <omp.h>
#include <algorithm>
#include <Eigen/Sparse>

// Using Eigen::SparseMatrix
// using SparseMatrix = Eigen::SparseMatrix<int>;

enum Resolution { RES_10BP, RES_100BP, RES_1KB };

Resolution decide_resolution(size_t distance) {
    if (distance <= 100000) return RES_10BP;
    if (distance <= 1000000) return RES_100BP;
    return RES_1KB;
}

// using only two rows
int levenshtein_distance_optimized(const std::string &s1, const std::string &s2) {
    size_t len1 = s1.size(), len2 = s2.size();
    if (len1 > len2) {
        return levenshtein_distance_optimized(s2, s1);
    }

    std::vector<int> current(len1 + 1), previous(len1 + 1);
    for (size_t i = 0; i <= len1; i++) {
        previous[i] = i;
    }

    for (size_t j = 1; j <= len2; j++) {
        current[0] = j;
        for (size_t i = 1; i <= len1; i++) {
            int cost = (s1[i - 1] == s2[j - 1]) ? 0 : 1;
            current[i] = std::min({previous[i] + 1, current[i - 1] + 1, previous[i - 1] + cost});
        }
        std::swap(current, previous);
    }

    return previous[len1];
}

std::string load_genome_from_fasta(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open FASTA file: " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line, genome_sequence;
    std::getline(file, line); // Skip header line
    while (std::getline(file, line)) {
        genome_sequence += line;
    }

    return genome_sequence;
}

//void save_chunk_to_csv(const SparseMatrix& matrix, const std::string& filename, size_t start, size_t end) {
//    std::ofstream out_file(filename);
//    if (!out_file.is_open()) {
//        throw std::runtime_error("Unable to open file for writing");
//    }
//
//    // iterate over nonzero elements of the sparse matrix within chunk
//    for (int k = start; k < end; ++k) {
//        for (SparseMatrix::InnerIterator it(matrix, k); it; ++it) {
//            if (it.row() >= start && it.row() < end) { // Only save elements within the chunk
//                out_file << it.row() << "," << it.col() << "," << it.value() << "\n";
//            }
//        }
//    }
//}
void save_result(size_t i, size_t j, int distance, const std::string& filename) {
    std::ofstream out_file(filename, std::ios::app); // Append to the file
    if (!out_file.is_open()) {
        throw std::runtime_error("Unable to open file for writing");
    }
    out_file << i << "," << j << "," << distance << "\n";
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input.fasta> <output.csv>\n";
        return 1;
    }

    int max_threads = omp_get_max_threads();
    omp_set_num_threads(max_threads);

    std::string chromosome_data = load_genome_from_fasta(argv[1]);
    size_t data_size = chromosome_data.size();

    size_t chunk_size = 10000; // genome chunks
    const size_t min_chunk_size = 1000;

    // Process the genome in chunks
    for (size_t start = 0; start < data_size; start += chunk_size) {
        size_t end = std::min(start + chunk_size, data_size);

        #pragma omp parallel for schedule(dynamic)
        for (size_t i = start; i < end; i++) {
            for (size_t j = i; j < data_size; j++) { 
                Resolution resolution = decide_resolution(j - i);
                int segment_size = (resolution == RES_10BP) ? 10 : (resolution == RES_100BP) ? 100 : 1000;

                // Avoid out-of-bounds access
                if (i + segment_size > data_size || j + segment_size > data_size) {
                    continue;
                }

                std::string segment1 = chromosome_data.substr(i, segment_size);
                std::string segment2 = chromosome_data.substr(j, segment_size);
                int distance = levenshtein_distance_optimized(segment1, segment2);

                // result directly to output
                save_result(i, j, distance, argv[2]);
            }
        }

        // Dynamic Chunk Size Adjustment (Optional - uncomment if needed)
        // // Monitor memory usage within the chunk
        // // ... (Implementation depends on your system/OS)

        // // Adjust chunk size based on memory usage
        // if (/* memory usage is high */) {
        //     chunk_size = std::max(chunk_size / 2, min_chunk_size);
        // } else if (/* memory usage is low */) {
        //     chunk_size = std::min(chunk_size * 2, data_size - start);
        // }
    }

    std::cout << "Processing complete. Results saved to: " << argv[2] << "\n";
    return 0;
}

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <Eigen/Sparse>
#include <omp.h>

enum Resolution { RES_10BP, RES_100BP, RES_1KB };

Resolution decide_resolution(size_t distance) {
    if (distance <= 100000) return RES_10BP;
    if (distance <= 1000000) return RES_100BP;
    return RES_1KB;
}

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
            current[i] = std::min({ previous[i] + 1, current[i - 1] + 1, previous[i - 1] + cost });
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
    std::getline(file, line);
    if (line.empty() || line[0] != '>') {
        std::cerr << "Invalid FASTA file format." << std::endl;
        exit(EXIT_FAILURE);
    }

    while (std::getline(file, line)) {
        genome_sequence += line;
    }

    return genome_sequence;
}

void save_to_csv(const Eigen::SparseMatrix<int>& matrix, const std::string& filename) {
    std::ofstream out_file(filename);
    if (!out_file.is_open()) {
        throw std::runtime_error("Unable to open file for writing");
    }

    for (int k = 0; k < matrix.outerSize(); ++k) {
        for (Eigen::SparseMatrix<int>::InnerIterator it(matrix, k); it; ++it) {
            out_file << it.row() << "," << it.col() << "," << it.value() << "\n";
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input.fasta> <output.csv>\n";
        return 1;
    }

    omp_set_num_threads(omp_get_max_threads());
    std::string chromosome_data = load_genome_from_fasta(argv[1]);

    Eigen::SparseMatrix<int> matrix(chromosome_data.size(), chromosome_data.size());

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < chromosome_data.size(); i++) {
        for (size_t j = i; j < chromosome_data.size(); j++) {
            Resolution resolution = decide_resolution(j - i);
            int segment_size = (resolution == RES_10BP) ? 10 : (resolution == RES_100BP) ? 100 : 1000;

            if (i + segment_size > chromosome_data.size() || j + segment_size > chromosome_data.size()) {
                continue;
            }

            std::string segment1 = chromosome_data.substr(i, segment_size);
            std::string segment2 = chromosome_data.substr(j, segment_size);
            int distance = levenshtein_distance_optimized(segment1, segment2);

            matrix.coeffRef(i, j) = distance;
            if (i != j) {
                matrix.coeffRef(j, i) = distance;
            }
        }
    }

    save_to_csv(matrix, argv[2]);
    std::cout << "Processing complete. Data saved to: " << argv[2] << "\n";
    return 0;
}

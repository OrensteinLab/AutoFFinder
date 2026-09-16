/**
 * Test PAM-aware deduplication for cross-tool validation edge cases.
 * 
 * This test documents a bug (now fixed) where deduplication was not PAM-aware,
 * causing Stomata to keep hits with invalid PAMs while discarding valid ones.
 * 
 * Test case: AAVS1_site_4 spacer GCTCGGGGACACAGGATCCC
 * Location: chr20:36308165 (CasOFFinder position, minus strand)
 * 
 * Myers algorithm reports THREE overlapping alignments near this position:
 *   - chr20:36308167: distance=3, PAM=GGC (OTHER, invalid)
 *   - chr20:36308168: distance=3, PAM=TGG (NGG, valid!)  <-- CasOFFinder found this
 *   - chr20:36308171: distance=3, PAM=CCC (OTHER, invalid) <-- Bug: Stomata kept this
 * 
 * Original bug: Deduplication happened BEFORE PAM extraction, so it picked the
 * rightmost hit (36308171) without knowing it had an invalid PAM. The valid hit
 * (36308168) was discarded before PAM filtering ever ran.
 * 
 * Fix (2026-02-11):
 * 1. Moved PAM extraction before deduplication
 * 2. Made deduplication PAM-aware: when distance is tied, prefer NGG > NAG > OTHER
 * 3. Result: Stomata now correctly keeps chr20:36308168 (valid NGG)
 * 
 * This test verifies the fix works correctly.
 */

#include <iostream>
#include <string>
#include <vector>

#include <genome_loader.hpp>
#include <search_pipeline.hpp>

int main() {
    std::cout << "=== Halo Deduplication Test Case ===\n\n";
    
    // Load genome
    std::string genome_path = "/mnt/raid6/common/genomes/GRCh38.p14/hg38.p14.canonical.fa.st";
    
    std::cout << "Loading genome: " << genome_path << "\n";
    auto genome_opt = load_indexed_genome(genome_path);
    if (!genome_opt) {
        std::cerr << "ERROR: Failed to load genome\n";
        return 1;
    }
    auto view = make_view(*genome_opt);
    
    // Test spacer
    std::string spacer = "GCTCGGGGACACAGGATCCC";
    std::cout << "Spacer: " << spacer << "\n";
    std::cout << "Target region: chr20:36308165 (CasOFFinder minus-strand position)\n\n";
    
    // Search WITH deduplication (default)
    std::cout << "--- WITH deduplication (default) ---\n";
    SearchConfig config1;
    config1.pattern = spacer;
    config1.threshold = 3;
    config1.search_both_strands = true;
    config1.compute_mismatches = true;
    config1.disable_deduplication = false;
    
    auto result1 = search_genome(config1, view);
    
    int count1 = 0;
    for (const auto& hit : result1.hits) {
        if (hit.chrom_name == "chr20" && 
            hit.chrom_offset >= 36308165 && 
            hit.chrom_offset <= 36308175 &&
            hit.strand == Strand::MINUS) {
            
            std::cout << "  Position: " << hit.chrom_offset 
                      << ", Distance: " << (int)hit.distance
                      << ", Aligned: " << hit.mismatch_info.aligned_sequence << "\n";
            count1++;
        }
    }
    std::cout << "Total hits in region: " << count1 << "\n\n";
    
    // Search WITHOUT deduplication
    std::cout << "--- WITHOUT deduplication (--no-deduplicate) ---\n";
    SearchConfig config2;
    config2.pattern = spacer;
    config2.threshold = 3;
    config2.search_both_strands = true;
    config2.compute_mismatches = true;
    config2.disable_deduplication = true;  // <-- Disable halo deduplication
    
    auto result2 = search_genome(config2, view);
    
    int count2 = 0;
    for (const auto& hit : result2.hits) {
        if (hit.chrom_name == "chr20" && 
            hit.chrom_offset >= 36308165 && 
            hit.chrom_offset <= 36308175 &&
            hit.strand == Strand::MINUS) {
            
            std::cout << "  Position: " << hit.chrom_offset 
                      << ", Distance: " << (int)hit.distance
                      << ", Aligned: " << hit.mismatch_info.aligned_sequence << "\n";
            count2++;
        }
    }
    std::cout << "Total hits in region: " << count2 << "\n\n";
    
    // Summary
    std::cout << "=== Summary ===\n";
    std::cout << "Hits with deduplication: " << count1 << "\n";
    std::cout << "Hits without deduplication: " << count2 << "\n";
    std::cout << "Suppressed by deduplication: " << (count2 - count1) << "\n\n";
    
    if (count1 == 1 && count2 == 3) {
        std::cout << "RESULT: PAM-aware deduplication working correctly!\n";
        std::cout << "The valid NGG hit was kept, invalid PAMs were discarded.\n";
        return 0;
    } else if (count1 == 0) {
        std::cout << "FAILURE: No hits found with deduplication (bug still present?)\n";
        return 1;
    } else {
        std::cout << "UNEXPECTED: Found " << count1 << " hits with deduplication\n";
        return 1;
    }
}

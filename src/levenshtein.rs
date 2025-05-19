/// Computes Levenshtein distance between two byte slices.
/// Uses O(min(m,n)) space and O(m*n) time.
pub fn levenshtein_distance(s1: &[u8], s2: &[u8]) -> u16 {
    if s1.is_empty() { return s2.len() as u16; }
    if s2.is_empty() { return s1.len() as u16; }

    // Ensure s1 is the shorter sequence for space optimization
    let (s1_effective, s2_effective) = if s1.len() <= s2.len() { (s1, s2) } else { (s2, s1) };

    let len1 = s1_effective.len();
    let len2 = s2_effective.len();

    let mut previous_row: Vec<u16> = (0..=len1 as u16).collect();
    let mut current_row: Vec<u16> = vec![0; len1 + 1];

    for j_idx in 1..=len2 {
        current_row[0] = j_idx as u16;
        for i_idx in 1..=len1 {
            let cost = if s1_effective[i_idx - 1] == s2_effective[j_idx - 1] { 0 } else { 1 };
            current_row[i_idx] = std::cmp::min(
                previous_row[i_idx] + 1, // Deletion from s1
                std::cmp::min(
                    current_row[i_idx - 1] + 1, // Insertion into s1
                    previous_row[i_idx - 1] + cost, // Substitution or match
                ),
            );
        }
        std::mem::swap(&mut previous_row, &mut current_row);
    }

    previous_row[len1]
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_levenshtein_empty() {
        assert_eq!(levenshtein_distance(b"", b""), 0);
        assert_eq!(levenshtein_distance(b"abc", b""), 3);
        assert_eq!(levenshtein_distance(b"", b"abc"), 3);
    }

    #[test]
    fn test_levenshtein_basic() {
        assert_eq!(levenshtein_distance(b"kitten", b"sitting"), 3);
        assert_eq!(levenshtein_distance(b"ACGT", b"AGGT"), 1);
    }

    #[test]
    fn test_levenshtein_identical() {
        assert_eq!(levenshtein_distance(b"abc", b"abc"), 0);
    }
     #[test]
    fn test_levenshtein_order_invariant() {
        assert_eq!(levenshtein_distance(b"longstring", b"short"), levenshtein_distance(b"short", b"longstring"));
    }
}
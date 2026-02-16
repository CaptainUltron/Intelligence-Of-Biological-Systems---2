from Bio.Seq import Seq

# DNA sequence
text = ("GTTCCGAAAGGCTAGCGCTAGGCGCCAAGCGGCCGGTTTCCTTGGCGACGGAGAGCGCGGGAATTTTAGA"
        "TAGATTGTAATTGCGGCTGCGCGGCCGCTGCCCGTGCAGCCAGAGGATCCAGCACCTCTCTTGGGGCTTC"
        "TCCGTCCTCGGCGCTTGGAAGTACGGATCTTTTTTCTCGGAGAAAAGTTCACTGGAACTGGAAGAAATGG")

# Convert DNA to protein
dna_seq = Seq(text)
protein_seq = str(dna_seq.translate(to_stop=False)).replace("*", "")

print("Protein Sequence:")
print(protein_seq)
print()

# Amino acid molecular weights
AA_MW = {
    'G':57, 'A':71, 'S':87, 'P':97, 'V':99,
    'T':101, 'C':103, 'I':113, 'L':113, 'N':114,
    'D':115, 'K':128, 'Q':128, 'E':129, 'M':131,
    'H':137, 'F':147, 'R':156, 'Y':163, 'W':186
}

# Length of peptide
k = 4

def branch_and_bound(protein, k):
    max_weight = 0
    best = ""
    
    # Maximum possible weight of any amino acid (for upper bound)
    MAX_AA = max(AA_MW.values())

    for i in range(len(protein) - k + 1):
        weight = 0

        for j in range(k):
            # Add current amino acid weight
            weight += AA_MW[protein[i + j]]

            # Remaining positions
            remaining = k - j - 1

            # Upper bound if all remaining were the heaviest amino acid
            upper_bound = weight + remaining * MAX_AA

            # Prune: cannot beat current best
            if upper_bound <= max_weight:
                break

        # If full length processed and better weight found
        if j == k - 1 and weight > max_weight:
            max_weight = weight
            best = protein[i:i+k]

    return best, max_weight


# Run Branch and Bound
result = branch_and_bound(protein_seq, k)
print("Branch & Bound Result:")
print("Best peptide:", result[0])
print("Weight:", result[1])
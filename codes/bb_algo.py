from Bio.Seq import Seq

# Your DNA input
text = ("GTTCCGAAAGGCTAGCGCTAGGCGCCAAGCGGCCGGTTTCCTTGGCGACGGAGAGCGCGGGAATTTTAGA"
    "TAGATTGTAATTGCGGCTGCGCGGCCGCTGCCCGTGCAGCCAGAGGATCCAGCACCTCTCTTGGGGCTTC"
    "TCCGTCCTCGGCGCTTGGAAGTACGGATCTTTTTTCTCGGAGAAAAGTTCACTGGAACTGGAAGAAATGG")

# Convert DNA → Protein using Bio
dna_seq = Seq(text)
protein_seq = str(dna_seq.translate(to_stop=False)).replace("*", "")

print("Protein Sequence:")
print(protein_seq)

AA_MW = {
'A':89,'R':174,'N':132,'D':133,'C':121,'Q':146,'E':147,'G':75,
'H':155,'I':131,'L':131,'K':146,'M':149,'F':165,'P':115,'S':105,
'T':119,'W':204,'Y':181,'V':117
}

k = 4

def branch_and_bound(protein, k):
    max_weight = 0
    best = ""

    for i in range(len(protein)-k+1):
        weight = 0

        for j in range(k):
            weight += AA_MW[protein[i+j]]

            # prune if already smaller
            if weight < max_weight and j == k-1:
                break

        if weight > max_weight:
            max_weight = weight
            best = protein[i:i+k]

    return best, max_weight

print("Branch & Bound:", branch_and_bound(protein_seq, k))
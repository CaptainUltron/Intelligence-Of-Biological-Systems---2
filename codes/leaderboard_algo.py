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
'G':57, 'A':71, 'S':87, 'P':97, 'V':99,
    'T':101, 'C':103, 'I':113, 'L':113, 'N':114,
    'D':115, 'K':128, 'Q':128, 'E':129, 'M':131,
    'H':137, 'F':147, 'R':156, 'Y':163, 'W':186
}

k = 4

def leaderboard(protein, k, N):
    peptides = []

    for i in range(len(protein)-k+1):
        peptide = protein[i:i+k]
        weight = sum(AA_MW[a] for a in peptide)
        peptides.append((peptide, weight))

    peptides.sort(key=lambda x: x[1], reverse=True)
    return peptides[:N]

print("Leaderboard Top 5:")
for p in leaderboard(protein_seq, k, 5):
    print(p)
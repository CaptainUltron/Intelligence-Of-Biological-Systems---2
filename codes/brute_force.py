from Bio.Seq import Seq

text = ("GTTCCGAAAGGCTAGCGCTAGGCGCCAAGCGGCCGGTTTCCTTGGCGACGGAGAGCGCGGGAATTTTAGA"
    "TAGATTGTAATTGCGGCTGCGCGGCCGCTGCCCGTGCAGCCAGAGGATCCAGCACCTCTCTTGGGGCTTC"
    "TCCGTCCTCGGCGCTTGGAAGTACGGATCTTTTTTCTCGGAGAAAAGTTCACTGGAACTGGAAGAAATGG")

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
def brute_force(protein, k):
    max_weight = 0
    best = ""

    for i in range(len(protein)-k+1):
        peptide = protein[i:i+k]
        weight = sum(AA_MW[a] for a in peptide)

        if weight > max_weight:
            max_weight = weight
            best = peptide

    return best, max_weight

print("Brute Force:", brute_force(protein_seq, k))
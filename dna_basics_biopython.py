from Bio.Seq import Seq
#sample dna sequence 
dna= Seq("ATCGCTAGCTATCGCGTACCGGTTAA")
#calculate nucleotide frequency
print("DNA Sequence:", dna)
print("Length:", len(dna))
print ("\nNucleotide Frequency:")
for nucleotide in "ATGC":
     print(nucleotide,":" ,dna.count(nucleotide))
#calculate GC content
gc_content = ((dna.count("G") + dna.count("C")) /len(dna)) * 100
print("\nGC Content:", round(gc_content , 2) , "%")
#generate reverse complement
print("Reverse Complement:" , dna.reverse_complement())
#Transcription 
print("RNA Sequence:" , dna.transcribe())
#Translation
complete_length = (len(dna) // 3) * 3
coding_sequence = dna[:complete_length]

protein = coding_sequence.translate()

print("\nProtein Sequence:", protein)

if len(dna) % 3 != 0:
    print("Warning: Incomplete codon at the end of DNA sequence.")
    print("Untranslated nucleotides:", dna[complete_length:])
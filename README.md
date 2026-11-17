\# Python Bioinformatics Toolkit



A collection of Python-based tools developed to perform basic bioinformatics analysis of DNA, RNA, and protein sequences.



This project is part of my bioinformatics learning and portfolio development.



\## Features



\### 1. FASTA Analyzer

\- Reads FASTA sequence files.

\- Calculates sequence length.

\- Performs nucleotide composition analysis.

\- Calculates GC content.

\- Generates sequence statistics.



\### 2. Pairwise Sequence Aligner

\- Performs global sequence alignment.

\- Performs local sequence alignment.

\- Calculates alignment score.

\- Reports sequence identity and mismatches.



\### 3. DNA Analysis Using Biopython

\- Calculates nucleotide composition.

\- Calculates GC content.

\- Generates reverse complement.

\- Performs DNA transcription into RNA.

\- Translates DNA into protein.



\### 4. Protein Analyzer

\- Calculates protein length.

\- Determines amino acid composition.

\- Estimates molecular weight.

\- Calculates isoelectric point (pI).

\- Estimates net charge.

\- Calculates GRAVY score.

\- Calculates instability index.

\- Detects ambiguous and invalid amino acids.



\## Technologies Used



\- Python 3.11

\- Biopython

\- Python Standard Library



\## Project Structure



```text

Python\_\_toolkit/

├── fasta\_analyzer.py

├── pairwise\_aligner.py

├── dna\_basics\_biopython.py

├── protein\_analyzer.py

├── sample.fasta

├── test\_data/

│   └── TPH1.fasta

└── tests/

&#x20;   ├── \_\_init\_\_.py

&#x20;   └── test\_protein\_analyzer.py

```



\## Installation



Install Python 3.11.



Install Biopython:



```bash

py -3.11 -m pip install biopython

```



\## Running the Tools



Run the individual scripts using:



```bash

py -3.11 fasta\_analyzer.py sample.fasta

```



```bash

py -3.11 pairwise\_aligner.py

```



```bash

py -3.11 dna\_basics\_biopython.py

```



```bash

py -3.11 protein\_analyzer.py

```



Check each script's command-line options for its specific input requirements.



\## Testing



Run the protein analyzer unit tests:



```bash

py -3.11 -m unittest discover -s tests -v

```



\*\*Current test status:\*\* 10 tests passed successfully.



\## Learning Objectives



\- Understand biological sequence analysis.

\- Apply Python programming to bioinformatics.

\- Work with FASTA files.

\- Perform DNA and protein sequence analysis.

\- Develop and test reusable bioinformatics scripts.



\## Author



Sreelakshmi S

MSc Bioinformatics| Bioinformatics Enthusiast

https://github.com/Sreelakshmipp-jpg|sreelakshmispp26@gmail.com





\## Screenshots



\### Protein Analyzer

Example output from analyzing the human TPH1 protein sequence.



!\[Protein Analyzer Output](images/protein\_analyzer\_output.png)



\### FASTA Analyzer

Example output showing FASTA sequence analysis.



!\[FASTA Analyzer Output](images/fasta\_analyzer\_output.png)


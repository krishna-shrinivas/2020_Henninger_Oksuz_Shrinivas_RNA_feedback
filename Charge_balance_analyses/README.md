## Charge balance analyses

This code computes the charge balance ratio, as defined in the Methods section of the accompanying manuscript.

The python code is called in the following manner:

`` python charged_analysis.py --i Input/test_data.xlsx --p Input/test_protein.fasta --r Input/test_RNA.fasta --partition rna --conc 1000.0 --output test_figure ``



 To run the code, you'll need 3 input files (see Input/ for details on formatting)

1. test_data - 1 sheet each for partition ratio of protein/RNA. Each column is for a distinct concentration of RNA, and rows represent different droplets.
2. test_protein - Sequence of protein
3. test_rna - Sequence of RNA

The protein concentration is passed as input under the --conc flag in units of nM.

--partition flag chooses whether to plot partition of rna/protein.

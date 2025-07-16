This is a new approach for DFS in MFI generation, based on compressed transactions representation. To use the MFI-CT algorithm, copy the dataset you want to process into the same folder as the .exe. To test it, the program needs 2 parameters: the name of the dataset and the relative support threshold. For example: .\MFI_CT.exe chess.dat 0.72

The executable returns a file MFI_CTOut.txt with the list of maximal frequent itemsets. Each row contains the index of the MFI, followed by an arrow, the corresponding maximal frequent itemset,  a colon, and its associated absolute support. This is an example:

1 -> 46 5 36 : 2302 

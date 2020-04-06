A few files that might be helpful:

For FUll Length Assignment

input files:
CS7_FTO: chemical shift table (has N, H, CO', CA, CA' CB and CB') of the peaks  for the full length FTO
FASTA_FTO.tab: AA TASTA sequence of the FTO protein. 
fix_ass16.tab: fixed assignment. this file contains some possible assignments from the previous runs. 
				In addition, some residue types have been fixed based on the info from the selective labeling and distinct chemical shifts (A, G, S and T)
psipred_FTO.tab: secondary structure prediction based on the Fasta sequences
mars.inp (inside Run folder): input file for mars. this file lists the chemical shift table, fixed assignment,FASTA and the psipred files.

Output files (inside Run folder): There are many files but here are the ones that might be most useful
assignment_AA: list of the final assigned peaks (on the second column) for each residue (first column)
assignment_AAs: this is the list of the assignment and the confidence of multiple possible peaks. for example   LYS_4 has two possible peaks, 315 (11) (i.e. ) peak #315 with 11% probability  and peak #1023 with 58% probability.
assignment_PR: this lists the each peak to the most probable residue. since it lists all of the peaks, not everything is accurate, sometimes same peak is assigned to multiple residues.

there are other multiple files that you may look at for the info you might want to check.




for the isolated C-temrinal domain (CTD) assignment:

I was able to purify isolated CTD only and run the assignment on it. by comparing the chemical shifts, the peak labels on CTD are same as FUll length TROSY. some of the CTD peaks that don't match with full length chemical shifts, are numbered 5001 and above. 

CTD starts from residue # 296. Please note that first S for both Full length and CTD is from the TEV cleavage site of the fusion tag we use. this is present on both full length and isolated CTD.

since CTD is smaller, the assignment is more reliable from this domain separately.
I have compiled the excel file with the peak # and the assigned residues. (this doesn't list the unassigned residues). The N-terminal assignment (reidue # 1-295) is from the full-length assignment. from 296-475 is from the latest CTD asignment. please note that this file also includes the confidence (H for high, M for medium and L for the low).
Most of the High and Medium confidence residues are more or less reliable. most of the low confidence are less so (except for the ones that are assigned to a relatively large stretch of the sequence).

in addi


			
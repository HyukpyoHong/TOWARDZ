README
I. File descriptions 
	1. TOWARDZ.m
	 	- The efficient search code for translated networks with WR and ZD for given numbers of species and reactions using the pre- and post-translational filters introduced in the paper. 
	 	- In the code, key steps are described in comment lines.

	2. TOWARDZ_for_given.m
		- The efficient search code for translated network for a given network.
		- In the code, the two examples in the paper are given. 

	3. NetworkTranslation_naive.m
	 	- The naive search code for translated network with WR and ZD for given number of species and reactions. 
	 	- This code does not use the pre- and post-translational filters introduced in the paper. 

	4. CRN_translation_v1.m
		- Network translation code for a given network.
		- It generates all possitble translated network and check WR and ZD.
		- This code uses the post-translational filter (Theorem 3.2)
	
	5. CRN_translation_naive.m
		- Network translation code for a given network.
		- It generates all possitble translated network and check WR and ZD.
		- This code does not use the post-translational filter (Theorem 3.2)

	6. CRN_countlinkage.m
		- This function identifies the linkageclasses for a given network using a graph search algorithm.
		- This is necessary to check WR and to compute a deficiency.

	7. CRN_postivie_solution_optim.m
		- This function checks the existence of the vector with postive entries in the kernel of a given matrix.
		- This function utilizes an equivalent linear programming problem (Theorem 3.1 - iii).

	8. de2bi_hh.m
		- This function converts the number in the decimal system (i.e., base 10) to the binary number (i.e., base 2).


II. Figure generations
	1. The relevant data for Figures 1 and 2 can be generated using TOWARDZ.m by specifying the numbers of species and reactions.

	2. The examples in Figures 3 and 5 are manually chosen, and then using TOWARDZ_for_given.m, one can obtain the translated networks in Figures 4 and 6. The generalized CRN representation in Figure 7 can be derived manually.
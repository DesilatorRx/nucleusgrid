HFBTHO_IO
l392	declare the logical (CONST) variables set_neck_constraint_r
			    	      set_dfrag_constraint_r
				      set_csi_constraint_r
				      
l396	declare the real(CONST) variables    neckLag_r, neckRequested
		    	 	      dfragLag_r, dfragRequested
				      csiLag_r, csiRequested
	respectively the lagrange multiplier and the target value of the constraint



l419
l461
l487	read the file HFBTHO_NAMELIST and import the values of the variables declared above

l609
l641	report the variables above in thoout.dat



HFBTHO_MAIN
l296	activate the constraints and set their requested values


HFBTHO_VARIABLES
l161	declare the HO representation of the constraint operators gaussian_neck(:)
		       		      	     			  ddfrag(:)
								  ccsi(:)
l163	declare the real variables			  	  neckLag, neckValue, neckRequested(AGAIN?), mixing_neck=0.5									             dfragLag, dfragValue, dfragRequested(AGAIN?), mixing_dfrag=0.5
		    	 				  	  csiLag, csiValue, csiRequested(AGAIN?), mixing_csi=0.5
l236	declare logical vartiables				  neck_constraints=.False.
								  dfrag_constraints=.False.
								  csi_constraints=.False.

l282	set up namelist


HFBTHO_LARGE_SCALE 
l361	extend the size (l306) mpi vector with real types (11->13), to include the new constraint values
l384	extend the size (l306) mpi vector with logical types(10->12), to include the new constraint logicals

l486
l507	update the vector above when the broadcasted variables change

HFBTHO_FISSION
l130	declare constraint value at each iteration	Q_NECK
							DFRAG
							CSI
l953	print the characteristics of the fission fragments
l1302	routines DFRAGCALC() and CSICALC() calculate the value of constraints (fragment distance and mass asymmetry)
l1323	routines dfrag_computeField(ib) and csi_computeField(ib) calculate the matrix representation of the relative constraint in the K-block

HFBTHO_SOLVER
l289	count the number of constraints in the run, and save it in numberCons. Set the corresponding Lagrange multiplier at zero (WHY?)
l324	set multLambda(1:ncons_eff) at zero for the constraints not related to multipole moments (WHY?)
l818	print headings to thoout.dat
l981	print the characteristics of the run
l2323	check for io errors in HFBTHO_NAMELIST
l3030	allocate the broyden in and out vectors according to the number of MAXIMUM NUMBER of active constraints
l4983	add the constraining potential terms to the field
l5359	load the lagrange parameters for the constraints into broyden out vector
l5423	save the lagrange parameters parameters into broyden in vector
l6085 	contribution of the constraints to the half-projected (?) energies
l7000	calculate the expectation values of the constraints after convergence
l7400	build the vector of the deviation of the current constraint values from the requested values (cnsvec(numnerCons))
l7506	matrix of the constraints in the HO basis(neutron)
l7676	matrix of the constraints in the HO basis(proton)
l7977 	update the lagrange multiplier in order to get the requested constraint value in the current iteration





				     

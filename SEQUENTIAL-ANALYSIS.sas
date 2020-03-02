
ods html close;
ods html path="C:\Users\psioda\Documents\GitHub\SUGI-Session-4692\files" file="sequential-analysis.html";


%let nSims	         = 10000;
%let seed            = 1;

%let evidence_crit   = 0.975;
%let pse_eff_crit    = 0.975;
%let pse_fut_crit1   = 0.800;
%let pse_fut_crit2   = 0.100;

%let theta0          = 0.40;
%let theta1          = 0.67;

%let nMin            = 20;
%let nMax            = 60;
%let nBy             =  2;


%let h_y1            = 162;
%let h_y0            =  80;
%let h_n             = 242;


%let outcome_dist_mean       = 2.00;
%let outcome_dist_sd         = 0.10;
%let poisson_process_mean    = 0.55; 


data start;
 start = time();
run;

proc IML;
	call streaminit(&seed.);
	call randseed(&seed.);

	start BetaBin (y,n) global(alpha,beta);
		return exp(lcomb(n,y) + logbeta(y+alpha,n-y+beta) - logBeta(alpha,beta));
	finish;

	rho = 3;

	prm1_enr = &poisson_process_mean.;
	prm1_asc = &outcome_dist_mean.;
	prm2_asc = &outcome_dist_sd.;

	theta0 = &theta0.;
	theta1 = &theta1.;
	thetam = 0.5*theta0+0.5*theta1;


	nSims = &nSims.;

	nBy   = &nBy.;
	maxN  = &nMax.;
	minN  = &nMin.;

	evidence_crit = &evidence_crit.;
	pse_eff_crit  = &pse_eff_crit.;
	pse_fut_crit1 = &pse_fut_crit1.;
    pse_fut_crit2 = &pse_fut_crit2.;


	h_y1   = &h_y1.;
	h_y0   = &h_y0.;
	h_n    = &h_n.;

	/** Precomputation of borrowing parameters **/

	c0_matrix = J(maxN,maxN+1,0);
	c1_matrix = J(maxN,maxN+1,0);
	c2_matrix = J(maxN,maxN+1,0);
	a0_matrix = J(maxN,maxN+1,0);

	do n  = 1 to maxN;
	do y1 = 0 to n;
		r              = n;
		c              = y1+1;

		alpha          = h_y1 + 0.5;
		beta           = h_y0 + 0.5;
		PredProb       = betaBin(t(do(0,n,1)),n);
		obsProb        = betaBin(y1,n); 
		c1_matrix[r,c] = ((PredProb<=obsProb)#PredProb)[+];

		PredProb       = pdf('binomial',t(do(0,n,1)),theta0,n);
		obsProb        = pdf('binomial',y1,theta0,n);
		c2_matrix[r,c] = ((PredProb<=obsProb)#PredProb)[+];

		c0_matrix      = c1_matrix - c2_matrix;
		c0_matrix      = c0_matrix >< 1;
		c0_matrix      = c0_matrix <> 0;

		a0_matrix[r,c] = min(1,c0_matrix[r,c]*rho*n);
	end;
	end;


	pp_skept_matrix_borrow   = J(maxN,maxN+1,0);
	pp_skept_matrix_noborrow = J(maxN,maxN+1,0);
	pp_enthu_matrix          = J(maxN,maxN+1,0);

	/** Precomputation of posterior probabilities **/
	do n  = 1 to maxN;
	do y1 = 0 to n;
	    y0                       = n-y1;
		r                        = n;
		c                        = y1+1;

		alpha                         = 5.830 + y1 + a0_matrix[r,c]*h_y1;
		beta                          = 8.245 + y0 + a0_matrix[r,c]*h_y0;
		pp_skept_matrix_borrow[r,c]   = sdf('beta',theta0,alpha,beta);

		alpha                         = 5.830 + y1;
		beta                          = 8.245 + y0;
		pp_skept_matrix_noborrow[r,c] = sdf('beta',theta0,alpha,beta);

		alpha                         = 9.790 + y1;
		beta                          = 5.329 + y0;
		pp_enthu_matrix[r,c]          = cdf('beta',thetam,alpha,beta);
	end;
	end;



    tpv = do(0.40,0.76,0.02);
	*tpv = 0.40||0.67;

	do borrow = 0 to 1;
    do t = 1 to ncol(tpv);
	true_pi = tpv[t];
	results = J(nSims,26,0);
 	do sim = 1 to nSims;


		stop_trial	     = 0;     /** Indicator for early stoppage **/
		final_analysis   = 0;     /** Indicator for final analysis  **/
		n			     = 0;     /** Number of patients currently enrolled **/
		analysis         = 0;     /** Number of analyses performed **/
		time_of_analysis = {0 0}; /** vector for time of the analysis **/
		Eff              = {0 0}; /** vector for efficacy criterion [interim/final] **/
		Fut1             = {0 0}; /** vector for futility criterion 1 [interim/final] **/
		Fut2             = {0 0}; /** vector for futility criterion 2 [interim/final] **/
		Fut3             = {0 0}; /** vector for futility criterion   [interim/final] **/

		nInt           = .;
		y1Int 		   = .;
		y0Int          = .;
		nMissInt       = .;


	/** Generate the complete hypothetical dataset **/
	r    = J(maxN,1,0);
	w    = J(maxN,1,0);
	y    = J(maxN,1,.);

	/** enrollment times via a Poisson Process **/
	call randgen(r,'exponential',prm1_enr);
	r = cusum(r);

	/** outcome ascertainment times via a normal distribution **/
	call randgen(w,'normal',prm1_asc,prm2_asc);

	/** calculated time from enrollment to outcome ascertainment **/ 
	e = r + w;

	/** simulate response outcome **/
	call randgen(y,'bernoulli',true_pi);

	/** Sort dataset by ascending order of outcome ascertainment **/
	dat = r||w||e||y;
	call sort(dat,3);

	r = dat[,1];
	w = dat[,2];
	e = dat[,3];
	y = dat[,4];

	free dat;


		/** Simulate the sequentially monitored trial **/
		do until(stop_trial=1 & final_analysis=1);

			/** if the trial has not yet provided substantial evidence **/
		    if stop_trial = 0 then do;
				/** increment number of ascenertained outcomes **/;
				n = n + nBy; 

				/** increment the anaysis number **/
				if n >= minN then analysis = analysis + 1;

				/** identify time of current analysis **/
				time_of_analysis[1] = e[n];

				/** accumulate the outcome data into sufficient statistics **/
				y0 = sum((y=0)[1:n]);
				y1 = sum((y=1)[1:n]);

				/** determine how many subjects are currently already enrolled **/
			    nMiss = nrow(y[loc(r<time_of_analysis[1])])-n;

				/** determine how many subjects coudl still be enrolled **/
				nLeft = max(0,maxN - n);
			end;
			else do; /** if the trial has yet provided substantial evidence **/
				final_analysis = 1;
				nInt           = n;
				y1Int 		   = y1;
				y0Int          = y0;
				nMissInt       = nMiss;

				/** if the trial was NOT stopped for futility **/
				if FUT3[1] = 0 then do;
					time_of_analysis[2]     = max(e[loc(r<time_of_analysis[1])]);
					y                       = y[loc(r<time_of_analysis[1])];
					n                       = nrow(y);
					y0 						= sum((y=0));
					y1 						= sum((y=1));
					nMiss                   = 0;								
				end;
				else do; /** if the trial was stopped for futility **/
					time_of_analysis[2]     = time_of_analysis[1];
				end;									  
			end;

			if 	(n>= minN) then do;

				/** compute borrowing parameter **/
				row = n;
				col = y1+1;
				c0  = c0_matrix[row,col];
				a0  = a0_matrix[row,col]*borrow;

				/** compute posterior probability for efficacy **/
				postprob_skept = pp_skept_matrix_borrow[row,col]*(borrow=1)
							   + pp_skept_matrix_noborrow[row,col]*(borrow=0);


				/** compute posterior probability for futility **/
				postprob_enthu = pp_enthu_matrix[row,col];


				/** compute PSSE for skeptical & enthusiastic priors **/
				if nMiss = 0 then do; PSE_skept = postprob_skept; PSE_enthu = postprob_enthu; end;
				else do;
					
					nMissTot = n + nMiss;
					nLeftTot = n + nLeft;

					/** Skeptical prior **/

						/** Compute predictive distribution for missing outcomes **/
						alpha        = 5.830 + y1 + a0*h_y1;
						beta         = 8.245 + y0 + a0*h_y0;
						MissPredProb = betaBin(t(do(0,nMiss,1)),nMiss);

						/** Compute the probability of sustained evidence **/
			 			PSE_skept = 0;
						do	y1Miss = 0 to nMiss;
						  	y0Miss = nMiss - y1Miss;

							y1MissTot = y1 + y1Miss;
			                y0MissTot = nMissTot - y1MissTot;

							/** Compute final borrowing parameter value **/
							row  = nMissTot;
							col  = y1MissTot+1;
				     		pps = pp_skept_matrix_borrow[row,col]*(borrow=1)
							    + pp_skept_matrix_noborrow[row,col]*(borrow=0);

							PSE_skept = PSE_skept + (pps>evidence_crit)*MissPredProb[y1Miss+1];
						end;

						/** Compute predictive distribution for possible remaining outcomes **/
						alpha        = 5.830 + y1 + a0*h_y1;
						beta         = 8.245 + y0 + a0*h_y0;
						LeftPredProb = betaBin(t(do(0,nLeft,1)),nLeft);

						/** Compute the probability of sustained evidence **/
			 			PSE_skept2 = 0;
						do	y1Left = 0 to nLeft;
						  	y0Left = nLeft - y1Left;

							y1LeftTot = y1 + y1Left;
			                y0LeftTot = nLeftTot - y1LeftTot;

							row  = nLeftTot;
							col  = y1LeftTot+1;
				     		pps2 = pp_skept_matrix_borrow[row,col]*(borrow=1)
							     + pp_skept_matrix_noborrow[row,col]*(borrow=0);

							PSE_skept2 = PSE_skept2 + (pps2>evidence_crit)*LeftPredProb[y1Left+1];
						end;

					/** Enthusiastic prior **/

						/** Compute predictive distribution for missing outcomes **/
						alpha     = 9.790 + y1;
						beta      = 5.329 + y0;
						MissPredProb = betaBin(t(do(0,nMiss,1)),nMiss);

						/** Compute the probability of sustained evidence **/
						PSE_enthu = 0;
						do	y1Miss = 0 to nMiss;
						  	y0Miss = nMiss - y1Miss;

							y1MissTot = y1 + y1Miss;
			                y0MissTot = nMissTot - y1MissTot;


							row  = nMissTot;
							col  = y1MissTot+1;
							ppe  = pp_enthu_matrix[row,col];

							/** Compute beta distribution paramters for predicted final data and PSSE (Skeptical) **/
							PSE_enthu = PSE_enthu + (ppe>evidence_crit)*MissPredProb[y1Miss+1];
						end;
				end;


				/** Stop for efficacy **/
				if	(PSE_skept>=pse_eff_crit  )   then do; 
					stop_trial = 1; 
					Eff[1+final_analysis] = 1; 
				end; 

				/** Stop for futility **/
				if	(PSE_enthu >= pse_fut_crit1  )   then do; 
					stop_trial = 1; 
					FUT1[1+final_analysis] = 1; 
				end; 			

				/** Stop for futility **/
				if  (PSE_skept2 < pse_fut_crit2 )   then do; 
					stop_trial = 1; 
					FUT2[1+final_analysis] = 1; 
				end; 	

				FUT3 = (FUT1 + FUT2)><1;

				/** Stop for maximum sample size **/				
				if ( n>= maxN )   then stop_trial = 1;

			end;
		end;


		


        if final_analysis = 1 then 
          results[sim,] =  borrow||true_pi||analysis||nInt||y1Int||y0Int||nMissInt||n||y1||y0||(y1/n)||a0||c0||
                           postprob_skept||PSE_skept||PSE_skept2||postprob_enthu||PSE_enthu||time_of_analysis||Eff||FUT3||FUT1[2]||FUT2[2] ;

	end;
		if true_pi = tpv[1] & borrow = 0 then all_results = results;
		else all_results=all_results//results;
	end;
	end;




 	colnames={	"borrow" "true_pi" "I" "nInt" "y1Int" "y0Int" "nMissInt" "nFin" "y1Fin" "y0Fin" "piHatFin" "a0" "c0"
				"postprob_skept" "PSE_skept" "PSE_skept2" "postprob_enthu" "PSE_enthu"  
				"durInt"  "durFin" "EffInt" "EffFin" "FutInt" "FutFin" "FutCrit1" "FutCrit2" };

 create simulation_results 
   from all_results[c=colnames];
   append from all_results;
 close simulation_results;


quit;

data stop;
 set start;
 stop = (time()-start)/60;
run;

proc means data = simulation_results noprint nway; 
 class borrow true_pi  / missing;
 var I nInt y1Int y0Int nMissInt nFin y1Fin y0Fin piHatFin a0 c0
     postprob_skept PSE_skept PSE_skept2 postprob_enthu PSE_enthu durInt durFin EffInt EffFin FutInt FutFin FutCrit1 FutCrit2;
 output out = simulation_results2
   mean = I nInt y1Int y0Int nMissInt nFin y1Fin y0Fin piHatFin a0 c0
      postprob_skept PSE_skept PSE_skept2 postprob_enthu PSE_enthu durInt durFin EffInt EffFin FutInt FutFin FutCrit1 FutCrit2;
run;



ods html close;
ods html file="C:\Users\psioda\Documents\GitHub\SUGI-Session-4692\files\sequential-analysis.html";

	proc print data = simulation_results2;
	by borrow;
	 var borrow true_pi _freq_ I nInt y1Int y0Int nMissInt nFin y1Fin y0Fin piHatFin a0 c0
	      postprob_skept PSE_skept PSE_skept2 postprob_enthu PSE_enthu durInt durFin EffInt EffFin FutInt FutFin FutCrit1 FutCrit2;
	run;

ods html close;


ods graphics / noborder;
ods pdf file="C:\Users\psioda\Documents\GitHub\SUGI-Session-4692\files\sequential-analysis.pdf" startpage=no dpi=350;
	proc sgplot data = simulation_results2;

	 format true_pi 6.2;
	 series x=true_pi y=EffFin / group=borrow markers markerattrs=(symbol=circlefilled);
	 refline 0.025/ axis=y lineattrs=(pattern=2);
	 refline 0.67 / axis=x lineattrs=(pattern=2);
	 yaxis values=(0 to 1 by 0.1) grid label='Probability of Substantial Evidence';
	 xaxis label="True Response Probability";
	run;

	proc sgplot data = simulation_results;
	 vbox nFin / category=true_pi group=borrow nooutliers meanattrs=(symbol=diamondFilled);
	 yaxis grid label='Sample Size';
	 xaxis label="True Response Probability";

	run;
	ods pdf startpage=now;

	proc freq data = simulation_results noprint;
	 by borrow true_pi;
	 table FutCrit1 / out =simulation_results3(where=(FutCrit1=1));
	run;

	proc freq data = simulation_results noprint;
	 by borrow true_pi;
	 table FutCrit2 / out =simulation_results4(where=(FutCrit2=1));
	run;

	proc sgplot data = simulation_results3;
	 series x=true_pi y=percent / group=borrow markers markerattrs=(symbol=circlefilled);
	run;

	proc sgplot data = simulation_results4;
	 series x=true_pi y=percent / group=borrow markers markerattrs=(symbol=circlefilled);
	run;

ods pdf close;

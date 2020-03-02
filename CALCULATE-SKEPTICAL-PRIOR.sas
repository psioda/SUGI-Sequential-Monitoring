
option nodate nonumber;
title;
ods html close;
ods pdf file="C:\Users\psioda\Documents\GitHub\SUGI-Session-4692\files\skeptical-prior.pdf" startpage=no;

proc iml;
	theta0        = 0.40;
	theta1        = 0.67;
	evidence_crit = 0.975;

	start calc_beta(alpha,mode);
		return ((1-mode)*alpha + 2*mode - 1)/mode; 
	finish;

	delta = 0.01;
	stop  = 0;
	do alpha = 1+delta to 10 by delta while(stop=0);
		beta = calc_beta(alpha,theta0);
		prob = cdf('beta',theta1,alpha,beta);
			if (beta>1) & 
	           ((evidence_crit-0.00005)<prob) & 
	           (prob<(evidence_crit+0.00005)) then do;
				print alpha beta prob evidence_crit;
				stop=1;
			end;
	end;


	pi  = t(do(0.001,0.999,0.001));
	pdf = pdf('beta',pi,alpha,beta);

	plotdata = pi||pdf||(pi>theta1);
	create plotdata from plotdata[c={"pi" "pdf" "ind"}];
	 append from plotdata;
	close plotdata;
quit;


data plotdata;
 set plotdata;
  low = 0;
   if ind = 0 then pi0 = pi;
   if ind = 1 then pi1 = pi;
run;


proc sgplot data = plotdata noautolegend;
 where pdf > 1e-3;
 band x=pi0 upper=pdf lower=low / lineattrs=(color=blue) fillattrs=(color=blue);
 band x=pi1 upper=pdf lower=low / lineattrs=(color=red)  fillattrs=(color=red);
 yaxis label='Density Value';
 xaxis label='Parameter Value';
run;


ods pdf close;


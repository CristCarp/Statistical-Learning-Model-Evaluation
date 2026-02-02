proc iml;

/* Si l'on veut Y indépendant des X */
start independance(n_samples);
    /* Y est juste du bruit pur, indépendant de toute variable X */
    new_y = normal(j(n_samples, 1, 0)); 
    return(new_y);
finish;
store module=(Independance);

/* Pour simuler de la multicolinéarité avec Toeplitz, Iman-Conover permet de conserver les outliers (pas utile dans notre cas) */
start ImanConoverToeplitz(X, rho); /* Le rho est le degré de corrélation de Toeplitz */
   p=ncol(X);
   /* Toeplitz */
   C = j(p,p,0);
   do i = 1 to p;
      do j = 1 to p;
         C[i,j] = rho##abs(i-j);
      end;
   end;
   /* Iman-Conover */
   N = nrow(X);
   S = J(N, ncol(X));
   do i = 1 to ncol(X);
      ranks = ranktie(X[,i], "mean");
      S[,i] = quantile("Normal", ranks/(N+1)); 
   end;
   CS = corr(S);
   Q = root(CS);
   P = root(C); 
   T = solve(Q,P);
   Y = S*T;
   W = X;
   do i = 1 to ncol(Y);
      rank = rank(Y[,i]);
      tmp = W[,i]; call sort(tmp);
      W[,i] = tmp[rank];
   end;
   return(W);
finish;
store module=(ImanConoverToeplitz);

/* Pour simuler une rupture au point break*n, size va définir la variation de nos beta */
start Rupture(break, size, Y, X, beta, eps);
	beta2=j(nrow(beta),1,0);
	do i=1 to nrow(beta);
		beta2[i]=beta[i]+size*(2*rand("uniform")-1);
	end;
	idx=int(break*nrow(Y));
	do i=idx to nrow(Y);
		Y[i]=X[i,5]*beta2[1]+X[i,6]*beta2[2]+X[i,20]*beta2[3]+X[i,30]*beta2[4]+X[i,40]*beta2[5]+eps[i];
	end;
	return(Y);
finish;
store module=(Rupture);

/* Simule des outliers sur les résidus de fréquence freq ce qui va décaler leur moyenne à droite de size*/
start Extremes(eps, freq, size);
	u=ranuni(j(1,nrow(eps),1));
	do i=1 to nrow(eps);
		if u[i]<freq then do;
			eps[i]=eps[i]+size;
		end;
	end;
	return(eps);
finish;
store module=(Extremes);

%macro MCglm(n, MC, select, choose, stop, inde, rupture, multi, extreme);


/* ===================== Clean ===================== */
proc datasets lib=work nolist;
    delete bigtrain /*bigtest*/;

quit;



proc iml;
    load module=(Extremes Rupture ImanConoverToeplitz Independance);
	/* On décide de Générer un X de taille n*MC puis on le découpera en MC bloc n */
    n = &n.*&MC.;
	p = 50;
	
    beta = j(p,1,0);
    beta[1]=0.9; beta[2]=-1.0; beta[3]=1.2; beta[4]=-0.8; beta[5]=1.1;

    X = normal(j(n,p,0));
    eps = normal(j(n,1,0))*0.25;

	/* ===================== Multvarié ===================== */
	/* donc X suit une loi normal(m,var) avec m un vecteur et var une matrice diagonale */
	m = 2*ranuni(j(1,p,1))-1;
	var = j(p,p,0);
	do i = 1 to p;
		var[i,i] = 3*ranuni(0)+1;
	end;
	X = X*sqrt(var) + m;

    /* ===================== Valeurs extrêmes ===================== */
    %if &extreme.=oui_extreme %then %do;
        eps = Extremes(eps, 0.05, 4);
    %end;

    /* ===================== Multi-colinéarité ===================== */
    %if &multi.=oui_multi %then %do;
        X = ImanConoverToeplitz(X, 0.8);
    %end;

    /* Modèle */
    Y = X[,5]*beta[1] + X[,6]*beta[2] + X[,20]*beta[3] + X[,30]*beta[4] + X[,40]*beta[5] + eps;

    /* ===================== Indépendance ===================== */
    %if &inde.=oui_inde %then %do;
        Y = Independance(n);
    %end;

    /* ===================== Rupture ===================== */
    %if &rupture.=oui_rupture %then %do;
        Y = Rupture(0.5, 1, Y, X, beta, eps);
    %end;
	
	train = {};
	test = {};
	%do r = 1 %to &MC.;

        /* ===================== Split train/test ===================== */
		/*u = ranuni(j(1,&n.,0));
        idxTrain = loc(u < 2/3)+&n.*(&r.-1);
        idxTest  = loc(u >= 2/3)+&n.*(&r.-1);

        train = train // (j(ncol(idxTrain),1,&r.) || Y[idxTrain] || X[idxTrain,]);
        test  = test // (j(ncol(idxTest),1,&r.) || Y[idxTest] || X[idxTest,]);*/
		train = train // (j(&n.,1,&r.) || Y[1+&n.*(&r.-1):&n.*&r.] || X[1+&n.*(&r.-1):&n.*&r.,]);

	%end;

	colnames = "rep" ||"Y" || ( "X1":"X50");

    create bigtrain from train[colname=colnames];
    append from train;
    close bigtrain;

    /*create bigtest from test[colname=colnames];
    append from test;
    close bigtest;*/


quit;

ods exclude all;
ods output ParameterEstimates = ParEst;

proc glmselect data=bigtrain /*testdata=bigtest*/ plots=all;
    by rep;
    model Y = X1-X50
        / selection=&select choose=&choose stop=&stop;


run;
ods exclude none;
/* ===================== MODELE CHOISI ===================== */


data ModeleChoisi;
    set ParEst;
    where Parameter ne "Intercept";
run;

/* ===================== FP FN ===================== */


proc iml;

	vrai={X5 X6 X20 X30 X40};
    use ModeleChoisi; read all var {"rep"} into rep; close;
	use ModeleChoisi; read all var {"Parameter"} into ChoisiMC; close;

	/* Si ModeleChoisi est vide */
	if nrow(rep)=0 then do;
    	rep = j(1,1,.);
    	ChoisiMC = j(1,1,"");
	end;

    results = j(&MC., 3, .);      /* rep, FP, FN */
	
    do i = 1 to &MC.;
        idx = loc(rep=i);
		choisi = {};
		if ncol(idx)^=0 then do;
        	choisi = ChoisiMC[idx]`;
		end;

        /* FP = choisi \ vrai */
        FP = setdif(choisi, vrai);

        /* FN = vrai \ choisi */
        FN = setdif(vrai, choisi);

        results[i,] = i || ncol(FP) || ncol(FN);

    end;

    create Indicateurs from results[colname={"rep" "FP" "FN"}];
    append from results;
    close Indicateurs;


quit;

/* ===================== INDICATEURS ===================== */

data Indicateurs;
    set Indicateurs;
	methode  = "&select.";
	choose   = "&choose.";
	stop     = "&stop.";
    Exact    = (FP=0 and FN=0);
    Overfit  = (FP>0 and FN=0);
    Underfit = (FP=0 and FN>0);
    Mixed    = (FP>0 and FN>0);

run;


proc sql;
    create table ResumeIndicateurs as
    select 
        methode,
		choose,
        stop,
        mean(Exact)    as ExactM,
        mean(Overfit)  as OverfitM,
        mean(Underfit) as UnderfitM,
        mean(Mixed)    as MixedM
    from Indicateurs
    group by methode, choose, stop;

quit;


proc append base=Indicateurs_All data=ResumeIndicateurs force; run;

%mend;

proc datasets lib=work nolist;
	delete Indicateurs_All;

quit;

title "Dépendance linéaire";
%MCglm(200,300, elasticnet, cv, cv, non_inde, non_rupture, non_multi, oui_extreme);
proc print data=Indicateurs_All; run;

proc datasets lib=work nolist;
	delete Indicateurs_All;

quit;

title "Multi_colinéarité";
%MCglm(200,300, lasso, cv, cv, non_inde, non_rupture, oui_multi, non_extreme);
proc print data=Indicateurs_All; run;

proc datasets lib=work nolist;
	delete Indicateurs_All;
quit;

title "Valeurs Extremes";
%MCglm(200,300, lasso, cv, cv, non_inde, non_rupture, non_multi, oui_extreme);
proc print data=Indicateurs_All; run;

proc datasets lib=work nolist;
	delete Indicateurs_All;
quit;

title "Rupture";
%MCglm(200,300, lasso, cv, cv, non_inde, oui_rupture, non_multi, non_extreme);
proc print data=Indicateurs_All; run;

proc datasets lib=work nolist;
	delete Indicateurs_All;
quit;

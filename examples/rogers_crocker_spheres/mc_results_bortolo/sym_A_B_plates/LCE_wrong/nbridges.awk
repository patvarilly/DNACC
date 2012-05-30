
{

# in nanometer unit 
    z[NR]=$1*5 ;
    Pe[NR]=$2/5 ;
    
}END{

    pi=3.1415926535897932385 ;
    rho0 = 0.602214 ;
#denities (number of strands/nm^2)
#    na = 4800/(4*pi*550*550);
#    nb = 4200/(4*pi*550*550);
# for the current case
    na = 4500/(4*pi*550*550);
    nb = 4500/(4*pi*550*550);

    for(iDG=1;iDG<30;iDG++){

	DG = (iDG-1)/2 ;
	Govrho=exp(DG)/rho0 ;

	for(i=1;i<=NR;i++){
	    cA=na*Pe[i]*Govrho;
	    cB=nb*Pe[NR-i+1]*Govrho;
	    cTot=cA+cB;
	    fint[i]=1+cTot-sqrt( (1+cTot)^2-4*cA*cB );
	    fint[i]=fint[i]/(2*Govrho);
	}

# integer over z
	Int=0

#	for(i=1;i<=NR-1;i++){
#	    Int=Int+(z[i+1]-z[i])*(fint[i]+fint[i+1])/2
#	}

	dz=z[2]-z[1]
	for(i=1;i<=NR;i++){
#	    Int=Int+dz*(fint[i]+fint[i+1])/2
	    Int=Int+dz*fint[i]
	}

#	Int = Int / ( x[NR] - x[1] )
#	{print x[NR]+x[2]-2*x[1], Int}


#                          to compare with the simulations
	{print -DG, Int*5*5*129.98*129.98}

    }

}


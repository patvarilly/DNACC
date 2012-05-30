BEGIN{

    NDG=6;
    for(i=1;i<=NDG;i++){
	V[i] = 0;
	NV[i]= 0;
    }

}{

    ind = NR%NDG ;
    if(ind == 0){
	ind = NDG ;
    }
    x[ind]=$1;
    V[ind]=V[ind]+$2;  
    NV[ind]=NV[ind]+1; 

}END{

    for(ind=1; ind<=NDG; ind++){
	{print x[ind], V[ind]/NV[ind]};
    }

}
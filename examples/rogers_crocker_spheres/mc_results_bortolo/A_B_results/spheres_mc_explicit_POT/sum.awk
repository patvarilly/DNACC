{
    x[NR] = $1;
    y[NR] = $2;
}END{

    for(i=1;i<=NR/2;i++){
	ip= i+NR/2;
	{print x[i], y[i]+y[ip]}
    }

 }

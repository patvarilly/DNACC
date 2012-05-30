	Program DNAPairPotentials     

C****************************************************************************************************
C****************************************************************************************************
C*		Copyright: Mirjam E. Leunissen, FOM Institute AMOLF, The Netherlands - 11 October 2010		*
C*		No restrictions apply to the usage of this program provided that users cite:     			*
C*		M.E. Leunissen and D. Frenkel, J. Chem. Phys. 134, 084702 (2011), doi: 10.1063/1.3557794	*
C****************************************************************************************************
C****************************************************************************************************

C****************************************************************************************************	
C		This program generates the DNA-mediated pair interaction potentials for planar surfaces		*
C		and spherical particles, functionalized with stiff DNA constructs of Ldna = 20 nm long		*
C		at the user-defined average strand spacing or surface coverage (S), the DNA binding			*
C		strength (G) and the particle radius (R).													*
C		The derivation of the analytical expressions and the working range of these pair			*
C		interactions can be found in Section IV and the rest of:									*
C		M.E. Leunissen and D. Frenkel, J. Chem. Phys. 134, 084702 (2011), doi: 10.1063/1.3557794	*
C****************************************************************************************************	

C	--- Declare variables

	Implicit None

	Integer NGLeX
	Double Precision GLeX(10),GLeW(10)
	  
	Integer Mode,TotalPoints
	Double Precision hmin,hmax,dh,S,G,R
	  
	Double Precision Pi,Rho,PlateFitConstants(8),SphereFitConstants(3)
	  
	Integer P,N,M,W
	Double Precision Points(10000),h,h2,rperp1,rperp2,rperpint,rperp,frperp,
     &	B1,B2,B3,B,Attraction(10000),Repulsion(10000),Total(10000),Fplate,
     &	FplateInt,hplate,hplate2,Attractionh1
	 
		   	   
C	--- Gauss-Legendre quadrature
	
	NGleX	=	10
       
	GLeX	=	(/	-0.9739065285d0, -0.8650633667d0, -0.6794095683d0,
     &				-0.4333953941d0, -0.1488743390d0,  0.1488743390d0, 
     &				 0.4333953941d0,  0.6794095683d0,  0.8650633667d0, 0.9739065285d0 /)        
                
	GLeW	=	(/	 0.0666713443d0, 0.1494513492d0, 0.2190863625d0, 
     &				 0.2692667193d0, 0.2955242247d0, 0.2955242247d0,
     &				 0.2692667193d0, 0.2190863625d0, 0.1494513492d0, 0.0666713443d0 /)

  	  
C	--- Read input data and check the number of data points

	Open	(15,FILE='input',FORM='FORMATTED')	   
	Read	(15, *) 					
	Read	(15, *)
	Read	(15, *) Mode, hmin, hmax, dh, S, G, R
	Close	(15)

	TotalPoints	=	INT((hmax-hmin)/dh)	
		
	If(TotalPoints.Gt.10000)Then
		Write (*,*) '!!!!!**Too many data points, increase dH**!!!!!'
	Else


C	---- Define constants        
	  
	Pi	=	3.1415926535897932d0
C	Rho: standard number density for Ldna = 20 nm
	Rho	=	4817.7d0 
	  
	PlateFitConstants	=	
     &(/	1.34d-5,	6.27d-6,	-7.03d-6,
     &		5.90d-6,	1.24d-4,	-1.05d-8,
     &		3.49d-13,	-3.68d-18	/)

	SphereFitConstants	=	(/ -29.15d0, 4.18d-5, -6.78d-10 /)

	Open (61,FILE='potential')


C	--- Calculate the attractive and repulsive parts of the plate-plate interaction

	If(Mode.Eq.1)Then
	
		Write(*,*)
		Write(*,*) '*********************************************************'
		Write(*,*) 'Calculating the DNA-mediated plate-plate interaction...'
		Write(*,*) 'S =', S
		Write(*,*) 'G =', G
		Write(*,*) 'R =', R
	
		h	=	hmin
		P	=	0

		Do While(h.Le.hmax)
		
			P			=	P + 1
			Points(P)	=	h
			h2			=	h*h
		
			If(h.Ge.1.0d0)Then		
		
				Attraction(P)	=	(1.0d0/S**2)*(-DLOG(1.0d0+(DEXP(-G)*(4.0d0/(3.0d0*S**2*Rho)))*(1.0d0-0.25d0*h2)**1.5) +		
     &			PlateFitConstants(1)*(1.0d0/S)*DEXP(-G)*(2.0d0-h) + PlateFitConstants(2)*(1.0d0/(S**3))*DEXP(-G)*(2.0d0-h) +
     &			PlateFitConstants(3)*(1.0d0/S)*DEXP(-G)*(2.0d0-h)**2 + PlateFitConstants(4)*(1.0d0/(S**3))*DEXP(-G)*(2.0d0-h)**2)
				Repulsion(P)	=	0.0d0

			Else
				
				B1			=	(Pi/3.0d0)*((4.0d0-h2)**1.5)
				rperp1		=	DSQRT(2.0d0-h2-2.0d0*DSQRT(1.0d0-h2))
				rperp2		=	DSQRT(2.0d0-h2+2.0d0*DSQRT(1.0d0-h2))
				rperpint	=	0.0d0	  
		
				Do N = 1, NGleX
		
					rperp		=	0.5d0*(rperp2-rperp1)*GLeX(N)+0.5d0*(rperp1+rperp2)												
					frperp		=	-4.0d0*rperp*DSQRT(1.0d0-0.25d0*(h2+rperp*rperp))*
     &								DACOS((h*DSQRT(h2+rperp*rperp))/(2.0d0*rperp*DSQRT(1.0d0-0.25d0*(h2+rperp*rperp))))	 
					rperpint	=	rperpint + GLeW(N)*frperp
			
				Enddo			
		 
				B2				=	0.5d0*(rperp2-rperp1)*rperpint			
				Attraction(P)	=	(1.0d0/S**2)*(-DLOG(1.0d0+(1.0d0/(2.0d0*Pi*Rho*S**2*h2))*DEXP(-G)*(B1+B2)) +
     &								PlateFitConstants(5)*(1.0d0/S**2)*DEXP(-G)*(1.0d0-h)**2 +
     &								PlateFitConstants(6)*(1.0d0/S**4)*((DEXP(-G))**2)*(1.0d0-h)**2 +
     &								PlateFitConstants(7)*(1.0d0/S**6)*((DEXP(-G))**3)*(1.0d0-h)**2 +
     &								PlateFitConstants(8)*(1.0d0/S**8)*((DEXP(-G))**4)*(1.0d0-h)**2)	
				Repulsion(P)	=	-(1.0d0/S**2)*2.0d0*DLOG(h)
		
			Endif
		
			h	=	h + dh
		
		Enddo


C	--- Calculate the total interaction and write to file

		Write (61,*) '# Column:'
		Write (61,*) '# 1 - Plate-plate separation [Ldna]'
		Write (61,*) '# 2 - Repulsive contribution [kT]'
		Write (61,*) '# 3 - Attractive contribution [kT]'
		Write (61,*) '# 4 - Total DNA-mediated interaction [kT]'
		Write (61,*) '############################################'

			Do W = 1, P
				Total(W)	=	Attraction(W) + Repulsion(W)
				Write (61,9900)	Points(W), Repulsion(W), Attraction(W), Total(W) 
			Enddo	
		
	Else
	
C	--- Calculate the sphere-sphere interaction
	
		Write(*,*)
		Write(*,*) '*********************************************************'
		Write(*,*) 'Calculating the DNA-mediated sphere-sphere interaction...'
		Write(*,*) 'S =', S
		Write(*,*) 'G =', G
		Write(*,*) 'R =', R


C	--- Calculate the interaction at h = 1.0 (used later in the calculation of all h < 1)	
	
		h			=	1.0d0
		FplateInt	=	0.0d0
		
      	  Do N = 1, NGLeX		
		
			hplate		=	0.5d0*(2.0d0-h)*GLeX(N)+0.5d0*(2.0d0+h)
			Fplate		=	DLOG(1.0d0+(DEXP(-G)*(4.0d0/(3.0d0*S**2*Rho)))*(1.0d0-0.25d0*hplate*hplate)**1.5)
			FplateInt	=	FplateInt + GLeW(N)*Fplate		
	
          Enddo
	  
		Attractionh1	=	-Pi*R*0.5d0*(2.0d0-h)*FplateInt/S**2

	
C	--- Calculate the attractive and repulsive parts of the interaction for all data points
	
		h	=	hmin
		P	=	0

		Do While(h.Le.hmax)
		
			P			=	P + 1
			Points(P)	=	h
		
			If(h.Ge.1.0d0)Then		
		
				FplateInt	=	0.0d0
		
				Do N = 1, NGLeX		
		
					hplate		=	0.5d0*(2.0d0-h)*GLeX(N)+0.5d0*(2.0d0+h)
					Fplate		=	DLOG(1.0d0+(DEXP(-G)*(4.0d0/(3.0d0*S**2*Rho)))*(1.0d0-0.25d0*hplate*hplate)**1.5)
					FplateInt	=	FplateInt + GLeW(N)*Fplate		
	
				Enddo
	  
				Attraction(P)	=	-Pi*R*0.5d0*(2.0d0-h)*FplateInt/S**2 +
     &								SphereFitConstants(2)*(R*DEXP(-G)/(S**4))*(2.0d0-h)**3 + 
     &								SphereFitConstants(3)*(R*(DEXP(-G)*DEXP(-G))/(S**6))*(2.0d0-h)**3
				Repulsion(P)	=	0.0d0
		
			Else
	  
				FplateInt	=	0.0d0
	  
				Do N = 1, NGLeX		
		
					hplate		=	0.5d0*(1.0d0-h)*GLeX(N)+0.5d0*(1.0d0+h)
					hplate2		=	hplate*hplate	
					B1			=	(Pi/3.0d0)*((4.0d0-hplate2)**1.5)		
					rperp1		=	DSQRT(2.0d0-hplate2-2.0d0*DSQRT(1.0d0-hplate2))
					rperp2		=	DSQRT(2.0d0-hplate2+2.0d0*DSQRT(1.0d0-hplate2))
					rperpint	=	0.0d0	  
		
					Do M = 1, NGleX
		
						rperp	=	0.5d0*(rperp2-rperp1)*GLeX(M)+0.5d0*(rperp1+rperp2)
						frperp	=	-4.0d0*rperp*DSQRT(1.0d0-0.25d0*(hplate2+rperp*rperp))*
     &								DACOS((hplate*DSQRT(hplate2+rperp*rperp))/(2.0d0*rperp*DSQRT(1.0d0-0.25d0*(hplate2+rperp*rperp))))
						rperpint=	rperpint + GLeW(M)*frperp
			
					Enddo			
		 
					B2			=	0.5d0*(rperp2-rperp1)*rperpint
					Fplate		=	-DLOG(1.0d0+(1.0d0/(2.0d0*Pi*Rho*S**2*hplate2))*DEXP(-G)*(B1+B2))
					FplateInt	=	FplateInt + GLeW(N)*Fplate		
	
				Enddo
	  
				Attraction(P)	=	(Pi*R*0.5d0*(1.0d0-h)*FplateInt)/S**2 + Attractionh1 +
     &								SphereFitConstants(2)*(R*DEXP(-G)/(S**4))*(2.0d0-h)**3 + 
     &								SphereFitConstants(3)*(R*(DEXP(-G)*DEXP(-G))/(S**6))*(2.0d0-h)**3
				Repulsion(P)	=	(2.0d0*Pi*R/S**2)*(h*DLOG(h)-h+1.0d0) + (SphereFitConstants(1)*(1.0d0-h)**4)/(R*S**2)
	
			Endif
	
			h	=	h + dh
		
		Enddo


C	--- Calculate the total interaction and write to file
  
		Write (61,*) '# Column:'
		Write (61,*) '# 1 - Sphere-sphere separation [Ldna]'
		Write (61,*) '# 2 - Repulsive contribution [kT]'
		Write (61,*) '# 3 - Attractive contribution [kT]'
		Write (61,*) '# 4 - Total DNA-mediated interaction [kT]'
		Write (61,*) '############################################'
	  
		Do W = 1, P
			Total(W)	=	Attraction(W) + Repulsion(W)
			Write (61,9900)	Points(W), Repulsion(W), Attraction(W), Total(W)
		Enddo	
			
	Endif

	Close (61)

	Write(*,*) 'Calculation has terminated'
	Write(*,*) '*********************************************************'
	
	Endif
	
 9900    Format (ES20.10,3(',',ES20.10))	 
            
	End

!projet_scientifique.f95
!OBLIN etienne, THIVENT Raphaël , EP3
!15/05/2026
!Equation de Navier-Stokes incompressible : Application a la cavite entrainee a Re = 100

!module contenant les variables globales
module global
	implicit none 

	integer, parameter :: N=51
	integer, parameter :: nb_ite_max = 10000
	INTEGER, PARAMETER :: DBL=SELECTED_REAL_KIND(10,100)
	REAL (KIND=DBL) ::  L,dx,dy,r_tol,dt,t,t_max,BCE,BCL,CFL_restrictif,Fo_restrictif,Re,pho,nu,vol_tot,tol_integrale_vol_P,&
	integraleP,critere_arret_jacobi
	integer::counter,frequence_ecriture,type_discretisation,counter_ite_max
	REAL (KIND=DBL), dimension(N) :: x,y
	REAL (KIND=DBL), dimension(N*N) :: b,first_guess,p_1D,tab_volume
	REAL (KIND=DBL), dimension(N*N,N*N) :: A
	REAL (KIND=DBL),dimension(N,N)::u_pred,v_pred,u_init,v_init,u,v,u_np1,v_np1,p
	
end module global
	

!programme principal
program etape6
	use global
	implicit none
	 !counter_ite_max = 0
	call parametreprob()
	print *, 'L = ', L
	print *, 'N = ', N
	print *, 'Re = ', Re
	print *, 'CFL = ', CFL_restrictif
	print *, 'Fo = ', Fo_restrictif
	print *, 'r_tol = ', r_tol
	print *, 'frequence_ecriture = ', frequence_ecriture
	
	call defmesh()
	call def_volume()
	
	print *, 'dx = ', dx
	print *, 'dy = ', dy
	print*,'Volume attendu : V = ',L*L
	print*,'Volume obtenu : Vtot = ',vol_tot
	call condinit()
	
	if (type_discretisation == 1) then
		print *, 'type discretisation prediction vitesse : upwind'
		call avancement_temp1()
	elseif (type_discretisation == 2) then
		print *, 'type discretisation prediction vitesse : ordre2'	
		call avancement_temp2()
	elseif (type_discretisation == 4) then
		print *, 'type discretisation prediction vitesse : ordre4'
		call avancement_temp4()
	end if
	!print*,'Le nombre d iteration max dans le calcul de la jacobi a ete atteint ',counter_ite_max,' fois'
	call ecriture()
	
	print *,'Temps final : ', t
	print*,'-------FIN PROGRAMME-------'
endprogram etape6

!subroutine de lecture des parametres du probleme dans un fichier 
subroutine parametreprob()
    	use global
    	implicit none
    	character(len=100) :: line

    	open(10,file="input.dat")

    	read(10,'(A)') line ; read(line,*) L !lit la ligne entiere et garde juste le premier nombre trouvé et ne garde pas le reste
    	read(10,'(A)') line ; read(line,*) r_tol
    	read(10,'(A)') line ; read(line,*) frequence_ecriture
    	read(10,'(A)') line ; read(line,*) Re
    	read(10,'(A)') line ; read(line,*) pho
    	read(10,'(A)') line ; read(line,*) CFL_restrictif
    	read(10,'(A)') line ; read(line,*) Fo_restrictif
	read(10,'(A)') line ; read(line,*) t_max
    	read(10,'(A)') line ; read(line,*) BCE
    	read(10,'(A)') line ; read(line,*) BCL
    	read(10,'(A)') line ; read(line,*) tol_integrale_vol_P
    	read(10,'(A)') line ; read(line,*) type_discretisation
    	close(10)
    	
    	nu = 1.0_DBL/Re
end subroutine parametreprob

!subroutine qui definie le maillage
subroutine defmesh()
	use global
	implicit none
	integer ::  i
	dx = L / real(N-1, DBL)
	dy = L / real(N-1, DBL)
	x(1)=0.0_DBL
	y(1)=0.0_DBL
	do i=2,N
    		x(i)=x(i-1) + dx
    		y(i)=y(i-1) + dy
    	enddo
end subroutine defmesh

subroutine def_volume() !Creer un tableau de volume initial et donne une valeur a chaque case
    use global
    implicit none

    integer :: i, j,k

vol_tot = 0.0_DBL
do i = 1, N
	do j = 1, N
	k = (j-1)*N + i
		if ((i==1 .or. i==N) .and. (j==1 .or. j==N)) then
			tab_volume(k)=0.25_DBL*dx*dy
			
		elseif (i==1 .or. i==N .or. j==1 .or. j==N) then
			tab_volume(k) = 0.5_DBL*dx*dy
			
		else
			tab_volume(k)=dx*dy

		end if
		
		vol_tot = vol_tot + tab_volume(k)
	end do
end do
print*,'Different valeur de volume : dans un coin : ',&
tab_volume(1+(j-2)*N),&
'sur un cote : ',tab_volume(1+10*N),&
'Ni un coin ni un cote', tab_volume(10+15*N)
end subroutine def_volume

!fonction qui créer un index pour simplifier la lecture
function idx(i, j)
    	use global 
    	implicit none 
    	integer :: i, j
    	integer :: idx
    	idx = N*(j-1) + i
end function idx


!subroutine qui definie la condition initiale
subroutine condinit()
	use global
	implicit none
	u_init(:,:) = 0.0_DBL
	v_init(:,:) = 0.0_DBL
	u_init(:,1)=BCE !def composante u_pred
	u_init(1,:)=BCE 
	u_init(N,:)=BCE
	u_init(:,N)=BCL
        v_init(:,1)=BCE !def composante v_pred 
        v_init(:,N)=BCE
        v_init(1,:)=BCE
        v_init(N,:)=BCE
        
        u = u_init
	v = v_init
	p_1D = 0.0_DBL
endsubroutine condinit

!subroutine qui fait avancer temporellement
subroutine avancement_temp4()
	use global
	implicit none
	t=0.0_DBL
	counter=0
	call initialisation_systeme()
	Do while (t<t_max)
		call pas_de_temps()
		t=t+dt
		counter=counter+1
		u_pred = u
		v_pred = v
		!print*,'colonne j=1 u_pred : ',u_pred(:,1)
		call prediction4()
		call terme_de_droite()
		call resolution_jacobi()
		call correction()
		
		
		print *, counter 
		if ((int(MOD(counter,frequence_ecriture))==0))then     ! envoie à ecriture toute les 10 ite
       			print *, counter
    			call ecriture()
    		endif
		
	enddo 
endsubroutine avancement_temp4

!Avancement temporel avec discretisation avancement vitesse ordre 2
subroutine avancement_temp2()
	use global
	implicit none
	t=0.0_DBL
	counter=0
	call initialisation_systeme()
	Do while (t<t_max)
		call pas_de_temps()
		t=t+dt
		counter=counter+1
		u_pred = u
		v_pred = v
		!print*,'colonne j=1 u_pred : ',u_pred(:,1)
		call prediction2()
		call terme_de_droite()
		call resolution_jacobi()
		call correction()
		
		
		print *, counter 
		if ((int(MOD(counter,frequence_ecriture))==0))then     ! envoie à ecriture toute les 10 ite
       			print *, counter
    			call ecriture()
    		endif
		
	enddo 
endsubroutine avancement_temp2

!Avancement temporel avec discretisation avancement vitesse upwind
subroutine avancement_temp1()
	use global
	implicit none
	t=0.0_DBL
	counter=0
	call initialisation_systeme()
	Do while (t<t_max)
		call pas_de_temps()
		t=t+dt
		counter=counter+1
		u_pred = u
		v_pred = v
		!print*,v_pred
		!print*,'colonne j=1 u_pred : ',u_pred(:,1)
		call prediction1()
		call terme_de_droite()
		call resolution_jacobi()
		call correction()
		
		
		print *, counter 
		if ((int(MOD(counter,frequence_ecriture))==0))then     ! envoie à ecriture toute les 10 ite
       			print *, counter
    			call ecriture()
    		endif
		
	enddo 
endsubroutine avancement_temp1

! trouver le pas de temps le plus restrictif
subroutine pas_de_temps()
	use global
	implicit none
	real (KIND=DBL) ::vmax,umax, dt_CFL_u, dt_CFL_v,dt_Fo
	umax=maxval(u_pred)
    	vmax=maxval(v_pred)
	dt_Fo=(Fo_restrictif*dx**2)/nu
    	dt_CFL_u=(CFL_restrictif*dx)/umax
    	dt_CFL_v=(CFL_restrictif*dy)/vmax
	dt=min(dt_CFL_u,dt_CFL_v,dt_Fo)
	!print*,'pas de temps Fo',dt_Fo
	!print*,'pas de temps CFL_u',dt_CFL_u
	!print*,'pas de temps CFL_v',dt_CFL_v
end subroutine pas_de_temps


!subroutine pour trouver u_pred* et v_pred*
subroutine prediction4()
	use global 
	implicit none 
	real (KIND=DBL), dimension(N,N) :: u_pred_np1, v_pred_np1!u_pred apres un instant t
	real (KIND=DBL) ::Fo_effectif,CFL_v_effectif,CFL_u_effectif
	integer :: i,j

	u_pred_np1(1,:)=BCE
	u_pred_np1(N,:)=BCE
	u_pred_np1(:,1)=BCE
	u_pred_np1(:,N)=BCL

	v_pred_np1(1,:)=BCE
	v_pred_np1(N,:)=BCE
	v_pred_np1(:,1)=BCE
	v_pred_np1(:,N)=BCE
	
	Fo_effectif=(nu*dt)/(dx**2)
	!print*,'Fo effectif= ',Fo_effectif
	
	j=2 !On degrade le schema (CD 2eme ordre) pour pouvoir calculer les points en (:,j=2)
	Do i=2,N-1
		CFL_u_effectif=(u_pred(i,j)*dt)/(dx)
		CFL_v_effectif=(v_pred(i,j)*dt)/(dx)
		
		u_pred_np1(i,j) = u_pred(i,j)&
		+Fo_effectif * ( u_pred(i+1,j) - 2.0_DBL*u_pred(i,j) + u_pred(i-1,j)&
		+ u_pred(i,j+1) - 2.0_DBL*u_pred(i,j) + u_pred(i,j-1) )&
		-CFL_u_effectif/2.0_DBL *(u_pred(i+1,j)-u_pred(i-1,j))&
		-CFL_v_effectif/2.0_DBL *(u_pred(i,j+1)-u_pred(i,j-1))
		
		v_pred_np1(i,j) = v_pred(i,j)&
		+ Fo_effectif * ( v_pred(i+1,j) - 2.0_DBL*v_pred(i,j) + v_pred(i-1,j)&
		+ v_pred(i,j+1) - 2.0_DBL*v_pred(i,j) + v_pred(i,j-1) )&
		-CFL_u_effectif/2.0_DBL *(v_pred(i+1,j)-v_pred(i-1,j))&
		-CFL_v_effectif/2.0_DBL *(v_pred(i,j+1)-v_pred(i,j-1))
	end DO
	
	j=N-1 !On degrade le schema (CD 2eme ordre) pour pouvoir calculer les points en (:,j=N-1)
	Do i=2,N-1
		CFL_u_effectif=(u_pred(i,j)*dt)/(dx)
		CFL_v_effectif=(v_pred(i,j)*dt)/(dx)
		
		u_pred_np1(i,j) = u_pred(i,j)&
		+Fo_effectif * ( u_pred(i+1,j) - 2.0_DBL*u_pred(i,j) + u_pred(i-1,j)&
		+ u_pred(i,j+1) - 2.0_DBL*u_pred(i,j) + u_pred(i,j-1) )&
		-CFL_u_effectif/2.0_DBL *(u_pred(i+1,j)-u_pred(i-1,j))&
		-CFL_v_effectif/2.0_DBL *(u_pred(i,j+1)-u_pred(i,j-1))
		
		v_pred_np1(i,j) = v_pred(i,j)&
		+ Fo_effectif * ( v_pred(i+1,j) - 2.0_DBL*v_pred(i,j) + v_pred(i-1,j)&
		+ v_pred(i,j+1) - 2.0_DBL*v_pred(i,j) + v_pred(i,j-1) )&
		-CFL_u_effectif/2.0_DBL *(v_pred(i+1,j)-v_pred(i-1,j))&
		-CFL_v_effectif/2.0_DBL *(v_pred(i,j+1)-v_pred(i,j-1))
	end DO
	
	i=2 !On degrade le schema (CD 2eme ordre) pour pouvoir calculer les points en (i=2,:)
	Do j=2,N-1
		CFL_u_effectif=(u_pred(i,j)*dt)/(dx)
		CFL_v_effectif=(v_pred(i,j)*dt)/(dx)
		
		u_pred_np1(i,j) = u_pred(i,j)&
		+Fo_effectif * ( u_pred(i+1,j) - 2.0_DBL*u_pred(i,j) + u_pred(i-1,j)&
		+ u_pred(i,j+1) - 2.0_DBL*u_pred(i,j) + u_pred(i,j-1) )&
		-CFL_u_effectif/2.0_DBL *(u_pred(i+1,j)-u_pred(i-1,j))&
		-CFL_v_effectif/2.0_DBL *(u_pred(i,j+1)-u_pred(i,j-1))
		
		v_pred_np1(i,j) = v_pred(i,j)&
		+ Fo_effectif * ( v_pred(i+1,j) - 2.0_DBL*v_pred(i,j) + v_pred(i-1,j)&
		+ v_pred(i,j+1) - 2.0_DBL*v_pred(i,j) + v_pred(i,j-1) )&
		-CFL_u_effectif/2.0_DBL *(v_pred(i+1,j)-v_pred(i-1,j))&
		-CFL_v_effectif/2.0_DBL *(v_pred(i,j+1)-v_pred(i,j-1))
	end DO
	
	i=N-1 !On degrade le schema (CD 2eme ordre) pour pouvoir calculer les points en (i=N-1,:)
	Do j=2,N-1
		CFL_u_effectif=(u_pred(i,j)*dt)/(dx)
		CFL_v_effectif=(v_pred(i,j)*dt)/(dx)
		
		u_pred_np1(i,j) = u_pred(i,j)&
		+Fo_effectif * ( u_pred(i+1,j) - 2.0_DBL*u_pred(i,j) + u_pred(i-1,j)&
		+ u_pred(i,j+1) - 2.0_DBL*u_pred(i,j) + u_pred(i,j-1) )&
		-CFL_u_effectif/2.0_DBL *(u_pred(i+1,j)-u_pred(i-1,j))&
		-CFL_v_effectif/2.0_DBL *(u_pred(i,j+1)-u_pred(i,j-1))
		
		v_pred_np1(i,j) = v_pred(i,j)&
		+ Fo_effectif * ( v_pred(i+1,j) - 2.0_DBL*v_pred(i,j) + v_pred(i-1,j)&
		+ v_pred(i,j+1) - 2.0_DBL*v_pred(i,j) + v_pred(i,j-1) )&
		-CFL_u_effectif/2.0_DBL *(v_pred(i+1,j)-v_pred(i-1,j))&
		-CFL_v_effectif/2.0_DBL *(v_pred(i,j+1)-v_pred(i,j-1))
	end DO

   	Do i=3,N-2
   		Do j=3, N-2
	   		
		   		CFL_u_effectif=(u_pred(i,j)*dt)/(dx)
		   		!print*,'Le CFL effectif de u vaut : ',CFL_u_effectif
		   		CFL_v_effectif=(v_pred(i,j)*dt)/(dx)
				u_pred_np1(i,j) = u_pred(i,j)                                      &
	     + Fo_effectif * ( u_pred(i+1,j) - 2.0_DBL*u_pred(i,j) + u_pred(i-1,j)  &
		            + u_pred(i,j+1) - 2.0_DBL*u_pred(i,j) + u_pred(i,j-1) ) &
	     - CFL_u_effectif/12.0_DBL *&
	      (-u_pred(i+2,j) + 8.0_DBL* u_pred(i+1,j)-8.0_DBL* u_pred(i-1,j)+ u_pred(i-2,j))&
	      - CFL_v_effectif/12.0_DBL *&
	       (-u_pred(i,j+2) + 8.0_DBL* u_pred(i,j+1)-8.0_DBL* u_pred(i,j-1)+ u_pred(i,j-2))

		v_pred_np1(i,j) = v_pred(i,j)                                      &
	     + Fo_effectif * ( v_pred(i+1,j) - 2.0_DBL*v_pred(i,j) + v_pred(i-1,j)  &
		            + v_pred(i,j+1) - 2.0_DBL*v_pred(i,j) + v_pred(i,j-1) ) &
	     - CFL_u_effectif/12.0_DBL *&
	      (-v_pred(i+2,j) + 8.0_DBL* v_pred(i+1,j)-8.0_DBL* v_pred(i-1,j)+v_pred(i-2,j))&
	      - CFL_v_effectif/12.0_DBL *&
	      (-v_pred(i,j+2) + 8.0_DBL* v_pred(i,j+1)-8.0_DBL* v_pred(i,j-1)+v_pred(i,j-2))
        	enddo
        enddo 
    	u_pred=u_pred_np1
    	v_pred=v_pred_np1
end subroutine prediction4

!prediction avec upwind pour avancement vitesse
subroutine prediction1()
	use global 
	implicit none 
	real (KIND=DBL), dimension(N,N) :: u_pred_np1, v_pred_np1!u_pred apres un instant t
	real (KIND=DBL) ::Fo_effectif,CFL_v_effectif,CFL_u_effectif
	integer :: i,j

	u_pred_np1(1,:)=BCE
	u_pred_np1(N,:)=BCE
	u_pred_np1(:,1)=BCE
	u_pred_np1(:,N)=BCL

	v_pred_np1(1,:)=BCE
	v_pred_np1(N,:)=BCE
	v_pred_np1(:,1)=BCE
	v_pred_np1(:,N)=BCE
	
	Fo_effectif=(nu*dt)/(dx**2)
   	Do i=2,N-1
   		Do j=2, N-1
   			if (u_pred(i,j)>=0) then
	   			CFL_u_effectif=(u_pred(i,j)*dt)/(dx)
				CFL_v_effectif=(v_pred(i,j)*dt)/(dx)
	   		
						u_pred_np1(i,j)=u_pred(i,j)&
			     + Fo_effectif * ( u_pred(i+1,j) - 2.0_DBL*u_pred(i,j) + u_pred(i-1,j)  &
					    + u_pred(i,j+1) - 2.0_DBL*u_pred(i,j) + u_pred(i,j-1) ) &
			     - CFL_u_effectif * (u_pred(i,j) - u_pred(i-1,j))&
			     - CFL_v_effectif * (u_pred(i,j) - u_pred(i,j-1))
			     
			else
				CFL_u_effectif=(-u_pred(i,j)*dt)/(dx)
				CFL_v_effectif=(-v_pred(i,j)*dt)/(dx)
	   		
						u_pred_np1(i,j) = u_pred(i,j)&
			     + Fo_effectif * ( u_pred(i+1,j) - 2.0_DBL*u_pred(i,j) + u_pred(i-1,j)  &
					    + u_pred(i,j+1) - 2.0_DBL*u_pred(i,j) + u_pred(i,j-1) ) &
			     - CFL_u_effectif * (u_pred(i+1,j) - u_pred(i,j))&
			     - CFL_v_effectif * (u_pred(i,j+1) - u_pred(i,j))
			     
			end if
			
			if (v_pred(i,j)>=0) then
			
				CFL_u_effectif=(u_pred(i,j)*dt)/(dx)
				CFL_v_effectif=(v_pred(i,j)*dt)/(dx)
				
						v_pred_np1(i,j) = v_pred(i,j)&
			     + Fo_effectif * ( v_pred(i+1,j) - 2.0_DBL*v_pred(i,j) + v_pred(i-1,j)  &
					    + v_pred(i,j+1) - 2.0_DBL*v_pred(i,j) + v_pred(i,j-1) ) &
			     - CFL_u_effectif * (v_pred(i,j) - v_pred(i-1,j))&
			     - CFL_v_effectif * (v_pred(i,j) - v_pred(i,j-1))
			    
			else
				
				CFL_u_effectif=(-u_pred(i,j)*dt)/(dx)
				CFL_v_effectif=(-v_pred(i,j)*dt)/(dx)
				
						v_pred_np1(i,j) = v_pred(i,j)&
			     + Fo_effectif * ( v_pred(i+1,j) - 2.0_DBL*v_pred(i,j) + v_pred(i-1,j)  &
					    + v_pred(i,j+1) - 2.0_DBL*v_pred(i,j) + v_pred(i,j-1) ) &
			     - CFL_u_effectif * (v_pred(i+1,j) - v_pred(i,j))&
			     - CFL_v_effectif * (v_pred(i,j+1) - v_pred(i,j))
		     
		     end if
		     
        	enddo
        enddo 
    	u_pred=u_pred_np1
    	v_pred=v_pred_np1
    	!print*,v_pred
end subroutine prediction1

subroutine prediction2()
	use global 
	implicit none 
	real (KIND=DBL), dimension(N,N) :: u_pred_np1, v_pred_np1!u_pred apres un instant t
	real (KIND=DBL) ::Fo_effectif,CFL_v_effectif,CFL_u_effectif
	integer :: i,j

	u_pred_np1(1,:)=BCE
	u_pred_np1(N,:)=BCE
	u_pred_np1(:,1)=BCE
	u_pred_np1(:,N)=BCL

	v_pred_np1(1,:)=BCE
	v_pred_np1(N,:)=BCE
	v_pred_np1(:,1)=BCE
	v_pred_np1(:,N)=BCE
	
	Fo_effectif=(nu*dt)/(dx**2)
   	Do i=2,N-1
   		Do j=2, N-1
	   		CFL_u_effectif=(u_pred(i,j)*dt)/(dx)
	   		CFL_v_effectif=(v_pred(i,j)*dt)/(dx)
			u_pred_np1(i,j) = u_pred(i,j)                                      &
     + Fo_effectif * ( u_pred(i+1,j) - 2.0_DBL*u_pred(i,j) + u_pred(i-1,j)  &
                    + u_pred(i,j+1) - 2.0_DBL*u_pred(i,j) + u_pred(i,j-1) ) &
     - CFL_u_effectif/2.0_DBL * (u_pred(i+1,j) - u_pred(i-1,j))               &
     - CFL_v_effectif/2.0_DBL * (u_pred(i,j+1) - u_pred(i,j-1))

v_pred_np1(i,j) = v_pred(i,j)                                      &
     + Fo_effectif * ( v_pred(i+1,j) - 2.0_DBL*v_pred(i,j) + v_pred(i-1,j)  &
                    + v_pred(i,j+1) - 2.0_DBL*v_pred(i,j) + v_pred(i,j-1) ) &
     - CFL_u_effectif/2.0_DBL * (v_pred(i+1,j) - v_pred(i-1,j))               &
     - CFL_v_effectif/2.0_DBL * (v_pred(i,j+1) - v_pred(i,j-1))
        	enddo
        enddo 
    	u_pred=u_pred_np1
    	v_pred=v_pred_np1
end subroutine prediction2

! calcul du terme de droite de l'equation 2D de Poisson
subroutine terme_de_droite()
	use global
	implicit none
	integer :: i,j,k
	integer, external :: idx
	b(:)=0.0_DBL
	Do j=1, N
		Do i=1, N
		k=(j-1)*N + i
			if ((j==1).OR.(j==N).OR.(i==1).OR.(i==N)) then   
				b(k)=0.0_DBL
			else 
				b(k) = pho/dt * ((u_pred(i+1,j)-u_pred(i-1,j))/(2*dx) +(v_pred(i,j+1)-v_pred(i,j-1))/(2*dy))
			endif
		enddo
	enddo
end subroutine terme_de_droite

! initialisation du A pour Jacobi
subroutine initialisation_systeme()
	use global
	implicit none
	REAL (KIND=DBL) :: C1,C2
	integer :: i,j,k
	C1=1.0_DBL/dx**2
	C2=1.0_DBL/dy**2
	A(:,:)=0.0_DBL
	Do i=1, N
		Do j=1, N
		k = (j-1)*N + i
			if ((j==1).or.(j==N).or.(i==1).or.(i==N)) then
				A(k,k)=1
				if (i==1) then				! Condition de Neuman
					A(k,k+1)=-1
				else if (i==N) then
					A(k,k-1)=-1
				else if (j==1) then
					A(k,k+N)=-1
				else if (j==N) then
					A(k,k-N)=-1
				endif
			else 
				A(k,k-1)= C1
				A(k,k)= -2.0_DBL*(C1+C2)
				A(k,k+1)=C1
				A(k,k-N)=C2
				A(k,k+N)=C2
			endif
		enddo
	enddo
end subroutine initialisation_systeme


subroutine resolution_jacobi()
	use global
	implicit none
	REAL (KIND=DBL) :: somme
	integer :: i,j,nb_ite,k
	real(KIND=DBL) :: r
	REAL (KIND=DBL), dimension(N*N) :: p_np1,residus
    	nb_ite=0
    	Do 
        	Do j=1,N
            		somme = 0.0_DBL
            		Do i=1, N 
            			k = (j-1)*N + i
            			if (i == 1) then          ! Paroi Gauche 
		    			p_np1(k) = p_1D(k+1)
				elseif (i == N) then ! Paroi Droite 
					p_np1(k) = p_1D(k-1)
				elseif (j == 1) then      ! Paroi Bas 
					p_np1(k) = p_1D(k+N)
				elseif (j == N) then ! Paroi Haut
					p_np1(k) = p_1D(k-N)
				else                      ! dedans
					somme = A(k, k-1)*p_1D(k-1) + A(k, k+1)*p_1D(k+1) + &
					    A(k, k-N)*p_1D(k-N) + A(k, k+N)*p_1D(k+N)
					p_np1(k) = (1.0/A(k,k)) * (b(k) - somme)
				endif

			enddo 
		enddo  
		nb_ite= nb_ite+1
		
		critere_arret_jacobi = (sqrt((1/vol_tot)*sum(tab_volume*(p_np1-p_1D)**2)))/&
		sqrt((1/vol_tot)*sum(tab_volume*p_np1**2))
		
		
		p_1D= p_np1  

		
		integraleP = sum(p_1D-(sum(p_1D))/(N*N))
		
		
		!print *,'Ce que l on compare a la tolerance vaut : ',critere_arret_jacobi!TEST
       		!print*,'L integrale volumique de la pression vaut:',integraleP
		
		
		!if (nb_ite == nb_ite_max) then
			!counter_ite_max = counter_ite_max + 1
			
		!end if
		
       		if (((critere_arret_jacobi < r_tol).AND.(integraleP<tol_integrale_vol_P)).OR.&
       		(nb_ite>nb_ite_max))exit

	Do i=1, N
		Do j=1, N
			k = (j-1)*N + i
			p(i,j)=p_1D(k)
		enddo
	enddo
       	
	enddo
endsubroutine resolution_jacobi

subroutine correction()
	use global 
	implicit none
	integer :: i,j 
	Do i=2, N-1
		Do j=2, N-1
			u_np1(i,j) = u_pred(i,j) - dt/pho * (p(i+1,j)-p(i-1,j))/(2*dx)
			v_np1(i,j) = v_pred(i,j) - dt/pho * (p(i,j+1)-p(i,j-1))/(dy*2)
		enddo
	enddo
	u_np1(:,1)=BCE !conditions limites u
	u_np1(1,:)=BCE 
	u_np1(N,:)=BCE
	u_np1(:,N)=BCL
        v_np1(:,1)=BCE !conditions limites v
        v_np1(:,N)=BCE
        v_np1(1,:)=BCE
        v_np1(N,:)=BCE
	
	u=u_np1
	v=v_np1
endsubroutine correction			

!subroutine d'ecriture des resultats
subroutine ecriture()
	use global 
	implicit none
	integer :: i, j

	character(len=40)::nom_fichier

	write(nom_fichier, '(A,I0,A)') 'etape7_ord4_',counter,'_tecplot.dat'
	open(11,file=nom_fichier)
	write(11,*)'TITLE = "ETAPE7"'
	write(11,*)'VARIABLES="X","Y","u","v","p"'
	write(11,*)'ZONE I=',N,'J=',N,'DATAPACKING=POINT'
	do i=1,N
		do j=1, N
			write(11,*) x(i),y(j),u(i,j),v(i,j),p(i,j)
		enddo
   	enddo
    	close(11)
end subroutine ecriture




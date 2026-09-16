!projet_scientifique.f95
!OBLIN etienne, THIVENT Raphaël , EP3
!05/06/2026
!Instabilite de Rayleigh-Benard (IRB) : Etude de la convection naturelle
!Application Cellules de convection de Benard : etude en nombre de Rayleigh

!module contenant les variables globales
module global
	implicit none 

	integer, parameter :: Nx=81, Ny=41
	integer, parameter :: nb_ite_max = 10000
	INTEGER, PARAMETER :: DBL=SELECTED_REAL_KIND(10,100)
	REAL (KIND=DBL) :: Lx,Ly,dx,dy,r_tol,dt,t,t_max,CFL_restrictif,Fo_restrictif &
			,Re,pho,nu,vol_tot,tol_integrale_vol_P,g,Ra,T0,T1,T_ref,beta,alpha,IC_v,IC_p,BC_v,BC_p 
	integer::counter,frequence_ecriture
	REAL (KIND=DBL), dimension(Nx) :: x
	REAL (KIND=DBL), dimension(Ny) :: y
	REAL (KIND=DBL), dimension(Nx*Ny) :: b,first_guess,p_1D,tab_volume
	REAL (KIND=DBL), dimension(Nx*Ny,Nx*Ny) :: A
	REAL (KIND=DBL),dimension(Nx,Ny)::u_pred,v_pred,u,v,u_np1,v_np1,p,Temp, T_init 
	
end module global
	

!programme principal
program irb3
	use global
	implicit none
	 
	call parametreprob()
	call determination_Ra()
	print *, '-------------------PARAMETRES PROBLEME------------------------- '
	print *, 'Nx = ', Nx
	print *, 'Ny = ', Ny

	
	print *, 'Lx = ', Lx 
    	print *, 'Ly = ',Ly
    	print *, 'r_tol = ',r_tol
    	print *, 'freq ecriture = ', frequence_ecriture
    	print *, 'Re = ', Re
    	print *, 'pho = ', pho
    	print *, 'CFL = ', CFL_restrictif
	print *, 'Fo = ', Fo_restrictif
	print *, 't_max = ', t_max
	print *, 'nu = ', nu
	print *, 'alpha = ',alpha
        print *, 'IC_v = ',IC_v
	print *, 'IC_p = ', IC_p
	print *, 'BC_v = ', BC_v
	print *, 'BC_p = ', BC_p
	print *, 'T_ref = ',T_ref
	print *, 'T0 = ', T0
	print *, 'T1 = ', T1
	print *, 'beta = ', beta
	print *, 'g = ', g
	print *, 'Ra = ', Ra
	

	call defmesh()
	call def_volume()
	
	print *, 'dx = ', dx
	print *, 'dy = ', dy
	print*,'Volume attendu : V = ',Lx*Ly
	print*,'Volume obtenu : Vtot = ',vol_tot
	call condinit()
		
	call avancement_temp()
	
	
	call ecriture()
	
	print *,'Temps final : ', t
	print*,'-------FIN PROGRAMME-------'
endprogram irb3

!subroutine de lecture des parametres du probleme dans un fichier 
subroutine parametreprob()
    	use global
    	implicit none
    	character(len=100) :: line

    	open(10,file="input.dat")

    	read(10,'(A)') line ; read(line,*) Lx !lit la ligne entiere et garde juste le premier nombre trouvé et ne garde pas le reste
    	read(10,'(A)') line ; read(line,*) Ly
    	read(10,'(A)') line ; read(line,*) r_tol
    	read(10,'(A)') line ; read(line,*) frequence_ecriture
    	read(10,'(A)') line ; read(line,*) Re
    	read(10,'(A)') line ; read(line,*) pho
    	read(10,'(A)') line ; read(line,*) CFL_restrictif
    	read(10,'(A)') line ; read(line,*) Fo_restrictif
	read(10,'(A)') line ; read(line,*) t_max
	read(10,'(A)') line ; read(line,*) nu
	read(10,'(A)') line ; read(line,*) alpha
        read(10,'(A)') line ; read(line,*) IC_v
	read(10,'(A)') line ; read(line,*) IC_p
	read(10,'(A)') line ; read(line,*) BC_v
	read(10,'(A)') line ; read(line,*) BC_p
	read(10,'(A)') line ; read(line,*) T_ref
	read(10,'(A)') line ; read(line,*) T0
	read(10,'(A)') line ; read(line,*) T1
	read(10,'(A)') line ; read(line,*) beta
	read(10,'(A)') line ; read(line,*) g
    	close(10)
    	

end subroutine parametreprob

!subroutine qui definie le maillage
subroutine defmesh()
	use global
	implicit none
	integer ::  i
	dx = Lx / real(Nx-1, DBL)
	dy = Ly / real(Ny-1, DBL)
	x(1)=0.0_DBL
	y(1)=0.0_DBL
	do i=2,Nx
    		x(i)=x(i-1) + dx
    	enddo
    	do i=2,Ny
    		y(i)=y(i-1) + dy
    	enddo
end subroutine defmesh

subroutine def_volume() !Creer un tableau de volume initial et donne une valeur a chaque case
	use global
        implicit none

        integer :: i, j,k

	vol_tot = 0.0_DBL
	do i = 1, Nx
		do j = 1, Ny
			k = (j-1)*Nx + i
			if ((i==1 .or. i==Nx) .and. (j==1 .or. j==Ny)) then
				tab_volume(k)=0.25_DBL*dx*dy
				
			elseif (i==1 .or. i==Nx .or. j==1 .or. j==Ny) then
				tab_volume(k) = 0.5_DBL*dx*dy
				
			else
				tab_volume(k)=dx*dy

			end if
			
			vol_tot = vol_tot + tab_volume(k)
		end do
	end do
	print*,'Different valeur de volume : dans un coin : ',&
	tab_volume(1+(j-2)*Nx),&
	'sur un cote : ',tab_volume(1+10*Nx),&
	'Ni un coin ni un cote', tab_volume(10+15*Nx)
end subroutine def_volume

!fonction qui créer un index pour simplifier la lecture
function idx(i, j)
    	use global 
    	implicit none 
    	integer :: i, j
    	integer :: idx
    	idx = Nx*(j-1) + i
end function idx


!subroutine qui definie la condition initiale
subroutine condinit()
	use global
	implicit none
        
        u = IC_v
	v = IC_v
	p_1D = IC_p
	call init_temperature()
	Temp = T_init
endsubroutine condinit

!subroutine qui fait avancer temporellement
subroutine avancement_temp()
	use global
	implicit none
	t=0.0_DBL
	counter=0
	call bc_periodique_x()
	call ecriture()
	call initialisation_systeme()
	Do while (t<t_max)
		call pas_de_temps()
		t=t+dt
		counter=counter+1
		u_pred = u
		v_pred = v
		call temperature()
		call bc_periodique_x()
		call prediction()
		call bc_periodique_x()
		call terme_de_droite()
		call resolution_jacobi()
		call correction()
		call bc_periodique_x()
		
		
		print *, counter 
		if ((int(MODULO(counter,frequence_ecriture))==0))then     ! envoie à ecriture toute les 10 ite
       			print *, counter
    			call ecriture()
    		endif
		
	enddo 
endsubroutine avancement_temp

! trouver le pas de temps le plus restrictif
subroutine pas_de_temps()
	use global
	implicit none
	real (KIND=DBL) ::vmax,umax, dt_CFL_u, dt_CFL_v,dt_Fo_V,dt_Fo_T
	umax=maxval(u_pred)
    	vmax=maxval(v_pred)
    	umax = max( maxval(abs(u_pred)), 1e-12_DBL )  ! pour ne pas div par 0
	vmax = max( maxval(abs(v_pred)), 1e-12_DBL )
    	dt_CFL_u=(CFL_restrictif*dx)/umax
    	dt_CFL_v=(CFL_restrictif*dy)/vmax
	dt_Fo_T = Fo_restrictif * min(dx,dy)**2 / alpha
	dt_Fo_V = Fo_restrictif * min(dx,dy)**2 / nu

	dt = min(dt_CFL_u, dt_CFL_v, dt_Fo_T, dt_Fo_V)
end subroutine pas_de_temps

! subroutine qui assure la periodicite en x
subroutine bc_periodique_x()
        use global
        implicit none

        integer :: j

        u(1,:)  = u(Nx-1,:)
        u(Nx,:) = u(2,:)

        v(1,:)  = v(Nx-1,:)
        v(Nx,:) = v(2,:)

        p(1,:)  = p(Nx-1,:)
        p(Nx,:) = p(2,:)

        u_pred(1,:)  = u_pred(Nx-1,:)
        u_pred(Nx,:) = u_pred(2,:)

        v_pred(1,:)  = v_pred(Nx-1,:)
	v_pred(Nx,:) = v_pred(2,:)
	
	Temp(1,:)  = Temp(Nx-1,:)
	Temp(Nx,:) = Temp(2,:)

end subroutine bc_periodique_x


subroutine temperature()
	use global 
	implicit none 
	real(KIND=DBL), dimension(Nx,Ny) :: T_np1
	real(KIND=DBL) :: Fo_effectif_x, Fo_effectif_y,adv_u, adv_v
	integer :: i, j, im1, ip1
 
	! Conditions aux limites en y (Dirichlet)
	T_np1(:,1)  = T0
	T_np1(:,Ny) = T1
 
	Fo_effectif_x = (alpha * dt) / (dx**2)
 	Fo_effectif_y = (alpha * dt) / (dy**2)
 	
	Do i = 1, Nx
		Do j = 2, Ny-1
 
			im1 = MODULO(i-2 + Nx, Nx) + 1   ! periodicite en x
			ip1 = MODULO(i   + Nx, Nx) + 1
 
			!Advection en x : upwind change avec le signe de u 
			if (u(i,j) >= 0.0_DBL) then
				adv_u = u(i,j) * (Temp(i,j) - Temp(im1,j)) / dx
			else
				adv_u = u(i,j) * (Temp(ip1,j) - Temp(i,j)) / dx
			end if
 
			!Advection en y : upwind change avec le signe de v
			if (v(i,j) >= 0.0_DBL) then
				adv_v = v(i,j) * (Temp(i,j) - Temp(i,j-1)) / dy
			else
				adv_v = v(i,j) * (Temp(i,j+1) - Temp(i,j)) / dy
			end if
 
			!diffusion + advection
			T_np1(i,j) = Temp(i,j) &
			+ Fo_effectif_x * ( Temp(ip1,j) - 2.0_DBL*Temp(i,j) + Temp(im1,j)) &
			                + (Temp(i,j+1) - 2.0_DBL*Temp(i,j) + Temp(i,j-1) )*Fo_effectif_y &
			- dt * adv_u &
			- dt * adv_v
 
		enddo
	enddo
	Temp = T_np1
end subroutine temperature

subroutine init_temperature()
	use global 
	implicit none 
	real (KIND=DBL)::random_val
	integer :: i,j
	T_init(:,1)  = T0
	T_init(:,Ny) = T1
	
	Do j = 2, Ny-1
		Do i = 1, Nx
			call random_number(random_val)
			! Profil lineaire de base + bruit aleatoire 2D 
			T_init(i,j) = T0 + (T1 - T0) * real(j-1, DBL) / real(Ny-1, DBL) &
			           + (2.0_DBL * random_val - 1.0_DBL) * abs(T1 - T0) / 10.0_DBL
		enddo
        enddo  
end subroutine init_temperature

subroutine determination_Ra()
	use global
	implicit none
	real (KIND=DBL) :: delta_T
	delta_T=abs(T0-T1)
	Ra = (g*Ly**3*beta*delta_T)/(nu*alpha)
endsubroutine determination_Ra

!subroutine pour trouver u_pred* et v_pred*
subroutine prediction()
	use global 
	implicit none 
	real (KIND=DBL), dimension(Nx,Ny) :: u_pred_np1, v_pred_np1!u_pred apres un instant t
	real (KIND=DBL) ::Fo_effectif_x,Fo_effectif_y,CFL_v_effectif,CFL_u_effectif
	integer :: i,j, im1,ip1,im2,ip2

	u_pred_np1(:,1)=BC_v
	u_pred_np1(:,Ny)=BC_v

	v_pred_np1(:,1)=BC_v
	v_pred_np1(:,Ny)=BC_v
	
	Fo_effectif_x=(nu*dt)/(dx**2)
	Fo_effectif_y=(nu*dt)/(dy**2)

	
	j=2 !On degrade le schema (CD 2eme ordre) pour pouvoir calculer les points en (:,j=2)
	Do i=2,Nx-1
		CFL_u_effectif=(u_pred(i,j)*dt)/(dx)
		CFL_v_effectif=(v_pred(i,j)*dt)/(dy)
		
		im1 = MODULO(i-2 + Nx, Nx) + 1   ! car quand on passe de l autre cote 
        	ip1 = MODULO(i   + Nx, Nx) + 1

        
		u_pred_np1(i,j) = u_pred(i,j)&
		+Fo_effectif_x * ( u_pred(ip1,j) - 2.0_DBL*u_pred(i,j) + u_pred(im1,j))&
		+ (u_pred(i,j+1) - 2.0_DBL*u_pred(i,j) + u_pred(i,j-1) )*Fo_effectif_y&
		-CFL_u_effectif/2.0_DBL *(u_pred(ip1,j)-u_pred(im1,j))&
		-CFL_v_effectif/2.0_DBL *(u_pred(i,j+1)-u_pred(i,j-1))
		
		v_pred_np1(i,j) = v_pred(i,j)&
		+ Fo_effectif_x * ( v_pred(ip1,j) - 2.0_DBL*v_pred(i,j) + v_pred(im1,j))&
		+ (v_pred(i,j+1) - 2.0_DBL*v_pred(i,j) + v_pred(i,j-1) )*Fo_effectif_y&
		-CFL_u_effectif/2.0_DBL *(v_pred(ip1,j)-v_pred(im1,j))&
		-CFL_v_effectif/2.0_DBL *(v_pred(i,j+1)-v_pred(i,j-1))+ dt*Beta*g*(Temp(i,j)-T_ref)
	end DO
	
	j=Ny-1 !On degrade le schema (CD 2eme ordre) pour pouvoir calculer les points en (:,j=N-1)
	Do i=2,Nx-1
		CFL_u_effectif=(u_pred(i,j)*dt)/(dx)
		CFL_v_effectif=(v_pred(i,j)*dt)/(dy)
		
		im1 = MODULO(i-2 + Nx, Nx) + 1   ! car quand on passe de l autre cote 
        	ip1 = MODULO(i   + Nx, Nx) + 1
        	
		u_pred_np1(i,j) = u_pred(i,j)&
		+Fo_effectif_x * ( u_pred(ip1,j) - 2.0_DBL*u_pred(i,j) + u_pred(im1,j))&
		+ (u_pred(i,j+1) - 2.0_DBL*u_pred(i,j) + u_pred(i,j-1) )*Fo_effectif_y&
		-CFL_u_effectif/2.0_DBL *(u_pred(ip1,j)-u_pred(im1,j))&
		-CFL_v_effectif/2.0_DBL *(u_pred(i,j+1)-u_pred(i,j-1))
		
		v_pred_np1(i,j) = v_pred(i,j)&
		+ Fo_effectif_x * ( v_pred(ip1,j) - 2.0_DBL*v_pred(i,j) + v_pred(im1,j))&
		+ (v_pred(i,j+1) - 2.0_DBL*v_pred(i,j) + v_pred(i,j-1) )*Fo_effectif_y &
		-CFL_u_effectif/2.0_DBL *(v_pred(ip1,j)-v_pred(im1,j))&
		-CFL_v_effectif/2.0_DBL *(v_pred(i,j+1)-v_pred(i,j-1))+ dt*Beta*g*(Temp(i,j)-T_ref)
	end DO
	


   	Do i=1,Nx
   		Do j=3, Ny-2
	   		im1 = MODULO(i-2 + Nx, Nx) + 1   ! car quand on passe de l autre cote 
        		ip1 = MODULO(i   + Nx, Nx) + 1
        		im2 = MODULO(i-3 + Nx, Nx) + 1   
        		ip2 = MODULO(i+1 + Nx, Nx) + 1
		   	CFL_u_effectif=(u_pred(i,j)*dt)/(dx)
		   	CFL_v_effectif=(v_pred(i,j)*dt)/(dy)
		   	
			u_pred_np1(i,j) = u_pred(i,j)                                      &
	     + Fo_effectif_x * ( u_pred(ip1,j) - 2.0_DBL*u_pred(i,j) + u_pred(im1,j))  &
		            + (u_pred(i,j+1) - 2.0_DBL*u_pred(i,j) + u_pred(i,j-1) )*Fo_effectif_y &
	     - CFL_u_effectif/12.0_DBL *&
	      (-u_pred(ip2,j) + 8.0_DBL* u_pred(ip1,j)-8.0_DBL* u_pred(im1,j)+ u_pred(im2,j))&
	      - CFL_v_effectif/12.0_DBL *&
	       (-u_pred(i,j+2) + 8.0_DBL* u_pred(i,j+1)-8.0_DBL* u_pred(i,j-1)+ u_pred(i,j-2))

			v_pred_np1(i,j) = v_pred(i,j)                                      &
	     + Fo_effectif_x * ( v_pred(ip1,j) - 2.0_DBL*v_pred(i,j) + v_pred(im1,j))  &
		            + (v_pred(i,j+1) - 2.0_DBL*v_pred(i,j) + v_pred(i,j-1) )*Fo_effectif_y &
	     - CFL_u_effectif/12.0_DBL *&
	      (-v_pred(ip2,j) + 8.0_DBL* v_pred(ip1,j)-8.0_DBL* v_pred(im1,j)+v_pred(im2,j))&
	      - CFL_v_effectif/12.0_DBL *&
	      (-v_pred(i,j+2) + 8.0_DBL* v_pred(i,j+1)-8.0_DBL* v_pred(i,j-1)+v_pred(i,j-2))+ dt*Beta*g*(Temp(i,j)-T_ref)
        	enddo
        enddo 
    	u_pred=u_pred_np1
    	v_pred=v_pred_np1
end subroutine prediction


! calcul du terme de droite de l'equation 2D de Poisson
subroutine terme_de_droite()
	use global
	implicit none

	integer :: i, j, k, im1, ip1

	b(:) = 0.0_DBL

	Do j = 1, Ny
		Do i = 1, Nx

			k = (j-1)*Nx + i

			!PERIODICITE EN X 
			im1 = MODULO(i-2 + Nx, Nx) + 1
			ip1 = MODULO(i   + Nx, Nx) + 1

			!BORDS EN Y 
			if (j == 1 .or. j == Ny) then

				! Neumann 
				b(k) = 0.0_DBL

			else

				b(k) = pho/dt * ( &
				(u_pred(ip1,j) - u_pred(im1,j)) / (2.0_DBL*dx) + &
				(v_pred(i,j+1) - v_pred(i,j-1)) / (2.0_DBL*dy) )

			endif

		end do
	end do

end subroutine terme_de_droite

! initialisation du A pour Jacobi
subroutine initialisation_systeme()
	use global
	implicit none
	
	REAL (KIND=DBL) :: C1, C2
	integer :: i, j, k, im1, ip1

	C1 = 1.0_DBL/dx**2
	C2 = 1.0_DBL/dy**2

	A(:,:) = 0.0_DBL

	Do i = 1, Nx
		Do j = 1, Ny

			k = (j-1)*Nx + i

			!PERIODICITE EN X 
			im1 = MODULO(i-2 + Nx, Nx) + 1
			ip1 = MODULO(i   + Nx, Nx) + 1

			!BORDS EN Y 
			if (j == 1 .or. j == Ny) then


				A(k,k) = 1.0_DBL
				
				! Neumann 
				if (j == 1) then
					A(k,k+Nx) = -1.0_DBL
				else
					A(k,k-Nx) = -1.0_DBL
				endif

			else

				!INTERIEUR 
				A(k,k)   = -2.0_DBL*(C1 + C2)
				A(k,im1) = C1
				A(k,ip1) = C1
				A(k,k-Nx)= C2
				A(k,k+Nx)= C2

			endif

		end do
	end do

end subroutine initialisation_systeme


subroutine resolution_jacobi()
	use global
	implicit none
	REAL (KIND=DBL) :: somme,somme_res
	integer :: i,j,nb_ite,k,kip1,kim1,ip1,im1
	real(KIND=DBL) :: integraleP,critere_arret_jacobi
	REAL (KIND=DBL), dimension(Nx*Ny) :: p_np1,residus
    	nb_ite=0
    	Do 
        	Do j=1,Ny
            		somme = 0.0_DBL
            		Do i=1, Nx 
            			k = (j-1)*Nx + i
            			im1 = MODULO(i-2 + Nx, Nx) + 1
				ip1 = MODULO(i   + Nx, Nx) + 1
				kim1 = (j-1)*Nx + im1
				kip1 = (j-1)*Nx + ip1

				if (j == 1) then      ! Paroi Bas 
					p_np1(k) = p_1D(k+Nx)
				elseif (j == Ny) then ! Paroi Haut
					p_np1(k) = p_1D(k-Nx)
				else                      ! dedans
					somme = A(k, kim1)*p_1D(kim1) + A(k, kip1)*p_1D(kip1) + &
        				A(k, k-Nx)*p_1D(k-Nx) + A(k, k+Nx)*p_1D(k+Nx)
					p_np1(k) = (1.0_DBL/A(k,k)) * (b(k) - somme)
				endif

			enddo 
		enddo  
		nb_ite= nb_ite+1
		
		critere_arret_jacobi = (sqrt((1/vol_tot)*sum(tab_volume*(p_np1-p_1D)**2)))/&
		sqrt((1/vol_tot)*sum(tab_volume*p_np1**2))
		
		
		p_1D= p_np1  

		
		integraleP = sum(p_1D-(sum(p_1D))/(Nx*Ny))

		
       		if ((critere_arret_jacobi < r_tol).OR.&
       		(nb_ite>nb_ite_max))exit
       		
       	enddo

	Do i=1, Nx
		Do j=1, Ny
			k = (j-1)*Nx + i
			p(i,j)=p_1D(k)
		enddo
	enddo
       	
	
endsubroutine resolution_jacobi

subroutine correction()
	use global 
	implicit none
	integer :: i,j,ip1,im1 
	Do i=1, Nx
		Do j=2, Ny-1
			im1 = MODULO(i-2 + Nx, Nx) + 1
			ip1 = MODULO(i   + Nx, Nx) + 1  ! car quand on passe de l autre cote 
			u_np1(i,j) = u_pred(i,j) - dt/pho * (p(ip1,j)-p(im1,j))/(2*dx)
			v_np1(i,j) = v_pred(i,j) - dt/pho * (p(i,j+1)-p(i,j-1))/(dy*2) 
		enddo
	enddo
	u_np1(:,1)=BC_v !conditions limites u
	u_np1(:,Ny)=BC_v
        v_np1(:,1)=BC_v !conditions limites v
        v_np1(:,Ny)=BC_v

	
	u=u_np1
	v=v_np1
endsubroutine correction			

!subroutine d'ecriture des resultats
subroutine ecriture()
	use global 
	implicit none
	integer :: i, j

	character(len=40)::nom_fichier

	write(nom_fichier, '(A,I0,A)') 'irb3',counter,'_tecplot.dat'
	open(11,file=nom_fichier)
	write(11,*)'TITLE = "IRB3"'
	write(11,*)'VARIABLES="X","Y","u","v","Temp","p"'
	write(11,*)'ZONE I=',Ny,'J=',Nx,'DATAPACKING=POINT'
	do i=1,Nx
		do j=1, Ny
			write(11,*) x(i),y(j),u(i,j),v(i,j), Temp(i,j), p(i,j)
		enddo
   	enddo
    	close(11)
end subroutine ecriture

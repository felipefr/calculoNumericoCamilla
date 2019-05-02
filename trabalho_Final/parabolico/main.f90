subroutine fonte(f,x,y,t,alpha,kappa)
	use precisao
	
	implicit none 
	real(wp) :: f,x,y,t,alpha,kappa
	
	f = (2.0d0*PI * kappa - alpha)*dexp(-alpha*t)*dsin(PI*x)*dsin(PI*y)

end subroutine


subroutine func_u_0(u0,x,y)
	use precisao
	
	implicit none
	
	real(wp) :: u0,x,y
	
	u0 = dsin(PI*x)*dsin(PI*y)

end subroutine



program main 
	use resolvedoresSistemaLinear
	
	implicit none
	
	integer :: Nx, Ny, N, i, j, k, Ntempos, nn
	
	! A : Matriz do sistema, F : lado direito do sistema, U : vetor solucao (temperaturas) 
	real(wp) , allocatable :: U(:), Un(:) ! U é a solucao atual
	
	! propriedades fisicas
	real(wp) :: hx, hy, Lx, Ly, kappa, alpha, Dt , sigma_x, sigma_y, xi, yj, tn, T, f
	
	integer , parameter :: saida = 14
	
	open(saida , file = 'saida.txt')
	
	! difusibilidade
	kappa = 1.0d0
	alpha = 10.0d0
	
	Nx = 21
	Ny = 21
	
	Lx = 1.0d0
	Ly = 1.0d0
	
	T = 0.01d0
	Dt = 0.001d0
	Ntempos = int(T/Dt)
	
	N = Nx*Ny
	
	
	hx = Lx/real(Nx - 1)
	hy = Ly/real(Ny - 1)
	
	
	sigma_x = kappa*Dt/hx**2.0d0
	sigma_y = kappa*Dt/hy**2.0d0
	
	allocate(U(N),Un(N))
	
	write(saida,*) Ntempos + 1
	write(saida,*) N
	
	
	! inicializa condicao inicial
	do i = 1, Nx
		xi = (i-1)*hx
		do j = 1, Ny
			yj = (j-1)*hy
			
			k = (j - 1)*Nx + i
			
			call func_u_0(Un(k),xi,yj)
			
			write(saida,*) Un(k)
			
		end do
	end do
	
	U = Un
	
!~ 	! calcula matriz sem condicao de contorno
!~ 	! i diz a coluna, varia horizontalmente da esquerda pra direita
!~ 	! j diz a linha, varia verticalmente de cima pra baixo
	do nn = 1, Ntempos + 1
		tn = (nn-1)*Dt
		
		do i = 2, Nx-1
			xi = (i-1)*hx
			do j = 2, Ny-1
				yj = (j-1)*hy
				
				k = (j - 1)*Nx + i
				
				call fonte(f,xi,yj,tn,alpha,kappa)
				
				U(k) = sigma_x*(Un(k+1) +  Un(k-1) -2.0d0*Un(k)) + sigma_y*(Un(k+Nx) + Un(k-Nx) - 2.0d0*Un(k)) + Un(k) + f*Dt
			
			end do
		end do

		do k = 1, N
			write(saida,*) U(k)
		end do
		
		Un = U
		
	end do
	
	close(saida)
	
end program




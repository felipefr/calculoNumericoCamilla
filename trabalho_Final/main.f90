program main 
	use resolvedoresSistemaLinear
	
	implicit none
	
	integer :: Nx, Ny, N, i, j, k
	
	! A : Matriz do sistema, F : lado direito do sistema, U : vetor solucao (temperaturas) 
	real(wp) , allocatable :: A(:,:), U(:), F(:)
	
	! propriedades fisicas
	real(wp) :: hx, hy, Lx, Ly, kappa, alpha_x, alpha_y, beta , theta_bar
	
	integer , parameter :: saida = 14
	
	open(saida , file = 'saida.txt')
	
	! difusibilidade
	kappa = 1.0d0
	theta_bar = 1.0d0
	
	Nx = 20
	Ny = 20
	
	Lx = 1.0d0
	Ly = 1.0d0
	
	N = Nx*Ny
	
	allocate(A(N,N),U(N),F(N))
	
	A = 0.0d0
	
	U = 0.0d0
	F = 0.0d0
	
	hx = Lx/real(Nx + 1)
	hy = Ly/real(Ny + 1)
	
	
	alpha_x = -kappa/hx**2.0d0
	alpha_y = -kappa/hy**2.0d0
	
	
	
	beta = - 2.0d0*(alpha_x + alpha_y)
	
	! calcula matriz sem condicao de contorno
	! i diz a coluna, varia horizontalmente da esquerda pra direita
	! j diz a linha, varia verticalmente de cima pra baixo
	do i = 1, Nx
		do j = 1, Ny
			k = (j - 1)*Nx + i
			A(k,k) = beta
			
			if(i>1) then
 				A(k,k-1) = alpha_x
			end if
			
			if(j>1) then
				A(k, k - Nx) = alpha_y
			end if

			if(i < Nx) then
				A(k,k + 1) = alpha_x
			end if
			if(j < Ny) then
				A(k, k + Nx) = alpha_y
			end if
			
		end do
	end do
	
	! Insere condicoes de contorno
	! face esquerda
	i = 1
	do j = 1, Ny
		k = (j - 1)*Nx + i
		A(k,:) = 0.0d0
		A(k,k) = 1.0d0
		F(k) = 0.0
	end do
	
	! face direita
	i = Nx
	do j = 1, Ny
		k = (j - 1)*Nx + i
		A(k,:) = 0.0d0
		A(k,k) = 1.0d0
		F(k) = 0.0
	end do
	
	! face superior
	j = 1
	do i = 1, Nx
		k = (j - 1)*Nx + i
		A(k,:) = 0.0d0
		A(k,k) = 1.0d0
		F(k) = theta_bar
	end do
	
	! face inferior
	j = Ny
	do i = 1, Nx
		k = (j - 1)*Nx + i
		A(k,:) = 0.0d0
		A(k,k) = 1.0d0
		F(k) = 0.0
	end do
	
	call resolvePorLU(A,F,U,N)
	
	write(saida,*) U
	
	close(saida)

	
end program




program main 
	use resolvedoresSistemaLinear
	
	implicit none
	
	integer :: Nx, Ny, N, i, j, k
	
	real(wp) , allocatable :: A(:,:), U(:), F(:)
	real(wp) :: hx, hy, Lx, Ly, kappa, alpha_x, alpha_y, beta , theta_bar
	
	kappa = 1.0d0
	theta_bar = 10.0d0
	
	Nx = 10
	Ny = 10
	
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
			
			if( j == 1) then
				F(k) = -alpha_y*theta_bar
			end if

		end do
	end do
	
	call resolvePorLU(A,F,U,N)
	
	write(0,*) U
	
	
end program




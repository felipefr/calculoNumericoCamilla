module resolvedoresSistemaLinear
	
	use precisao
	
	implicit none
	
	
	contains
	
	! resolve sistema linear por decomposição LU, implementada na rotina fatoraLU
	subroutine inverteMatriz_porLU(A,Ainv,n) 
		
		integer, intent(in) :: n
		real(wp) , intent(inout) :: A(n,n)
		real(wp) , intent(out) :: Ainv(n,n)
		
		real(wp) :: ei(n), x(n)
		integer :: ipiv(n) , i

		! fatoração de A = LU , propriamente dita
		call fatoraLU(A, ipiv)
		
		ei = 0.0d0 ! canonical vector
		do i = 1,n
			ei(i) = 1.0d0
			call substituicoes(A,ipiv,ei,x,n)	
			Ainv(:,i) = x(:)
			x = 0.0d0
			ei(i) = 0.0d0
		end do
				
	end subroutine

	subroutine resolvePorLU(A,b,x,n) 
		
		integer, intent(in) :: n
		real(wp) , intent(inout) :: b(n), A(n,n), x(n)

		integer :: ipiv(n)

		! fatoração de A = LU , propriamente dita
		call fatoraLU(A, ipiv)
				
		call substituicoes(A,ipiv,b,x,n)

	end subroutine
	
	subroutine substituicoes(LU,ipiv,b,x,n)
		
		integer, intent(in) :: n,  ipiv(n)
		real(wp) , intent(inout) :: b(n), LU(n,n), x(n)

		integer :: i
		real(wp) :: y(n)

		! substituição para frente (foward-substitution)
		do i=1,n
		  y(i)= b(ipiv(i)) - dot_product(LU(ipiv(i),1:i-1),y(1:i-1)) 
		end do
		
		! substituição para trás (backward-substitution)
		do i = n,1,-1
		  x(i)=(y(i)-dot_product(LU(ipiv(i),i+1:),x(i+1:))) / LU(ipiv(i),i)
		end do
	
	end subroutine

	! fatoração LU de A, armazena em A
	subroutine fatoraLU(a,p)

		real(wp), intent(inout) :: A(:,:)
		integer, intent(out) :: p(:)
		integer                :: n, i,j,k,kmax
		n = size(a,1)
		p = [ ( i, i=1,n ) ]
		do k = 1,n-1
			kmax = maxloc(abs(a(p(k:),k)),1) + k-1
			if (kmax /= k ) p([k, kmax]) = p([kmax, k])
			A(p(k+1:),k) = A(p(k+1:),k) / A(p(k),k)
			forall (j=k+1:n) a(p(k+1:),j) = A(p(k+1:),j) - A(p(k+1:),k) * A(p(k),j)
		end do

	end subroutine
	
	
	!Gauss seidel
	subroutine resolvePorGaussSiedel(A,b,x,n)

		real(wp), intent(inout):: A(n,n), b(n)
		integer , intent(in) :: n
		real(wp), intent(out) :: x(n)
		
		integer :: ITER,ITERMAX,i,j 
	
		real(wp):: soma, TOL, erro
		real(wp) :: xold(n)
		
		ITERMAX=1e5
		TOL =1.0e-6
	
		ITER=0
		erro = 999.9

		do while( ITER < ITERMAX .and. erro > tol) 
		  
		  xold = x ! x inicia com um chute inicial
			do i = 1,n
				soma=0.0
				do j = 1, n
				  if ( i/=j ) then
					soma=soma+A(i,j)*x(j)
				  end if
				end do
							
			x(i)=(b(i)-soma)/A(i,i)
			end do
		
		  ITER=ITER+1
		  erro = sqrt((dot_product(x - xold,x - xold)))/ sqrt((dot_product(x,x))) ! solução não pode ser zero
		  write(0,*) erro
		end do 
		
		
	end subroutine
	
		! resolve sistema linear por funções da lapack LU
	subroutine resolveSistemaLapack(A,b,x,n)
	
		integer , intent(in) :: n
		real(wp) , intent(inout) :: A(n,n) 
		real(wp) , intent(in) :: b(n) 
		real(wp) , intent(out) :: x(n) 
 		
		
		integer :: lda, ipiv(n), info
		
		lda = n

		x = b ! copia solução

		!  fatora a matriz em LU
		if(wp == sp) then 
			call sgetrf ( n, n, A, lda, ipiv, info )
		else if(wp == dp) then
			call dgetrf ( n, n, A, lda, ipiv, info )
		end if
		
		if ( info .ne. 0 ) then
		write ( *, '(a)' ) ' '
		write ( *, '(a,i8)' ) '  Matrix is singular, INFO = ', info
		return
		end if

		!  resolve sistema linear.

		if(wp == sp) then 
			call sgetrs ( 'N', n, 1, A, lda, ipiv, x, n, info )
		else if(wp == dp) then
			call dgetrs ( 'N', n, 1, A, lda, ipiv, x, n, info )
		end if
		
		

		if ( info .ne. 0 ) then
		write ( *, '(a)' ) ' '
		write ( *, '(a,i8)' ) '  Solution procedure failed, INFO = ', info
		return
		end if
	
	
	end subroutine
	
end module

subroutine func (fx,x)
implicit none
real*8 :: x,fx

fx= 2.0d0*dcos(x) - 0.5d0*dexp(x)
end subroutine

subroutine dfunc (dfx,x)
implicit none
real*8 :: dfx,x

dfx= (-2.0d0)*dsin(x) - 0.5d0*dexp(x)

end subroutine

program main 
	implicit none 
	integer:: k,kmax
	real*8:: xk1, xk,fk,dfk,tol,erro !
	! valores inicias
	tol=1.0e-5
	kmax=100
	
	k=0
	erro=999.9
	
	xk=0.0d0
	
	do while(erro>tol .and. k<kmax)
	
		call func (fk,xk)
		call dfunc (dfk,xk)
		
		xk1 = xk - (fk/dfk)
		
		erro = dabs (xk1-xk)
		
		k=k+1
		xk=xk1
		
		write(*,*) 'k= ' , k, 'xk = ', xk, 'fk = ' , fk, 'erro = ' , erro 
		
	end do
	
	
end program


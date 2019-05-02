subroutine func (fx,x)
implicit none
real*8 :: x,fx

fx= 2.0 + 4.0*x
end subroutine

program main 
	implicit none 
	real*8::b,a
	a=1.0
	call func(b,a) 
	write(*,*) b

end program

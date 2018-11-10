subroutine Euler(x_in, y_in, h, y_out)
	implicit none
	real(8), intent(in)  :: x_in, y_in, h
	real(8), intent(out) :: y_out
	y_out = y_in + f(x_in, y_in)*h
	return
	contains
	real(8) function f(x, y)
		real(8), intent(in) :: x, y
		f = 
		return
	end function f
end subroutine Euler
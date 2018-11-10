subroutine RK2(x_in, y_in, h, y_out)
	implicit none
	real(8), intent(in)  :: x_in, y_in, h
	real(8), intent(out) :: y_out
	real(8)              :: k1, k2
	k1    = h*f(x_in, y_in)
	k2    = h*f(x_in+h/2.0d0, y_in+k1/2.0d0)
	y_out = y_in + k2
	return
	contains
	real(8) function f(x, y)
		real(8), intent(in) :: x, y
		f = 
		return
	end function f
end subroutine RK2
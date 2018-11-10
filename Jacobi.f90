subroutine Jacobi(N, A, eigenvalue)
	implicit none
	integer, intent(in)     :: N
	real(8), intent(inout)  :: A(1:N,1:N)
	real(8), intent(out)    :: eigenvalue(1:N)
	integer                 :: pivot_Rows, pivot_Columns, i
	real(8)                 :: Uans(1:N,1:N), total, pivot_a, phi, t, c, s, tau, temp_mat1, temp_mat2
	real(8),parameter       :: ep = 1.0d-6

	do i = 1, N
		Uans(1:i-1, i) = 0.0d0
		Uans(i, i)     = 1.0d0
		Uans(i+1:N, i) = 0.0d0
	end do
	
	call component_sum(N, A, total)
	do while(total > ep)
	
	
		call pivot_max(N, A, pivot_a, pivot_Rows, pivot_Columns)

		
		phi = 0.50d0*(A(pivot_Columns, pivot_Columns)-A(pivot_Rows, pivot_Rows)) / pivot_a
		t = sign(1.0d0, phi) / (abs(phi) + sqrt(phi**2 + 1.0d0))
  
  
		c = 1.0d0 / sqrt(t**2 + 1.0d0)
		s = t*c
		tau = s/(1.0d0 + c)

		
		A(pivot_Rows, pivot_Rows)       = A(pivot_Rows, pivot_Rows) - t*pivot_a
		A(pivot_Columns, pivot_Columns) = A(pivot_Columns, pivot_Columns) + t*pivot_a
		A(pivot_Rows, pivot_Columns)    = 0.0d0
		A(pivot_Columns, pivot_Rows)    = 0.0d0
		do i = 1, N
			if(i /= pivot_Rows .AND. i /= pivot_Columns) then
				temp_mat1           = A(pivot_Rows, i) 
				temp_mat2           = A(pivot_Columns, i) 
				A(pivot_Rows, i)    = temp_mat1 - s*(temp_mat2 + tau*temp_mat1)
				A(pivot_Columns, i) = temp_mat2 + s*(temp_mat1 - tau*temp_mat2)
				A(i, pivot_Rows)    = A(pivot_Rows, i)
				A(i, pivot_Columns) = A(pivot_Columns, i)
			end if
		end do
		
		do i = 1, N
			temp_mat1              = Uans(i, pivot_Rows)
			temp_mat2              = Uans(i, pivot_Columns)
			Uans(i, pivot_Rows)    = temp_mat1 - s*(temp_mat2 + tau*temp_mat1)
			Uans(i, pivot_Columns) = temp_mat2 + s*(temp_mat1 - tau*temp_mat2)
		end do

		call component_sum(N, A, total)		

	end do
	do i = 1, N
		eigenvalue(i) = A(i, i)
	end do
	A(1:N, 1:N) = Uans(1:N, 1:N)
	
	contains
	subroutine pivot_max(N, A, temp, temp_Rows, temp_Columns)
		integer, intent(in)  :: N
		real(8), intent(in)  :: A(1:N, 1:N)
		integer, intent(out) :: temp_Rows, temp_Columns
		real(8), intent(out) :: temp
		integer              :: i, j
		temp         = 0.0d0
		temp_Rows    = 2
		temp_Columns = 1
		do i = 1, N
			do j = 1, N
				if((i /= j) .AND. abs(temp) < abs(A(j, i))) then
					temp         = A(j, i)
					temp_Rows    = j
					temp_Columns = i
				end if
			end do
		end do
	end subroutine pivot_max

	subroutine component_sum(N, A, ans)
		integer, intent(in)  :: N
		real(8), intent(in)  :: A(1:N, 1:N)
		real(8), intent(out) :: ans
		integer              :: i, j
		ans = 0.0d0
		do i = 1, N
			do j = 1, N
				if((i /= j)) then
					ans = ans + abs(A(i, j))
				end if
			end do
		end do
	end subroutine component_sum
end subroutine Jacobi
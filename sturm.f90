subroutine sturm(N, Amat, time_max, eps, lambda)
	implicit none
	integer, intent(in)  :: N, time_max
	real(8), intent(in)  :: Amat(1:N, 1:N), eps
	real(8), intent(out) :: lambda(1:N)
	integer              :: k, j, counter, time
	real(8)              :: left, right, f1, f2, f3, temp
	
	temp      = lambda_max(N, Amat)
	right     = temp
	do k = 1, N
		left  = -temp
		do time = 1, time_max
			! 固有値を更新
			lambda(k) = 0.50d0*(left + right)

			! Nを計算
			counter = 0
			f1 = 1.0d0
			f2 = lambda(k) - Amat(1,1)
			if(f1*f2 <= 0.0d0) then
				counter = counter + 1
			end if
			do j = 2, N
				f3 = (lambda(k) - Amat(j,j))*f2 - Amat(j-1,j)**2*f1
				if(f2*f3 <= 0.0d0) then
					counter = counter + 1
				end if
				f1 = f2
				f2 = f3
			end do
			! 探索範囲を更新
			if(counter < k) then
				right = lambda(k)
			else
				left  = lambda(k)
			end if
			
			! 収束判定
			if(ABS(left - right) < eps) then
				exit
			end if
		end do

		if(time == time_max+1) then
			stop "Not converged!"
		end if
	end do
	
	contains
	real(8) function lambda_max(N, Amat)
		integer, intent(in)  :: N
		real(8), INTENT(IN)  :: Amat(1:N, 1:N)
		integer              :: i
		real(8)              :: temp
		lambda_max = max(ABS(Amat(1,1))+ABS(Amat(2,1)), ABS(Amat(N-1,N))+ABS(Amat(N,N)))
		do i = 2, N-1
			temp = ABS(Amat(i-1,i)) + ABS(Amat(i,i)) + ABS(Amat(i+1,i))
			if(lambda_max < temp) then
				lambda_max = temp
			end if
		end do
		return
	end function lambda_max
	
end subroutine sturm
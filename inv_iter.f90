  ! 逆反復法で固有値を改良し、対応する固有ベクトルを求めるサブルーチン
  ! Amat は実3重対角行列を仮定
subroutine inv_iter(N, Amat, lambda, time_max, eps, Eigenvector)
	implicit none
	integer, intent(in)    :: N, time_max
	real(8), intent(in)    :: Amat(1:N, 1:N), eps
	real(8), intent(inout) :: lambda(1:N)
	real(8), intent(out)   :: Eigenvector(1:N, 1:N)

	integer              :: i, j, time
	real(8)              :: res
	real(8), allocatable :: v(:), y(:), Bmat(:, :)
	
	allocate(v(1:N), y(1:N), Bmat(1:N, 1:N))
	do j = 1, N
		! 行列コピー
		Bmat = Amat
		! 固有ベクトルを単位ベクトルで初期化
		v(1)   = 1.0d0
		v(2:N) = 0.0d0

		do time = 1, time_max

			! B = A - lambda*I
			do i = 1, N
				Bmat(i,i) = Amat(i,i) - lambda(j)
			end do

			! B*y = v を y について解く
			call solve_tridiag(N, Bmat, v, y)

			! 固有値を更新
			lambda(j) = lambda(j) + 1.0d0/dot_product(v,y)

			! 固有ベクトルを更新
			v = y/NORM2(y)

			! 残差 |A*v - lambda*v| を計算
			y   = matmul(Amat,v) - lambda(j)*v
			res = NORM2(y);

			! 収束判定
			if(res < eps) then
				exit
			end if

		end do
		
		if(time == time_max+1) then
			stop "inv_iter: Not converged!"
		end if
		
		Eigenvector(1:N, j) = v(1:N)

	end do
	deallocate(v, y, Bmat)
	
	contains
	! 実3重対角行列についての連立一次方程式 A*x = b を解くサブルーチン
	! Amat は実3重対角行列を仮定
	subroutine solve_tridiag(N, Amat, b, x)
		integer, intent(in)  :: N
		real(8), intent(in)  :: Amat(1:N, 1:N), b(1:N)
		real(8), intent(out) :: x(1:N)
		integer              :: i
		real(8)              :: tmp
		real(8), allocatable :: P(:), Q(:)

		allocate(P(N-1), Q(N))

		! 係数 P, Q を漸化式で計算
		P(1) = -Amat(2,1)/Amat(1,1)
		Q(1) = b(1)/Amat(1,1)
		do i = 2, N-1
			tmp = Amat(i,i) + Amat(i-1,i)*P(i-1)
			P(i) = -Amat(i+1,i)/tmp
			Q(i) = (b(i) - Amat(i-1,i)*Q(i-1))/tmp
		end do
		Q(N) = (b(N) - Amat(N-1,N)*Q(N-1))/(Amat(N,N) + Amat(N-1,N)*P(N-1))

		! 連立方程式の解を漸化式で計算
		x(N) = Q(N)
		do i = N-1, 1, -1
			x(i) = P(i)*x(i+1) + Q(i)
		end do

		deallocate(P, Q)

	end subroutine solve_tridiag

end subroutine inv_iter

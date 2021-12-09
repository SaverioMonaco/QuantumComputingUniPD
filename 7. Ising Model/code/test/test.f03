program test
  integer :: ii, N, jj


  ii = 1
  N = 2
  do jj=1,2**N,1
    print*, -2*(modulo( (jj-1)/int(2**(N-ii)),2) ) +1
  end do
end program test

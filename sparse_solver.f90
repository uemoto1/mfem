module sparse_solver
  implicit none
      
  type sparse_matrix_csr
    integer :: ndim, nz
    real(8), allocatable :: value(:)
    integer, allocatable :: icolm(:)
    integer, allocatable :: irptr(:)
  end type sparse_matrix_csr
  
contains
  
  subroutine sparse_matmul(a, x, y)
    implicit none
    type(sparse_matrix_csr), intent(in) :: a
    real(8), intent(in) :: x(:)
    real(8), intent(out) :: y(:)
    integer :: i, ij_min, ij_max, ij, j
    real(8) :: temp
    ! compute matrix-vector multification: y = a x
    do i = 1, a%ndim
      temp = 0.
      ij_min = a%irptr(i)
      ij_max = a%irptr(i+1) - 1
      do ij = ij_min, ij_max
        j = a%icolm(ij)
        temp = temp + a%value(ij) * x(j)
      end do
      y(i) = temp
    end do
    return
  end subroutine sparse_matmul
  
  
  integer function getindex(a, i, j) result(ij_index)
    implicit none
    type(sparse_matrix_csr), intent(in) :: a
    integer, intent(in) :: i, j
    integer :: ij_min, ij_max, ij_mid
    integer :: j_min, j_max, j_mid
    ! j_min: the lower-end of the i-th row
    ij_min = a%irptr(i)
    j_min = a%icolm(ij_min)
    if (j < j_min) then
      ij_index = 0
      return ! out of range
    else if (j == j_min) then
      ij_index = ij_min ! find
      return
    end if
    ! j_max: the upper-end of the i-th row
    ij_max = a%irptr(i+1) - 1
    j_max = a%icolm(ij_max)
    if (j > j_max) then
      ij_index = 0 ! out of range
      return
    else if (j == j_max) then
      ij_index = ij_max ! find
      return
    end if
    ! j_min < j < j_max: use binary search
    do while (1 < (ij_max - ij_min))
      ij_mid = (ij_min + ij_max) / 2
      j_mid = a%icolm(ij_mid)
      if (j > j_mid) then
        ij_min = ij_mid
      else if (j < j_mid) then
        ij_max = ij_mid
      else
        ij_index = ij_mid
        return
      end if
    end do
    ij_index = 0 ! not nonzero element
    return
  end function getindex
    
    
  real(8) function getitem(a, i, j) result(value)
    implicit none
    type(sparse_matrix_csr), intent(in) :: a
    integer, intent(in) :: i, j
    integer :: ij_index
    ! search 
    ij_index = getindex(a, i, j)
    if (0 < ij_index) then
      value = a%value(ij_index)
    else
      value = 0
    end if
    return
  end function getitem
  
  
  subroutine allocate_sparse(ndim, nz, a)
    implicit none
    type(sparse_matrix_csr), intent(inout) :: a
    integer, intent(in) :: ndim
    integer, intent(in) :: nz
    a%ndim = ndim
    a%nz = nz
    allocate(a%value(a%nz))
    allocate(a%irptr(a%ndim + 1))
    allocate(a%icolm(a%nz))
    return
  end subroutine
  
  
  subroutine deallocate_sparse(a)
    implicit none
    type(sparse_matrix_csr), intent(inout) :: a
    if (allocated(a%value)) deallocate(a%value)
    if (allocated(a%irptr)) deallocate(a%irptr)
    if (allocated(a%icolm)) deallocate(a%icolm)
    a%ndim = 0
    a%nz = 0
    return
  end subroutine
  
  
  subroutine empty_sparse_like(src, dst)
    implicit none
    type(sparse_matrix_csr), intent(in) :: src
    type(sparse_matrix_csr), intent(inout) :: dst
    call allocate_sparse(src%ndim, src%nz, dst)
    dst%icolm(:) = src%icolm(:)
    dst%irptr(:) = src%irptr(:)
    return
  end subroutine
  
  
  subroutine imcomplete_lu0(a, lu)
    implicit none
    type(sparse_matrix_csr), intent(in) :: a
    type(sparse_matrix_csr), intent(inout) :: lu
    integer :: i, ij_min, ij_max, ij, j, k_max, ik, k
    real(8) :: temp, l_ii
    do i = 1, lu%ndim
      ij_min = lu%irptr(i)
      ij_max = lu%irptr(i+1)-1
      do ij = ij_min, ij_max
        j = lu%icolm(ij)
        temp = getitem(a, i, j)
        k_max = min(i, j) - 1
        do ik = ij_min, ij_max
          k = a%icolm(ik)
          if (k_max < k) exit
          temp = temp - lu%value(ik) * getitem(lu, k, j)
        end do
        
        if (j <= i) then
          lu%value(ij) = temp
          if (j == i) l_ii = temp
        else
          lu%value(ij) = temp / l_ii
        endif
      end do
    end do
    return
  end subroutine imcomplete_lu0
  
  
  subroutine solve_precondion(lu, x)
    implicit none
    type(sparse_matrix_csr), intent(in) :: lu
    real(8), intent(inout) :: x(:)
    integer :: i, j, ij, ij_min, ij_max
    real(8) :: temp, l_ii
    
    ! Back substitution
    do i = lu%ndim, 1, -1
      temp = x(i)
      ij_min = lu%irptr(i)
      ij_max = lu%irptr(i+1) - 1
      do ij = ij_max, ij_min, -1
        j = lu%icolm(ij)
        if (j > i) then
          temp = temp - lu%value(ij) * x(j)
        else
          exit
        end if
      end do
      x(i) = temp
    end do
    
    ! Forward substitution
    do i = 1, lu%ndim
      temp = x(i)
      ij_min = lu%irptr(i)
      ij_max = lu%irptr(i+1) - 1
      do ij = ij_min, ij_max
        j = lu%icolm(ij)
        if (j < i) then
          temp = temp - lu%value(ij) * x(j)
        else if (j == i) then
          l_ii = lu%value(ij)
        else
          exit
        end if
      end do
      x(i) = temp / l_ii
    end do
    return
  end subroutine solve_precondion
    
  
  integer function gmres_csr(a, b, x, nrst, maxiter, tol, lu)
    implicit none
    type(sparse_matrix_csr), intent(in) :: a
    real(8), intent(in) :: b(:)
    real(8), intent(inout) :: x(:)
    integer,  intent(in) :: nrst, maxiter
    real(8), intent(in) :: tol
    type(sparse_matrix_csr), intent(in), optional :: lu
    integer :: n, m, i, j
    real(8) :: v(a%ndim), q(a%ndim, nrst+1)
    real(8) :: h(nrst+1, nrst), r(nrst+1, nrst), e(nrst+1)
    real(8) :: c(nrst), s(nrst), y(nrst)
    real(8) :: beta, rtemp, d
    logical :: flag
    
    do m = 1, maxiter
       call sparse_matmul(a, x(:), v(:))
       v(:) = b(:) - v(:)
       if (present(lu)) call solve_precondion(lu, v) ! preconditioning
       beta = sqrt(dot_product(v(:), v(:)))
       e(1) = beta
       q(:, 1) = v(:) / beta
       do n = 1, nrst
          ! arnordi process
          call sparse_matmul(a, q(:, n), v(:))
          if (present(lu)) call solve_precondion(lu, v) ! preconditioning
          do j = 1, n
             h(j, n) = dot_product(q(:, j), v(:))
             v(:) = v(:) - h(j, n) * q(:, j)
          end do
          h(n+1, n) = sqrt(dot_product(v(:), v(:)))
          q(:, n + 1) = v(:) / h(n + 1, n)
          
          ! givens rotation
          r(1, n) = h(1, n)
          do j = 2, n
             rtemp = c(j-1) * r(j-1, n) + s(j-1) * h(j, n)
             r(j, n) = -s(j-1) * r(j-1, n) + c(j-1) * h(j, n)
             r(j - 1, n) = rtemp
          end do
          d = sqrt(r(n, n) ** 2 + h(n + 1, n) ** 2)
          c(n) = r(n, n) / d
          s(n) = h(n+1, n) / d
          r(n, n) = c(n) * r(n, n) + s(n) * h(n+1, n)
          e(n+1) = - s(n) * e(n)
          e(n) = c(n) * e(n)
       end do
       
       ! back substitution
       do j = nrst, 1, -1
          y(j) = e(j)
          do i = j + 1, nrst
             y(j) = y(j) - r(j, i) * y(i)
          end do
          y(j) = y(j) / r(j, j)
       end do
       
       ! update x
       do i = 1, nrst
          x(:) = x(:) + y(i) * q(:, i)
       end do
       
       if (abs(e(nrst+1)) < tol) then
         gmres_csr = 0
         return
       end if
    end do
    gmres_csr = maxiter
    return
  end function gmres_csr
end module sparse_solver
      
program test
  use sparse_solver
implicit none
  integer, parameter :: ndim = 1000
  integer, parameter :: nz = 3 * ndim
  real(8), parameter :: pi = 3.14159265
  real(8) :: dx

  type(sparse_matrix_csr) :: a, lu
  real(8) :: b(ndim), x(ndim)
  integer :: i, j, n
  
  b = 0.
  x = 0.
  
  dx = pi / ndim
  n = 1
  a%ndim = ndim
  a%nz = nz
  allocate(a%value(nz), a%icolm(nz), a%irptr(ndim+1))
  do i = 1, ndim
    a%irptr(i) = n
    if (i == 1) then
      b(i) = 1.
      a%value(n) = 1.
      a%value(n+1) = 0.
      a%value(n+2) = 0.
      a%icolm(n) = i
      a%icolm(n+1) = i+1
      a%icolm(n+2) = i+2
    endif
    if ((1 < i) .and. (i < ndim)) then
      a%value(n) = 1.00 / dx ** 2
      a%value(n+1) = 1.00 - 2.00 / (dx ** 2)
      a%value(n+2) = 1.00 / dx ** 2
      a%icolm(n) = i-1
      a%icolm(n+1) = i
      a%icolm(n+2) = i+1
    endif
    if  (i == ndim) then
      b(i) = -1.
      a%value(n) = 0.
      a%value(n+1) = 0.
      a%value(n+2) = 1.
      a%icolm(n) = i-2
      a%icolm(n+1) = i-1
      a%icolm(n+2) = i
    endif
    n = n + 3
  enddo
  a%irptr(ndim+1) = n
  call empty_sparse_like(a, lu)
  call imcomplete_lu0(a, lu)
  


  i = gmres_csr(a, b, x, 50, 1000, 1d-20, lu=lu)
  
  do i = 1, ndim
    write(*,*)  x(i)
  enddo
  
  stop 
end program test
  
  
        
        
      

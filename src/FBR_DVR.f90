module FBR_DVR
use machina_basic
implicit none

real(kind=f8), parameter :: pi = 4.0_f8*atan(1.0_f8)
integer(kind=i4) :: ngrid
real(kind=f8) :: lefta,rightb   ! [a,b]
real(kind=f8) :: weight
real(kind=f8),allocatable,dimension(:) :: grid
real(kind=f8),allocatable,dimension(:,:) ::  Xmn
real(kind=f8),allocatable,dimension(:) :: diagX_eigen

type,public :: FBR 
    real(kind=f8),allocatable,dimension(:,:) :: V_FBR
    real(kind=f8),allocatable,dimension(:,:) :: T_FBR
    real(kind=f8),allocatable,dimension(:,:) :: H_FBR
    real(kind=f8),allocatable,dimension(:,:) :: C_FBR
    real(kind=f8),allocatable,dimension(:) :: E_FBR

    contains
    procedure,public :: FBR_init
    procedure,public :: FBR_matrix_calc
    procedure,public :: diagH_FBR

end type

type,public :: DVR
    real(kind=f8),allocatable,dimension(:,:) :: B
    real(kind=f8),allocatable,dimension(:,:) :: V_DVR
    real(kind=f8),allocatable,dimension(:,:) :: T_DVR
    real(kind=f8),allocatable,dimension(:,:) :: H_DVR
    real(kind=f8),allocatable,dimension(:,:) :: C_DVR
    real(kind=f8),allocatable,dimension(:) :: E_DVR

    contains
    procedure,public :: DVR_init
    procedure,public :: DVR_matrix_calc
    procedure,public :: diagH_DVR
    
end type

contains

subroutine get_ngrid
    implicit none
    integer :: i

    write(*,'(/,a)') 'Enter the number of grids(must > 10):'
    read(*,*) ngrid 
    rightb = 10.0    
    lefta = -10.0
    allocate(grid(ngrid))
    allocate(Xmn(ngrid,ngrid))
    allocate(diagX_eigen(ngrid))

    do i = 1, ngrid
        grid(i) = lefta + i*(rightb-lefta)/(ngrid+1)
    end do
    weight = (rightb-lefta)/(ngrid+1)

end subroutine

function sin_bais(n,x) result(value)
    implicit none 
    integer,intent(in) :: n
    real(kind=f8),intent(in) :: x
    real(kind=f8) :: value

    value = sqrt(2.0_f8/(rightb-lefta))*sin(real(n,kind=f8)*pi*(x-lefta)/(rightb-lefta))

end function

!--------------FBR--------------!
subroutine FBR_init(this)
    implicit none
    class(FBR) :: this 

    allocate(this%V_FBR(ngrid,ngrid))
    allocate(this%T_FBR(ngrid,ngrid))
    allocate(this%H_FBR(ngrid,ngrid))
    allocate(this%C_FBR(ngrid,ngrid))
    allocate(this%E_FBR(ngrid))

end subroutine 



subroutine FBR_matrix_calc(this)
    implicit none 
    class(FBR) :: this
    integer :: i,j,l

    associate(n=>ngrid, T=>this%T_FBR, V=>this%V_FBR, w=>weight, H=>this%H_FBR)
    V = 0.0_f8
    do i = 1, n
        do j = 1, n
            if (i==j) then
                T(i,j) =  0.5_f8*(i*pi/(rightb-lefta))**2
            else
                T(i,j) = 0.0_f8
            end if 
            do l = 1, ngrid
                V(i,j) = V(i,j) + w*sin_bais(i,grid(l))*0.5_f8*grid(l)**2*sin_bais(j,grid(l))
            end do 
        end do 
    end do 
    H = T + V 
    end associate

end subroutine

subroutine diagH_FBR(this)
    implicit none
    class(FBR) :: this 
    real(kind=f8),dimension(3*ngrid-1) :: w
    real(kind=f8),allocatable,dimension(:,:) :: tempH
    integer :: info 

    allocate(tempH(ngrid,ngrid))
    associate(n=>ngrid,H=>this%H_FBR, C=>this%C_FBR, E=>this%E_FBR)
    tempH = H
    call dsyev('V', 'L', n, tempH, n, E, w, 3*n-1, info)
    C = tempH 
    deallocate(tempH)
    end associate
end subroutine

!------------------DVR------------------!
subroutine DVR_init(this)
    implicit none
    class(DVR) :: this 

    allocate(this%B(ngrid,ngrid))
    allocate(this%V_DVR(ngrid,ngrid))
    allocate(this%T_DVR(ngrid,ngrid))
    allocate(this%H_DVR(ngrid,ngrid))
    allocate(this%C_DVR(ngrid,ngrid))
    allocate(this%E_DVR(ngrid))

end subroutine 

subroutine DVR_matrix_calc(this)
    implicit none
    class(DVR) :: this 
    integer :: i,j

    associate(n=>ngrid, T=>this%T_DVR, V=>this%V_DVR, H=>this%H_DVR)
    do i = 1, n 
        do j = 1, n 
            this%B(i,j) = sqrt(2.0_f8/real(n+1,kind=f8))*sin(i*j*pi/real(n+1,kind=f8))
            if (i==j) then 
                T(i,j) = (pi**2)/(4.0_f8*(rightb-lefta)**2) &
                         *((2.0_f8*(real(n+1,kind=f8)**2)+1.0_f8)/3.0_f8 &
                         -1.0_f8/(sin(pi*i/real(n+1,kind=f8))**2))

                V(i,j) = 0.5_f8*grid(i)**2
            else 
                T(i,j) = ((-1.0)**(i-j))*pi**2/(4.0_f8*(rightb-lefta)**2) &
                         *(1.0_f8/(sin(pi*(i-j)/(2.0_f8*real(n+1,kind=f8)))**2) &
                         -1.0_f8/(sin(pi*(i+j)/(2.0_f8*real(n+1,kind=f8)))**2))
                
                V(i,j) = 0.0_f8
            end if 
        end do 
    end do 
    H = T + V
    end associate

end subroutine

subroutine diagH_DVR(this)
implicit none
    class(DVR) :: this 
    real(kind=f8),dimension(3*ngrid-1) :: w
    real(kind=f8),allocatable,dimension(:,:) :: tempH
    integer :: info 

    allocate(tempH(ngrid,ngrid))
    associate(n=>ngrid,H=>this%H_DVR, C=>this%C_DVR, E=>this%E_DVR)
    tempH = H
    call dsyev('V', 'L', n, tempH, n, E, w, 3*n-1, info)
    C = tempH 
    deallocate(tempH)
    end associate
end subroutine


!----------------diag<m|X|n>----------------!
subroutine diagX()
implicit none
    real(kind=f8),dimension(3*ngrid-1) :: w
    integer :: i,j,info 

    Xmn = 0.0_f8
    do i = 1, ngrid
        do j = 1, ngrid
            if (i==j) then
                Xmn(i,j) = 0.5*(lefta**2-rightb**2)
            else 
                Xmn(i,j) = 4.0_f8*(rightb-lefta)*i*j*(-1.0_f8+cos(i*pi)*cos(j*pi))/((i**2-j**2)**2*pi**2)
            end if
        end do
    end do 

    call dsyev('V', 'L', ngrid, Xmn, ngrid, diagX_eigen, w, 3*ngrid-1, info)

end subroutine



end module



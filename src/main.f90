program main
use machina_basic
use FBR_DVR
    implicit none

    type(FBR) :: FBR0 
    type(DVR) :: DVR0
    real(kind=f8) :: E_exact,error_FBR,error_DVR
    integer :: i,j

    call get_ngrid()

    call FBR0%FBR_init()
    call FBR0%FBR_matrix_calc()
    call FBR0%diagH_FBR()

    call DVR0%DVR_init()
    call DVR0%DVR_matrix_calc()
    call DVR0%diagH_DVR()

    write(*,'(/,a)') 'm = 1, omega = 1, interval [-10,10]'
    write(*,'(a)') 'Only show error <= 0.0001' 
    write(*,'(/,a4,5x,a8,5x,a8,5x,a8,5x,a8,5x,a8)') 'n','E_FBR','E_DVR','E_exact','errorFBR','errorDVR'
    write(*,'(a4,5x,a8,5x,a8,5x,a8,5x,a8,5x,a8)') '====','========','========','========','========','========'
    do i = 1, ngrid
        E_exact = 0.5_f8+real(i-1,kind=f8)
        error_FBR = abs(E_exact - FBR0%E_FBR(i))
        error_DVR = abs(E_exact - DVR0%E_DVR(i))
        if (error_DVR <= 0.0001 .or. error_FBR <= 0.0001) then
            write(*,'(i4,5x,f8.5,5x,f8.5,5x,f8.5,5x,f8.5,5x,f8.5)') i-1,FBR0%E_FBR(i),DVR0%E_DVR(i),E_exact,error_FBR,error_DVR
        else 
            exit
        end if
    end do

    call diagX()
    write(*,'(/,a)') 'diag<m|X|n> vs Bml for 10*10 matrix element'
    write(*,'(a)') 'Eigenvector of <m|X|n>: '
    do i = 1, 10
        write(*,'(*(f13.6))') (Xmn(i,j),j=1,10)
    end do 
    write(*,'(/,a)') 'Bml: '
    do i = 1, 10
        write(*,'(*(f13.6))') (DVR0%B(i,j),j=1,10)
    end do 

    write(*,'(/,a)') 'Grid(xl)'
    write(*,'(50(f13.6))') grid
    write(*,'(a)') 'Eigenvalue of <m|X|n>'
    write(*,'(50(f13.6))') diagX_eigen

    open(99,file='FBR_coeff.log',status='replace')
    do i = 1, ngrid
        write(99,'(*(f13.6))') (FBR0%C_FBR(i,j),j=1,ngrid)
    end do 
    open(100,file='DVR_coeff.log',status='replace')
    do i = 1, ngrid
        write(100,'(*(f13.6))') (DVR0%C_DVR(i,j),j=1,ngrid)
    end do 
    close(99)
    close(100)
    
    

end program main

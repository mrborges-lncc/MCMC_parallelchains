program quartile
! Fortran program to test the 3 methods of Quartile Calculation.
! Quartiles are calculated differently in EXCEL, Minitab and SAS method 5
! This program tests the 3 methods.
! by Sukhbinder Singh (aerogeek)
! Background info http://sukhbinder.wordpress.com/2011/01/22/fortran-program-to-calculate-quartiles-as-calculated-by-minitab-and-excel/ 

implicit none
interface
 subroutine sasm5(x,quart)
  real*4 :: x(:),quart
 end subroutine

 subroutine sasm4(x,quart)
  real*4 :: x(:),quart
 end subroutine

 subroutine excel(x,quart)
  real*4 :: x(:),quart
 end subroutine
end interface

real*4 :: a,b,c,x(8),quart,eps
!data x/1.0d0,2.0d0,3.0d0,4.0d0/
data x/-4.50707786E-02,  0.408986449,  -3.95143293E-02, 8.60971287E-02, 0.118125968,0.216884896, 0.246499196, 0.322567642/
integer :: indx(8)

indx = 0
eps = 1e-6
write(*,*)x
call hpsort_eps_epw (8, x, indx, eps)
write(*,*)x
write(*,*)indx
Print *, "Excel Method"
quart=0.25d0
call excel(x,quart)
print *,'  q1=',quart

quart=0.500d0
call excel(x,quart)
print *,'  Median=',quart

quart=0.750d0
call excel(x,quart)
print *,'  q3=',quart



Print *, "Minitab Method"
quart=0.25d0
call sasm4(x,quart)
print *,'  q1=',quart

quart=0.500d0
call sasm4(x,quart)
print *,'  Median=',quart

quart=0.750d0
call sasm4(x,quart)
print *,'  q3=',quart


Print *, "SAS Method 5"
quart=0.25d0
call sasm5(x,quart)
print *,'  q1=',quart

quart=0.500d0
call sasm5(x,quart)
print *,'  Median=',quart

quart=0.750d0
call sasm5(x,quart)
print *,'  q3=',quart


end

subroutine excel(x,quart)
implicit none
! Excel's method to calculate quartiles.
! Based on discussion in this paper http://www.haiweb.org/medicineprices/manual/quartiles_iTSS.pdf
!

real*4 :: x(:),quart,a,b,c
integer :: n,ib

n=size(x)

a=(n-1)*quart
call getgp(a,b,c)

ib=int(c)
!print *,n,a,b,c,ib


quart= (1-b)*x(ib+1) +b*x(ib+2)

end 

subroutine sasm4(x,quart)
implicit none
! Minitab’s method of calculating is the same as the SAS Method 4.
! Based on discussion in this paper http://www.haiweb.org/medicineprices/manual/quartiles_iTSS.pdf
!

real*4 :: x(:),quart,a,b,c
integer :: n,ib


n=size(x)

a=(n+1)*quart
call getgp(a,b,c)

ib=int(c)
!print *,n,a,b,c,ib

if((ib+1)>n) then
quart=(1-b)*x(ib) +b*x(n)
else
quart=(1-b)*x(ib) +b*x(ib+1)
end if

end 


subroutine sasm5(x,quart)
! Calculate Quartiles using SAS Method 5
! This method is the default method of SAS and is based on the empirical distribution function. 
! Based on discussion in this paper http://www.haiweb.org/medicineprices/manual/quartiles_iTSS.pdf
!
implicit none
real*4 :: x(:),quart,a,b,c,tol,diff
integer :: n,ib

tol=1.e-8
n=size(x)

a=n*quart
call getgp(a,b,c)

ib=int(c)
!print *,n,a,b,c,ib

diff=b-0.0d0
if(diff <=tol) then
  quart=(x(ib+1)+x(ib))/2.0d0
else
 quart=x(ib+1)
end if

end subroutine sasm5

subroutine getgp(a,b,c)
! Subroutine to that returns the Right hand and Left hand side digits of a decimal number
real*4 :: a,b,c

b=mod(a,1.0d0)
c=a-b

end subroutine getgp
!                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from flib/hpsort_eps
  !---------------------------------------------------------------------
  subroutine hpsort_eps_epw (n, ra, ind, eps)
  !---------------------------------------------------------------------
  ! sort an array ra(1:n) into ascending order using heapsort algorithm,
  ! and considering two elements being equal if their values differ
  ! for less than "eps".
  ! n is input, ra is replaced on output by its sorted rearrangement.
  ! create an index table (ind) by making an exchange in the index array
  ! whenever an exchange is made on the sorted data array (ra).
  ! in case of equal values in the data array (ra) the values in the
  ! index array (ind) are used to order the entries.
  ! if on input ind(1)  = 0 then indices are initialized in the routine,
  ! if on input ind(1) != 0 then indices are assumed to have been
  !                initialized before entering the routine and these
  !                indices are carried around during the sorting process
  !
  ! no work space needed !
  ! free us from machine-dependent sorting-routines !
  !
  ! adapted from Numerical Recipes pg. 329 (new edition)
  !
  !use kinds, ONLY : DP
  implicit none  
  !-input/output variables
  integer, intent(in)   :: n  
  real(4), intent(in)  :: eps
  integer :: ind (n)  
  real(4) :: ra (n)
  !-local variables
  integer :: i, ir, j, l, iind  
  real(4) :: rra  
!
  ! initialize index array
  IF (ind (1) .eq.0) then  
     DO i = 1, n  
        ind (i) = i  
     ENDDO
  ENDIF
  ! nothing to order
  IF (n.lt.2) return  
  ! initialize indices for hiring and retirement-promotion phase
  l = n / 2 + 1  

  ir = n  

  sorting: do 
  
    ! still in hiring phase
    IF ( l .gt. 1 ) then  
       l    = l - 1  
       rra  = ra (l)  
       iind = ind (l)  
       ! in retirement-promotion phase.
    ELSE  
       ! clear a space at the end of the array
       rra  = ra (ir)  
       !
       iind = ind (ir)  
       ! retire the top of the heap into it
       ra (ir) = ra (1)  
       !
       ind (ir) = ind (1)  
       ! decrease the size of the corporation
       ir = ir - 1  
       ! done with the last promotion
       IF ( ir .eq. 1 ) then  
          ! the least competent worker at all !
          ra (1)  = rra  
          !
          ind (1) = iind  
          exit sorting  
       ENDIF
    ENDIF
    ! wheter in hiring or promotion phase, we
    i = l  
    ! set up to place rra in its proper level
    j = l + l  
    !
    DO while ( j .le. ir )  
       IF ( j .lt. ir ) then  
          ! compare to better underling
          IF ( hslt( ra (j),  ra (j + 1) ) ) then  
             j = j + 1  
          !else if ( .not. hslt( ra (j+1),  ra (j) ) ) then
             ! this means ra(j) == ra(j+1) within tolerance
           !  if (ind (j) .lt.ind (j + 1) ) j = j + 1
          ENDIF
       ENDIF
       ! demote rra
       IF ( hslt( rra, ra (j) ) ) then  
          ra (i) = ra (j)  
          ind (i) = ind (j)  
          i = j  
          j = j + j  
       !else if ( .not. hslt ( ra(j) , rra ) ) then
          !this means rra == ra(j) within tolerance
          ! demote rra
         ! if (iind.lt.ind (j) ) then
         !    ra (i) = ra (j)
         !    ind (i) = ind (j)
         !    i = j
         !    j = j + j
         ! else
             ! set j to terminate do-while loop
         !    j = ir + 1
         ! endif
          ! this is the right place for rra
       ELSE
          ! set j to terminate do-while loop
          j = ir + 1  
       ENDIF
    ENDDO
    ra (i) = rra  
    ind (i) = iind  

  END DO sorting    
contains 

  !  internal function 
  !  compare two real number and return the result

  logical function hslt( a, b )
    REAL(4) :: a, b
    IF( abs(a-b) <  eps ) then
      hslt = .false.
    ELSE
      hslt = ( a < b )
    end if
  end function hslt

  !
end subroutine hpsort_eps_epw

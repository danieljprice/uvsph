module readuv
 !
 ! read the data from the uv data file
 !
 implicit none

contains

subroutine read_uv_data(filename,npts,u,v,re,im,weights,ierr)
 character(len=*), intent(in) :: filename
 integer, intent(out) :: npts,ierr
 real,    intent(out), allocatable :: u(:),v(:),re(:),im(:),weights(:)
 integer :: iu,i,j,npts_orig
 real :: dum,umin,umax,vmin,vmax
 logical, parameter :: hermitian=.false.

 ! open file
 open(newunit=iu,file=filename,status='old',form='formatted',iostat=ierr)
 if (ierr /= 0) then
    print*,'ERROR: could not open '//trim(filename)
    return
 endif

 ! skip header lines
 do i=1,2
    read(iu,*,iostat=ierr)
 enddo

 ! count entries
 npts = 0
 do while (ierr == 0)
    read(iu,*,iostat=ierr)
    npts = npts + 1
 enddo
 npts = npts - 1
 print*,' got npts = ',npts

 npts_orig = npts
 if (hermitian) npts = 2*npts_orig
 ! allocate memory
 allocate(u(npts),v(npts),re(npts),im(npts),weights(npts))
 rewind(iu)
 do i=1,2
    read(iu,*,iostat=ierr)
 enddo

 ! read the data
 print*,' reading data'
 do i=1,npts_orig
    read(iu,*) u(i),v(i),dum,dum,re(i),im(i),weights(i)
 enddo
 close(iu)

 if (hermitian) then
    print*,'hermitian fill-in ->',npts,' points'
    do i=1,npts_orig
       u(npts_orig+i) = -u(i)
       v(npts_orig+i) = -v(i)
       re(npts_orig+i) = re(i)
       im(npts_orig+i) = -im(i)
       weights(npts_orig+i) = weights(i)
    enddo
 endif

 umin = minval(u)
 umax = maxval(u)
 vmin = minval(v)
 vmax = maxval(v)
 print*,' u = [',umin,'->',umax,']'
 print*,' v = [',vmin,'->',vmax,']'

end subroutine read_uv_data

!-----------------------------------------------------------------
!   output pixmap as an ascii file
!-----------------------------------------------------------------
subroutine write_pix(filename,datpix,npixx,npixy,xmin,ymin,xmax,ymax)
 integer, intent(in) :: npixx,npixy
 real,    intent(in), dimension(npixx,npixy) :: datpix
 real,    intent(in) :: xmin,ymin,xmax,ymax
 character(len=*), intent(in) :: filename
 character(len=10) :: stringx,stringy
 character(len=30) :: fmtstring
 integer :: ierr,j
 integer, parameter :: iunit = 166
!
!--write ascii file
!
 open(unit=iunit,file=filename,status='replace',form='formatted',iostat=ierr)
 if (ierr /=0) then
    print*,'error opening '//trim(filename)
    return
 endif
 write(*,"(a)",ADVANCE='NO') '> writing pixel map to file '//trim(filename)//' ...'

 write(stringx,"(i10)") npixx
 write(stringy,"(i10)") npixy
 write(iunit,"(a)",err=66) '# '//trim(adjustl(filename))//' created by uvsph'
 write(iunit,"(a)",err=66) '# Contains 2D pixel array '//trim(adjustl(stringx))//' x '//trim(adjustl(stringy))//' written as '
 write(iunit,"(a)",err=66) '#   do j=1,'//trim(adjustl(stringy))
 write(iunit,"(a)",err=66) '#      write(*,*) dat(1:'//trim(adjustl(stringx))//',j)'
 write(iunit,"(a)",err=66) '#   enddo'
 write(iunit,"(a,1pe14.6,a,1pe14.6)",err=66) '# image: min = ',minval(datpix),' max = ',maxval(datpix)
 write(iunit,"(a,1pe14.6,a,1pe14.6)",err=66) '# x axis: min = ',xmin,' max = ',xmax
 write(iunit,"(a,1pe14.6,a,1pe14.6)",err=66) '# y axis: min = ',ymin,' max = ',ymax
 write(iunit,"(a)",err=66) '# '//trim(adjustl(stringx))//' '//trim(adjustl(stringy))

 write(fmtstring,"(a,i6,a)",iostat=ierr) '(',npixx,'(1pe14.6))'
 if (ierr /= 0) then
    do j=1,npixy
       write(iunit,*,err=66) datpix(1:npixx,j)
    enddo
 else
    do j=1,npixy
       write(iunit,fmtstring,err=66) datpix(1:npixx,j)
    enddo
 endif
 close(iunit)
 print "(a)",'OK'
 return

66 continue
 print "(a)",' ERROR during write '
 close(iunit)
 return

end subroutine write_pix


end module readuv

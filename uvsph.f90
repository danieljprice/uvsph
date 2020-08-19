program uvsph
!
! program to interpolate uv plane data to a set of pixels
!
 use readuv,           only:read_uv_data,write_pix
 use interpolations2D, only:interpolate2D_pixels
 integer, parameter :: nx = 512
 integer :: npts,ierr
 real, allocatable :: u(:),v(:),re(:),im(:)
 integer, allocatable :: mask(:)
 real :: image_real(nx,nx),image_im(nx,nx),rho(nx,nx),r2,uvtaper
 character(len=120) :: uvfile

! get command line arguments
 nargs = command_argument_count()
 if (nargs /= 1) then
    print "(a)",' Usage: uvsph uvfile.txt'
    stop
 endif
 call get_command_argument(1,uvfile)

! step 1: read uv data from file
 call read_uv_data(uvfile,npts,u,v,re,im,ierr)

! step 2: apply uv taper
 uvtaper = 800.
 print*,' applying uv taper with lam=',uvtaper
 do i=1,npts
    ! apply uv taper
    r2 = u(i)**2 + v(i)**2
    re(i) = re(i)*exp(-r2/(2.*uvtaper**2))
    im(i) = im(i)*exp(-r2/(2.*uvtaper**2))
 enddo

! step 3: set up the uv grid
 umin = -4000.
 umax = 4000.
 hfac = 50. ! h = hfac*(1/n)^0.5, where n = points per pixel^2

! step 4: interpolate visibilities to the set of pixels
! in my first attempts at uv-plane interpolation, I just used a
! smoothing length per point h(r) = max(0.15*r,20.), where r=sqrt(u^2 + v^2)
! In the following routines, we compute the local point density
! in the uv plane and use this to set the smoothing length
!
 allocate(mask(npts))
 mask = 1
 call interpolate2D_pixels(u,v,mask,npts, &
      umin,umin,umax,umax,image_real,nx,nx,&
      normalise=.true.,adaptive=.true.,dat=re,datpix2=rho,fac=hfac)
 call interpolate2D_pixels(u,v,mask,npts, &
      umin,umin,umax,umax,image_im,nx,nx,&
      normalise=.true.,adaptive=.true.,dat=im,fac=hfac)
 deallocate(mask)

 ! check total "mass" is conserved by the interpolation, i.e. number of points
 print*,' total number of points is ',sum(rho)

! step 5: write interpolated uv-plane images to file
 call write_pix('uv-img-re.pix',image_real,nx,nx,umin,umin,umax,umax)
 call write_pix('uv-img-im.pix',image_im,nx,nx,umin,umin,umax,umax)
 call write_pix('rho.pix',rho,nx,nx,umin,umin,umax,umax)

 deallocate(u,v,re,im)

end program uvsph

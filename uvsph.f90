program uvsph
!
! program to interpolate uv plane data to a set of pixels
!
 use readuv,           only:read_uv_data,write_pix
 use interpolations2D, only:interpolate2D_pixels
 use system_utils,     only:get_command_flag,get_command_option
 implicit none
 integer, parameter :: nx = 512
 integer :: npts,ierr,i,nargs,j
 integer, allocatable :: mask(:)
 real, allocatable :: u(:),v(:),re(:),im(:),weights(:)
 real :: image_real(nx,nx),image_im(nx,nx),rho(nx,nx),psf(nx,nx)
 real :: r2,uvtaper,hfac,umin,umax,points_per_beam,du,taper,ui,vi
 character(len=120) :: uvfile,string
 logical :: adaptive,get_psf
 real, parameter :: pi = 4.*atan(1.)

! get command line arguments
 nargs = command_argument_count()
 if (nargs < 1) then
    print "(a)",' Usage: uvsph uvfile.txt'
    print "(a)",'        --uvtaper=3000    (apply uv taper)'
    print "(a)",'        --hfac=50         (ratio between beam size and 1/sqrt(points per pixel)^0.5)'
    print "(a)",'        --fixed           (use fixed h, applying hfac)'
    print "(a)",'        --skip_psf        (do not compute psf)'
    stop
 endif
 ! filename is first argument that is not a flag
 do i=1,nargs
    call get_command_argument(1,string)
    if (string(1:1) /= '-') then
       uvfile = string
       exit
    endif
 enddo

 adaptive = .not.get_command_flag('fixed')  ! want adaptive=.true. by default
 get_psf = .not.get_command_flag('skip_psf')              ! psf=false by default
 if (.not.get_psf) print*,' psf = ',get_psf

! step 1: read uv data from file
 call read_uv_data(uvfile,npts,u,v,re,im,weights,ierr)

! step 2: set up the uv grid
 umax = (int(max(-minval(u),maxval(u),-minval(v),maxval(v)))/1000 + 1)*1000.
 umin = -umax
 print*,' using [umin,umax]=',umin,umax

! points_per_beam = 100000. !5000.
! hfac = sqrt(points_per_beam/pi)/2.  ! h = hfac*(1/n)^0.5, where n = number density of points
 hfac = get_command_option('hfac',30.)
 print*,' beam size / uv spacing = ',hfac,' uv points per beam = ',pi*(2.*hfac)**2

! step 3: apply uv taper
 uvtaper = get_command_option('uvtaper',umax/4.)
 print*,' applying uv taper with lam=',uvtaper, ' or t=',uvtaper/206265.,' arcsec'
 do i=1,npts
    ! apply uv taper
    r2 = u(i)**2 + v(i)**2
    taper = exp(-r2/(2.*uvtaper**2))
    re(i) = re(i)*taper
    im(i) = im(i)*taper
 enddo

! step 4: interpolate visibilities to the set of pixels
! In the following routines, we compute the local point density
! in the uv plane and use this to set the smoothing length
!
 allocate(mask(npts))
 mask = 1
 call interpolate2D_pixels(u,v,mask,npts, &
      umin,umin,umax,umax,image_real,nx,nx,&
      normalise=.true.,adaptive=adaptive,dat=re,datpix2=rho,fac=hfac)
 ! check total "mass" is conserved by the interpolation, i.e. number of points
 print*,' number of points (integral of number density) is ',sum(rho),' should be ',npts

 call interpolate2D_pixels(u,v,mask,npts, &
      umin,umin,umax,umax,image_im,nx,nx,&
      normalise=.true.,adaptive=adaptive,dat=im,datpix2=rho,fac=hfac)

 if (get_psf) then
    re = 1.0
    call interpolate2D_pixels(u,v,mask,npts, &
         umin,umin,umax,umax,psf,nx,nx,&
         normalise=.true.,adaptive=adaptive,dat=re,datpix2=rho,fac=hfac)
 endif
 deallocate(mask)

! reapply uv taper in interpolated image
 du = (umax-umin)/nx
 do j=1,nx
    vi = umin + (j-0.5)*du
    do i=1,nx
       ui = umin + (i-0.5)*du
       r2 = ui*ui + vi*vi
       taper = exp(-r2/(2.*uvtaper**2))
       image_real(i,j) = image_real(i,j)*taper
       image_im(i,j) = image_im(i,j)*taper
       psf(i,j) = psf(i,j)*taper
    enddo
 enddo

! step 5: write interpolated uv-plane images to file
 call write_pix('uv-img-re.pix',image_real,nx,nx,umin,umin,umax,umax)
 call write_pix('uv-img-im.pix',image_im,nx,nx,umin,umin,umax,umax)
 call write_pix('rho.pix',rho,nx,nx,umin,umin,umax,umax)
 if (get_psf) call write_pix('psf.pix',psf,nx,nx,umin,umin,umax,umax)

 deallocate(u,v,re,im)

end program uvsph

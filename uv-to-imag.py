#!/usr/bin/env python
#
# Script to inverse fourier transform an "image" in the u-v plane
# to obtain a true image in real space
#
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import os, sys

def read_pix( file ):
    array = np.loadtxt(file,skiprows=9)
    print (file,array.shape)
    return array

def write_fits( img, filename ):
    from astropy.io import fits
    img32 = np.float32(img)
    hdu = fits.PrimaryHDU(img32)
    hdu.writeto(filename,overwrite=True)
    return

def get_uvmap_uvsph():
    # read amplitude and phase from file
    re = read_pix("uv-img-re.pix")
    im = read_pix("uv-img-im.pix")
    c = re + 1j*im
    return c

def image_to_uv( img ):
    # perform Fourier transform
    img_shift = np.fft.fftshift(img)
    uvmap = np.fft.fft2(img_shift)
    uvmap = np.fft.ifftshift(uvmap)
    return uvmap

def uv_to_image( uvmap ):
    # perform inverse Fourier transform
    uv_shift = np.fft.ifftshift(uvmap)
    img = np.fft.ifft2(uv_shift)
    img = np.fft.fftshift(img)
    img = np.real(img)
    return img

def plot_rho():
    img = read_pix("rho.pix")
    print("max density is ",img.max()," points per pixel^2")
    plt.imshow(img)
    #plt.imshow(np.log10(img),vmax=np.log10(img.max()),vmin=-6)
    plt.show()

def plot_psf():
    fig = plt.figure(figsize=(12,6))
    fig.suptitle('psf in uv plane (left) and image plane (right)')
    # plot psf in uv plane (should be just 1 everywhere)
    fig.add_subplot(1,2,1)
    img = read_pix("psf.pix")
    plt.imshow(img,vmin=0.95,vmax=1.05,cmap='RdBu')
    # plot psf in image plane
    img_psf = uv_to_image(img)
    fig.add_subplot(1,2,2)
    plt.imshow(np.log10(img_psf),cmap='inferno')
    plt.savefig('psf.pdf')
    plt.show()

def post_uvsph_pipeline():
    # plot interpolated point density in the uv plane
    plot_rho()

    # plot psf
    #plot_psf()

    # get uv image from file
    c = get_uvmap_uvsph()

    job='hvar-uvsph'
    prefix="channel20-"+job

    # plot it
    fig = plt.figure(figsize=(12,6))
    fig.suptitle('uv plane and channel map - '+job)
    fig.add_subplot(1,2,1)
    uvimg = np.log(np.abs(c))
    plt.imshow(c.imag,cmap='inferno',vmax=1.25,vmin=-1.)
    print(uvimg.min(),uvimg.max())

    # perform the operation
    img = uv_to_image(c)

    # plot the interpolated uv plane and the image plane side by side
    fig.add_subplot(1,2,2)
    plt.imshow(img,cmap='inferno',vmax=img.max())
    print('max intensity is ',img.max())
    #plt.savefig(prefix+'.pdf')
    #plt.savefig(prefix+'.png')
    #plt.show()

    # plot just the image, in a large format
    fig = plt.figure(figsize=(12,12))
    vmin=0. #img.min()
    #vmax=3.e-4
    vmax=img.max()
    plt.imshow(img,cmap='inferno',interpolation='gaussian',vmin=vmin,vmax=vmax)
    plt.savefig('uvsph.png')
    plt.show()
    plt.imshow(np.log10(img),cmap='inferno',interpolation='gaussian',vmin=np.log10(img.max())-4)
    plt.show()

    # write fits
    write_fits(img,'uvsph.fits')

def post_uvsph_noninteractive():
    c = get_uvmap_uvsph()          # get uv image from file
    img = uv_to_image(c)           # perform the operation
    write_fits(img,'uvsph.fits')   # write fits

def read_test_img( imgfile ):
    import imageio
    img = imageio.imread(imgfile)
    img = img/255. # normalise intensity to 1.0
    print(img.shape)
    uvimg = image_to_uv(img)
    # plot real and imaginary parts
    re = np.real(uvimg)
    im = np.imag(uvimg)
    vmax= np.max(im)
    #plt.imshow(re,cmap='inferno',vmax=vmax)
    #plt.show()
    #plt.imshow(im,cmap='inferno',vmax=vmax)
    #plt.show()
    img2 = uv_to_image(uvimg)
    plt.imshow(img,cmap='inferno')
    plt.show()
    #write_fits(re,'test-re.fits')
    #write_fits(im,'test-im.fits')
    return re,im

def sample_img_on_uv_points(file,re,im):
    data = np.loadtxt(file)
    u = data[:,0]
    v = data[:,1]
    #umin = np.min(np.min(u),np.min(v))
    umin=-30000
    umax=30000
    npix = len(re[:,1])
    xx = np.linspace(umin,umax,npix)
    #print("npix=",npix,np.shape(re),np.shape(im),len(xx))
    du = (umax-umin)/npix
    for i,u in enumerate(data):
        ui = data[i,0]
        vi = data[i,1]
        #print(ui,vi)
        ipix = np.int((ui-umin)/du)
        jpix = np.int((vi-umin)/du)
        if (ipix > 0 and jpix > 0 and ipix < npix and jpix < npix):
           #data[i,4] = interpolate.interp2d(xx,xx,re,kind='linear')
           #data[i,5] = interpolate.interp2d(xx,xx,im,kind='linear')
           data[i,4] = re[jpix,ipix]
           data[i,5] = im[jpix,ipix]
           data[i,6] = 1.0
        else:
           data[i,4:6] = 0.

    outfile='uv-test-data.txt'
    print("saving output to ",outfile)
    np.savetxt(outfile,data)
    return outfile

def create_test_img(imgfile,tracksfile):
    re, im = read_test_img(imgfile)
    filename = sample_img_on_uv_points(tracksfile,re,im)
    return filename

def run_uvsph(uvfile):
    print("running uvsph... on ",uvfile)
    os.system('uvsph '+uvfile+' --hfac=10. --uvtaper=6000.')
    return

#uvfile = create_test_img('tests/dan-big-bw.png','Elias24_COcube_channel21.txt')
if (len(sys.argv) == 3):
   imgfile = sys.argv[1]
   trackfile = sys.argv[2]
   uvfile = create_test_img(imgfile,trackfile)
   run_uvsph(uvfile)
elif (len(sys.argv) == 2):
   uvfile = sys.argv[1]
   run_uvsph(uvfile)

post_uvsph_noninteractive()

#post_uvsph_pipeline()

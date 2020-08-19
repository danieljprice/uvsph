#!/usr/bin/env python
#
# Script to inverse fourier transform an "image" in the u-v plane
# to obtain a true image in real space
#
import numpy as np
import matplotlib.pyplot as plt
import os
#import cv2
def generate_maps():
    os.system('splash -ndspmhd imag_00000.dat -r 3 -o ascii -dev re.png')
    os.system('splash -ndspmhd imag_00000.dat -r 4 -o ascii -dev im.png')

def read_pix( file ):
    array = np.loadtxt(file,skiprows=9)
    print (file,array.shape)
    return array

def write_fits( img ):
    from astropy.io import fits
    img32 = np.float32(img)
    hdu = fits.PrimaryHDU(img32)
    hdu.writeto('uvimag.fits',overwrite=True)
    return

def get_test_img():
    array = np.zeros((512,512))
    x0 = 256.
    y0 = 256.
    h2 = 256.
    for j in range(0,512):
        yi = j-y0
        for i in range(0,512):
            xi = i-x0
            r2 = xi**2 + yi**2
            if (r2 < 50.):
               array[i][j]= 1.
            else:
               array[i][j]= 0.
            #array[i][j] = np.exp(-r2/h2)
    return array

def get_uvmap():
    # read amplitude and phase from file
    re = read_pix("imag_00000.dat_vx_proj.pix")
    im = read_pix("imag_00000.dat_vy_proj.pix")
    c = re + 1j*im
    return c

def get_uvmap_cont():
    # read amplitude and phase from file
    re = read_pix("cont_00000.dat_vx_proj.pix")
    im = read_pix("cont_00000.dat_vy_proj.pix")
    c = re + 1j*im
    return c

def get_uvmap_uvsph():
    # read amplitude and phase from file
    re = read_pix("uv-img-re.pix")
    im = read_pix("uv-img-im.pix")
    c = re + 1j*im
    return c

def uv_to_image( uvmap ):
    # perform inverse Fourier transform
    uv_shift = np.fft.ifftshift(uvmap)
    img = np.fft.ifft2(uv_shift)
    img = np.fft.fftshift(img)
    img = np.abs(img)
    return img

def run_test():
    img = get_test_img()
    freq_img = np.fft.fft2(img)
    freq_img = np.fft.fftshift(freq_img)
    print(freq_img.dtype,freq_img.shape)
    plt.imshow(np.abs(freq_img))
    plt.show()
    freq_img = np.fft.ifftshift(freq_img)
    img2 = np.fft.ifft2(freq_img)
    plt.imshow(img2.real)
    plt.show()

def plot_rho():
    img = read_pix("rho.pix")
    print("max density is ",img.max()," points per pixel^2")
    plt.imshow(img)
    plt.show()

#run_test()
#exit()

plot_rho()
#generate_maps()

# get uv image from file
#c = get_uvmap()
#c = get_uvmap_cont()
c = get_uvmap_uvsph()

job='hvar-uvsph'
prefix="channel20-"+job

# plot it
fig = plt.figure(figsize=(12,6))
fig.suptitle('uv plane and channel map - '+job)
fig.add_subplot(1,2,1)
uvimg = np.log(np.abs(c))
#plt.imshow(uvimg,cmap='inferno',vmax=1.25,vmin=-8.)
plt.imshow(c.imag,cmap='inferno',interpolation='nearest',vmax=1.25,vmin=-1.)
print(uvimg.min(),uvimg.max())
#uv.savefig(prefix+'-uv.png')
#plt.savefig(prefix+'-uv.png')
#plt.show()

# perform the operation
img = uv_to_image(c)
print(img.shape,img.dtype)

# write fits
#write_fits(img)

# plot it
#img = np.log(img)
fig.add_subplot(1,2,2)
plt.imshow(img,cmap='inferno',vmax=0.0015)
print(img.max())
#plt.imshow(img,vmax=0.003)
#plt.imshow(img,vmin=np.log(1.e-10),vmax=np.log(1.e-6))
#xy.show()
plt.savefig(prefix+'.pdf')
plt.savefig(prefix+'.png')
plt.show()

fig = plt.figure(figsize=(12,12))
plt.imshow(img,cmap='inferno',interpolation='gaussian',vmax=0.0015)
plt.show()

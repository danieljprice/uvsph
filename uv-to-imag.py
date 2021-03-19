#!/usr/bin/env python
#
# Script to inverse fourier transform an "image" in the u-v plane
# to obtain a true image in real space
#
import numpy as np
import matplotlib.pyplot as plt

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

def plot_rho():
    img = read_pix("rho.pix")
    print("max density is ",img.max()," points per pixel^2")
    plt.imshow(img)
    plt.show()

# plot interpolated point density in the uv plane
plot_rho()

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
plt.imshow(c.imag,cmap='inferno',vmax=1.25,vmin=-1.)
print(uvimg.min(),uvimg.max())

# perform the operation
img = uv_to_image(c)

# plot the interpolated uv plane and the image plane side by side
fig.add_subplot(1,2,2)
plt.imshow(img,cmap='inferno',vmax=img.max())
print('max intensity is ',img.max())
plt.savefig(prefix+'.pdf')
plt.savefig(prefix+'.png')
plt.show()

# plot just the image, in a large format
fig = plt.figure(figsize=(12,12))
vmin=0.
vmax=2.e-4
vmax=img.max()
plt.imshow(img,cmap='inferno',interpolation='gaussian',vmin=vmin,vmax=vmax)
plt.savefig('uvsph.png')
plt.show()

# write fits
#write_fits(img)

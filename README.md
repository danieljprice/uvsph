uvsph
=======
A code using kernel interpolation for image reconstruction in radio astronomy using interferometers


Install
-------
To install and compile the code, use::

```
git clone --recurse-submodules https://github.com/danieljprice/uvsph
cd uvsph
make
```

Run
---
To run the code, use::
```
./uvsph PDS70_1channel.txt
python uv-to-imag.py
```

UV data format
--------------
Currently the file format is assumed to be a simple ascii file with the following format::

```
# x,y,z,date,real,im,weight
# u [m]  v [m]  w [m]  date  real [Jy]  img [Jy]  weight
315.460082333 115.873969448 -257.28136949 4977934957.58 -1.76917896396 0.761573538228 0.470307826996
53.5598112749 -35.1843797917 -32.183281389 4977934957.58 2.15451703883 0.0258628865339 0.552058488131
```

i.e. where u and v points are in columns 1 and 2, and the real and imaginary visibility data is in columns 5 and 6

Daniel Price (2020)

from astropy.io import fits
import re
import numpy

def readfrom(filename):
    # Open the file and get HDU list
    hdulist = fits.open(filename)

    # Column-oriented FITS tables have empty primary HDU
    # The actual data is located inside the first HDU
    if len(hdulist) == 1:
        raise RuntimeError('this is not a column-oriented FITS table')
    hdr = hdulist[1].header
    data = hdulist[1].data[0]
    colnames = hdulist[1].data.dtype.names

    ncol = len(data)
    tbl = dict()

    # Create columns one by one
    for i in range(0, ncol):
        # Store that into the dictionary
        tbl[colnames[i].lower()] = data[i]

    return tbl

def _get_format_code(v):
    vtype = type(v[0])
    if vtype is numpy.ndarray:
        tmp = v
        vlen = 1
        vdim = []
        while vtype is numpy.ndarray:
            vdim.append(len(tmp))
            vlen = vlen*len(tmp)
            tmp = tmp[0]
            vtype = type(tmp)
    else:
        vlen = len(v)
        vdim = [vlen]

    if vtype is numpy.float32:
        ftype='E'
    elif vtype is numpy.float64:
        ftype='D'
    elif vtype is numpy.uint8:
        ftype='B'
    elif vtype is numpy.int32:
        ftype='J'
    elif vtype is numpy.int64:
        ftype='K'
    else:
        raise RuntimeError('unsupported type '+str(vtype))

    return str(vlen)+ftype, '('+','.join(str(d) for d in reversed(vdim))+')'

def writeto(filename, data, output_verify='exception', clobber=False, checksum=False):
    # Build the ColDef array
    cols = fits.ColDefs([])
    for colname, value in data.iteritems():
        # Figure out which FITS format to use to store this column
        vformat, vdim = _get_format_code(value)

        # Add new column to the list
        cols.add_col(fits.Column(
            name=colname.upper(),
            format=vformat,
            dim=vdim))

    # Create the FITS table in memory
    hdu = fits.BinTableHDU.from_columns(cols, nrows=1, fill=True)
    for colname, value in data.iteritems():
        hdu.data[colname.upper()] = value
    # Write it on the disk
    hdu.writeto(filename, output_verify=output_verify, clobber=clobber, checksum=checksum)



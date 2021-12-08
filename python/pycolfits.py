from astropy.io import fits
import numpy as np
import re

def readfrom(filename, lower_case=False, upper_case=False, **kwargs):
    # Open the file and get HDU list
    hdulist = fits.open(filename, **kwargs)

    # Column-oriented FITS tables have empty primary HDU
    # The actual data is located inside the first HDU
    if len(hdulist) == 1:
        raise RuntimeError('this is not a column-oriented FITS table')

    hdr = hdulist[1].header
    data = hdulist[1].data[0]

    # Create columns one by one
    tbl = dict()
    for c in hdulist[1].columns:
        # Store that into the dictionary
        if lower_case:
            cname = c.name.lower()
        elif upper_case:
            cname = c.name.upper()
        else:
            cname = c.name

        # Check the length of columns to avoid any empty ones,
        # as astropy is not able to read empty columns.
        if re.search(r'[(,]0[),]', c.dim) is None:
            tbl[cname] = data[c.name]
        else:
            tbl[cname] = np.array(())

    return tbl

def _get_format_code(v):
    # Convert numpy type code into FITS
    if isinstance(v, str):
        dims = ()
        tdtype = 'S'+str(len(v))
    else :
        if isinstance(v, np.ndarray):
            dims = v.shape
        else:
            dims = ()

        tdtype = v.dtype.str.replace('>','').replace('<','').replace('|','')

        if isinstance(v, np.chararray):
            tdtype = tdtype.replace('U', 'S')

    if tdtype[0] == 'S':
        # Strings are handled in a peculiar way in FITS
        # They are stored as multi-dimensional array of bytes
        # and the last dimension is the maximum length of all
        # strings in the array
        strlen = int(tdtype[1:])
        tdtype = 'a'
        dims = dims + (strlen,)

    if len(dims) == 0:
        repeat = ''
        dims = (1,)
    else:
        repeat = str(np.prod(dims))

    tform = repeat+fits.column.NUMPY2FITS[tdtype]
    tdim = '('+','.join(str(d) for d in reversed(dims))+')'

    return tform, tdim

def writeto(filename, data, lower_case=False, upper_case=False, **kwargs):
    # Build the ColDef array
    cols = fits.ColDefs([])
    for colname, value in data.items():
        # Figure out which FITS format to use to store this column
        tform, tdim = _get_format_code(value)

        if lower_case:
            cname = colname.lower()
        elif upper_case:
            cname = colname.upper()
        else:
            cname = colname

        # Add new column to the list
        cols.add_col(fits.Column(
            name=cname,
            format=tform,
            dim=tdim))

    # Create the FITS table in memory
    hdu = fits.BinTableHDU.from_columns(cols, nrows=1, fill=True)
    for colname, value in data.items():
        if lower_case:
            cname = colname.lower()
        elif upper_case:
            cname = colname.upper()
        else:
            cname = colname

        hdu.data[cname] = value

    # Write it on the disk
    hdu.writeto(filename, **kwargs)



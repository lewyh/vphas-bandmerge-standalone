import sys
import warnings
import dateutil.parser
import numpy as np
from astropy.io import fits

__author__ = 'hfarnhill'


def print_colour_warning(message=None):
    print("+----------------------------------------------------------------+")
    if message is not None:
        lenzero = 63 - len(message)
        print('| {0}{1:<{2}}|'.format(message, ' ', lenzero))
    print("| Please supply EITHER:                                          |")
    print("| * a full or partial set of RED concatenation source lists, OR  |")
    print("| * a full or partial set of BLUE concatenation source lists.    |")
    print("|                                                                |")
    print("| Do NOT:                                                        |")
    print("| * Mix & match from red and blue concatenations                 |")
    print("+----------------------------------------------------------------+")


def runmerge(filelist, radius, aperture):
    if len(filelist) > 8:
        print_colour_warning()
        sys.exit(0)

    # Create empty table which will store key FITS header information
    dtype = [('Filename', '|S256'), ('FieldID', '|S10'), ('obsdate', '|S8'), ('Filter', '|S6'), ('Exp', int),
             ('concat', '|S3')]
    table = np.zeros((len(filelist)), dtype=dtype)

    # Run through the supplied filenames, populating the table with header keywords
    for i, fn in enumerate(filelist):
        table[i]['Filename'] = fn
        h = fits.getheader(fn, 0)
        table[i]['FieldID'] = h["OBJECT"].rstrip()
        table[i]['Filter'] = h["HIERARCH ESO INS FILT1 NAME"].rstrip()
        table[i]['Exp'] = h["HIERARCH ESO TPL EXPNO"]

        obsdate = dateutil.parser.parse(h['DATE'])
        table[i]['obsdate'] = "{0}{1:02}{2:02}".format(obsdate.year, obsdate.month, obsdate.day)

        # Use fact that OBS NAME follows format "p88vphas_1294_hhna" where first character
        # after field identifier can distinguish between red and blue concatenations
        concat = h["HIERARCH ESO OBS NAME"].split('_')[2][0]
        if concat == 'h':
            table[i]['concat'] = 'red'
        elif concat == 'u':
            table[i]['concat'] = 'blu'

    # Ensure that red and blue r-band exposures are not being mixed.
    if len(np.unique(table['concat'])) > 1:
        print_colour_warning("ERROR: Mix of red and blue concat r-band filters.")
        sys.exit(0)

    # Check that only one unique VPHAS field is being bandmerged
    if len(np.unique(table['FieldID'])) > 1:
        print_colour_warning("ERROR: Mix of VPHAS fields.")
        sys.exit(0)

    # Check that there are not too many source lists from any given filter, and that
    # no pointing is repeated (e.g. no 2x r_1 - should be r_1 and r_2)
    if table[0]['concat'] == 'red':
        filter_set = ['r_SDSS', 'i_SDSS', 'NB_659']
    else:
        filter_set = ['r_SDSS', 'u_SDSS', 'g_SDSS']
    filter_max = [2, 2, 3]
    for i, filt in enumerate(filter_set):
        mask = np.where(table['Filter'] == filt)
        if len(table[mask]) > filter_max[i]:
            print_colour_warning("ERROR: Too many {0}-band source lists".format(filt))
            sys.exit(0)
        exposures = np.unique(table[mask]['Exp'], return_counts=True)
        if np.any(exposures[1] > 1):
            print_colour_warning("ERROR: Repeated {0}-band pointing source lists".format(filt))
            sys.exit(0)

    for row in table:
        print(
            "Calculating magnitudes from fluxes inside APERTURE {0} for field {1}, filter {2}, exposure {3}.""".format(
                aperture, row[1], row[3], row[4]))
        single_band(row, aperture)

    # TODO Merge: put together a STILTS command that uses the right column names
    # TODO Launch a subprocess which runs stilts (check java works first!)
    # TODO Clean up intermediate files


def single_band(fieldinfo, aperture):
    outfn = "{1}_{2}_{3}_{4}.fits".format(*fieldinfo)
    filts = {'NB_659': 'Ha', 'r_SDSS': 'r', 'i_SDSS': 'i', 'g_SDSS': 'g', 'u_SDSS': 'u'}
    filtnam = filts[fieldinfo[3]]
    expno = fieldinfo[4]

    f = fits.open(fieldinfo[0])
    EXPTIME = float(f[0].header['EXPTIME'])

    totrows = 0
    for ccd in range(1, 33):
        totrows += int(f[ccd].header["naxis2"])

    outcols = [59, 60, 61, 58, 55, 63, 64, 65, 66]
    colnames = ['RA', 'DEC', "Class_{0}_{1}".format(filtnam, expno), "Av_conf_{0}_{1}".format(filtnam, expno),
                "badpix_{0}_{1}".format(filtnam, expno), "CCD_{0}_{1}".format(filtnam, expno),
                "OID_{0}_{1}".format(filtnam, expno), "{0}_{1}".format(filtnam, expno),
                "err_{0}_{1}".format(filtnam, expno)]
    colunits = ["RADIANS", "RADIANS", "Flag", "Number", "Number", "Number", "Number", "mag", "mag"]
    colp = []
    for i in range(len(outcols)):
        outcols[i] -= 1
        f[1].columns[outcols[i]].name = colnames[i]
        f[1].columns[outcols[i]].unit = colunits[i]
        colp.append(f[1].columns[outcols[i]])
    nrows = 0
    for ccd in range(1, 33):
        ZP = float(f[ccd].header["MAGZPT"])
        APCOR = float(f[ccd].header["APCOR{0}".format(aperture)])
        data = f[ccd].data
        rows = f[ccd].header["naxis2"]
        for i in range(rows):
            data[i][62] = float(ccd)
            data[i][63] = float(nrows + i + 1)
            if data[i][23] > 0.:
                data[i][64] = ZP - APCOR - 2.5 * np.log10(data[i][23] / EXPTIME)
                magerr = 2.5 * np.log10(1. + data[i][24] / data[i][23])
                data[i][65] = max(0.001, magerr)
            else:
                data[i][64] = 99.99
                data[i][65] = 99.99
        if ccd == 1:
            newtab = fits.new_table(colp, nrows=totrows)
        else:
            for i in range(len(outcols)):
                newtab.data.field(i)[nrows:nrows + rows] = f[ccd].data.field(outcols[i])
        nrows = nrows + rows

    # fill in information for new columns
    newtab.header["TTYPE6"] = "CCD_{0}_{1}".format(filtnam, expno)
    newtab.data.names[5] = "CCD_{0}_{1}".format(filtnam, expno)
    newtab.header["TUNIT6"] = "Number"
    newtab.header["TTYPE7"] = "OID_{0}_{1}".format(filtnam, expno)
    newtab.data.names[6] = "OID_{0}_{1}".format(filtnam, expno)
    newtab.header["TUNIT7"] = "Number"
    newtab.header["TTYPE8"] = "{0}_{1}".format(filtnam, expno)
    newtab.data.names[7] = "{0}_{1}".format(filtnam, expno)
    newtab.header["TUNIT8"] = "mag"
    newtab.header["TTYPE9"] = "err_{0}_{1}".format(filtnam, expno)
    newtab.data.names[8] = "err_{0}_{1}".format(filtnam, expno)
    newtab.header["TUNIT9"] = "mag"
    # rename copied columns
    newtab.header["TTYPE3"] = "Class_{0}_{1}".format(filtnam, expno)
    newtab.data.names[2] = "Class_{0}_{1}".format(filtnam, expno)
    newtab.header["TTYPE4"] = "Av_conf_{0}_{1}".format(filtnam, expno)
    newtab.data.names[3] = "Av_conf_{0}_{1}".format(filtnam, expno)
    newtab.header["TTYPE5"] = "badpix_{0}_{1}".format(filtnam, expno)
    newtab.data.names[4] = "badpix_{0}_{1}".format(filtnam, expno)

    # copy primary HDU, add file name information and merge with table
    newhdu = fits.PrimaryHDU()
    newhdu.header = f[0].header
    newhdu.header.set("CASUFILE", fieldinfo[0], comment="CASU File Name", after="ARCFILE")
    newhdu.header.set("VPHAFILE", outfn, comment="VPHAS File Name", after="CASUFILE")

    # copy information from first header
    clin = f[1].header.cards
    for h in (
            "MAGZPT", "MAGZRR", "EXTINCT", "APCOR3", "MED_PA", "NEBULISD", "CROWDED", "APASSZPT", "APASSZRR",
            "APASSNUM"):
        newhdu.header.set(clin[h].keyword, clin[h].value, clin[h].comment)

    # Obtain ellipticity and seeing for each CCD append to header.
    # Calculate average value over whole field, append to header
    ellipticity = np.zeros(32)
    seeing = np.zeros(32)
    skylevel = np.zeros(32)
    for i in range(1, 33):
        ellipticity[i - 1] = f[i].header['ELLIPTIC']
        seeing[i - 1] = f[i].header['SEEING']
        skylevel[i - 1] = f[i].header['SKYLEVEL']
        newhdu.header.set("SEE_{0}".format(i), f[i].header['SEEING'],
                          comment="Average pixel FWHM from CCD{0}".format(i))
        newhdu.header.set("ELL_{0}".format(i), f[i].header['ELLIPTIC'],
                          comment="Average ellipticity from CCD{0}".format(i))
        newhdu.header.set("SKY_{0}".format(i), f[i].header['SKYLEVEL'],
                          comment="Median sky brightness from CCD{0}".format(i))
    newhdu.header.set("SEEING", np.mean(seeing) * 0.21, comment="Average FWHM (arcsec)")
    newhdu.header.set("ELLIPTIC", np.mean(ellipticity), comment="Average ellipticity")
    newhdu.header.set("SKYLEVEL", np.median(skylevel), comment="Average sky level (counts/pixel)")

    current_module = sys.modules['bandmerge']
    newhdu.header.set("HISTORY", "created with vphas-bandmerge-standalone v" + current_module.__version__)

    warnings.filterwarnings('ignore', category=fits.verify.VerifyWarning, append=True)
    hdulist = fits.HDUList([newhdu, newtab])
    # verify output table and save
    hdulist.writeto(outfn, output_verify='silentfix')


def process(args=None):
    import argparse
    import glob
    import sys

    class ThrowingArgumentParser(argparse.ArgumentParser):
        # pass
        def error(self, message):
            print("+----------------------------------------------------------------+")
            print("| Please supply EITHER:                                          |")
            print("| (-f) a list of filenames of source lists, OR                   |")
            print("| (-d) a directory containing ONLY the source lists to be merged |")
            print("|                                                                |")
            print("| Optionally:                                                    |")
            print("| (-r) Crossmatching radius in arcseconds (default 0.5)          |")
            print("| (-a) Aperture in which to measure flux (int: 1-7) (default 3)  |")
            print("+----------------------------------------------------------------+")
            sys.exit(0)

    parser = ThrowingArgumentParser(description='Bandmerge VPHAS+ source lists.')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-f', '--files', action='store', nargs='+', dest='files', help='Filenames of source lists.')
    group.add_argument('-d', '--dir', action='store', nargs=1, dest='dir',
                       help='Path of directory containing ONLY the source lists to be merged.')
    parser.add_argument('-r', '--radius', action='store', nargs=1, dest='radius', default=0.5, type=float,
                        help='Crossmatching radius (arcsec). Default = 0.5.')
    parser.add_argument('-a', '--aperture', action='store', nargs=1, dest='aperture', default=3, type=int,
                        help='Aperture in which to measure flux. Default = 3. Defined at www.vphas.eu')
    args = parser.parse_args(args)

    if not 0 < args.aperture < 8:
        parser.error("")

    if args.files is None and args.dir is not None:
        filelist = glob.glob("{0}*.fits".format(args.dir[0]))
        runmerge(filelist, args.radius, args.aperture)
    elif args.files is not None and args.dir is None:
        runmerge(args.files, args.radius, args.aperture)


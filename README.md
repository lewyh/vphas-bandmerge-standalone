# vphas-bandmerge-standalone
Standalone script that converts reduced VPHAS+ source lists obtained from the ESO archive into bandmerged catalogues.

Prerequisites:
* numpy (v1.9 or greater)
* astropy (v0.3 or greater)

* Java

This script can be used to produce only single-band catalogues (i.e. ignore the bandmerging step), produce both single-band and band-merged catalogues, or only a band-merged catalogue.

Either a list of files, or a directory path which contains only the files to be bandmerged should be provided.

The radii of apertures used to measure flux from VPHAS+ images are defined relative to a "core" radius of 1 arcsecond. For apertures 1-7 the radii are:  
* 1 : 1/2 * rcore  
* 2 : 1/sqrt(2) * rcore  
* 3 : rcore  
* 4 : sqrt(2) * rcore  
* 5 : 2 * rcore  
* 6 : 2*sqrt(2) * rcore  
* 7 : 4 * rcore  
By default, magnitudes corresponding to the flux measured in aperture 3 are returned.

The crossmatching radius determines the maximum distance (in arcsec) between two sources in different filters before they will not be associated.

optional arguments:  
  -h, --help            show this help message and exit  
  -f FILES [FILES ...]  Filenames of source lists.  
  -d DIR                Path of directory containing source lists for merging.  
  -a APERTURE           Aperture # for measuring flux. Default = 3.  
  -r RADIUS             Crossmatching radius (arcsec). Default = 0.5.  
  -x XMX                Memory (in MB) for crossmatching. Default = 2048.  
  -k                    Keep intermediate single-band catalogues.  
  -n                    No bandmerging step (only generate single-band files)  
  -c                    Overwrite pre-existing files (no warnings!)
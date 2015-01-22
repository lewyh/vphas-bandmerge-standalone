__author__ = 'hfarnhill'


def process(args=None):
    import argparse
    import glob
    import sys

    class ThrowingArgumentParser(argparse.ArgumentParser):
        def error(self, message):
            print("+----------------------------------------------------------------+")
            print("| Please supply EITHER:                                          |")
            print("| (-f) a list of filenames of source lists, OR                   |")
            print("| (-d) a directory containing ONLY the source lists to be merged |")
            print("+----------------------------------------------------------------+")
            sys.exit(0)

    parser = ThrowingArgumentParser(description='Bandmerge VPHAS+ source lists.')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-f', '--files', action='store', nargs='+', dest='files', help='Filenames of source lists.')
    group.add_argument('-d', '--dir', action='store', nargs=1, dest='dir',
                       help='Path of directory containing ONLY the source lists to be merged.')
    args = parser.parse_args(args)

    if args.files is None and args.dir is not None:
        filelist = glob.glob(args.dir)
        runmerge(filelist)
    elif args.files is not None and args.dir is None:
        runmerge(args.files)

import sys
process(sys.argv)
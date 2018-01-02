#! /usr/bin/env python3

from collections import namedtuple
from operator import sub
from functools import reduce
import xml.etree.ElementTree as ET
import argparse

# saf are xml from MA5
# convert to text for gnuplot and generic plotting programs

# parser
parser = argparse.ArgumentParser(description="create a run file from a \
           directory containing input files")

parser.add_argument('input',
                    type=str,
                    action='store',
                    help="saf file to process")
parser.add_argument('-o', '--output',
                    dest='output',
                    type=str,
                    action='store',
                    help="currently unimplemented")
parser.add_argument('-x', '--xsec',
                    type=float,
                    action='store',
                    default=1.0,
                    help="Parton level cross section to be used in scaling")
parser.add_argument('-n', '--nevents',
                    type=int,
                    action='store',
                    default=1,
                    help="Number of events to be used in scaling")
parser.add_argument('-b', '--rebin',
                    action='store_true',
                    help="""Should the bins be scaled so they integrate to the
                    total cross section""")
parser.add_argument('--fb',
                    action='store_true',
                    help="""Convert final cross section into fb from pb""")

args = parser.parse_args()

# hardcode test
# test = "/home/luke/Documents/Physics/Research/saf2text/test/histos.saf"

# with open(test) as f:
    # xml = f.read()
# root = ET.fromstring("<root>" + xml + "</root>")

# simple container for data we want
class safhisto(namedtuple("saf_histogram", "obs bins xsec")):
    def __str__(self):
        # header
        pretty_histo = "# " + self.obs + "\n"
        pretty_histo += "# {}    {}\n".format("binmid [GeV]", "xsec [pb]")

        # body
        for bin, xsec in zip(self.bins, self.xsec):
            pretty_histo += "{:6.1f}    {:12.14f}\n".format(bin, xsec)

        return pretty_histo

# structure of histogram

# histogram
#   description
#   statistics
#   data

def histo(histtree, sigma=1., nevents=1, rebin=False, fb=False):
    """"""
    des, stat, dat = histtree

    obs_name, bins, binsize = description(des)
    stat = statistics(stat)
    uflow, bdata, oflow = data(dat)

    # rescale bdata
    rescale = sigma/nevents
    if rebin:
        rescale = rescale / binsize
    if fb:
        # conversion factor for pb -> fb
        rescale = rescale * 1000.

    bdata = [dat * rescale for dat in bdata]

    # return obs_name, bins, bdata
    return safhisto(obs_name, bins, bdata)

def description(elem, bin_alignment="mid"):
    """"""
    obs, _, binning, *other = elem.text.strip().split("\n")

    # observable cleanup
    obs = obs.strip("\"")
    obs = obs.strip("\'")

    # binning
    nbins, xmin, xmax = binning.split()
    # cast to appropriate form
    nbins = int(nbins)
    xmin = float(xmin)
    xmax = float(xmax)

    binsize = (xmax - xmin)/nbins
    if bin_alignment == "left":
        bin_fn = lambda x: x
    elif bin_alignment == "mid":
        bin_fn = lambda x: x+1/2
    elif bin_alignment == "right":
        bin_fn = lambda x: x+1
    else:
        pass

    bins = [xmin + binsize*bin_fn(n) for n in range(nbins)]

    return obs, bins, binsize

def statistics(elem):
    pass

def data(elem):
    """"""
    hline = elem.text.strip().split("\n")

    # handle saf with negative event weights
    try:
        hdata = [reduce(sub, map(float, line.split(" ")[:2])) for line in hline]
    except:
        hdata = [float(line.split(" ")[:1]) for line in hline]

    uflow, *bdata, oflow = hdata

    return uflow, bdata, oflow


if __name__ == "__main__":

    # 
    with open(args.input, "r") as f:
        xml = f.read()
    root = ET.fromstring("<root>" + xml + "</root>")

    for enum, hist_elem in enumerate([elem for elem in root if elem.tag == "Histo"]):
        histogram = histo(hist_elem, sigma=args.xsec, nevents=args.nevents, rebin=args.rebin, fb=args.fb)
        with open("histogram-" + histogram.obs, "w") as f:
            f.write(str(histogram))

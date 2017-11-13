#! /usr/bin/env python3

from collections import namedtuple
import xml.etree.ElementTree as ET

# saf are xml from MA5
# convert to text for gnuplot and generic plotting programs

# hardcode test
test = "/home/luke/Documents/Physics/Research/saf2text/test/histos.saf"

with open(test) as f:
    xml = f.read()
root = ET.fromstring("<root>" + xml + "</root>")

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

def histo(histtree, sigma=1., nevents=1, rebin=False):
    """"""
    des, stat, dat = histtree

    obs_name, bins, binsize = description(des)
    stat = statistics(stat)
    uflow, bdata, oflow = data(dat)

    # rescale bdata
    rescale = sigma/nevents
    if rebin:
        rescale = rescale / binsize

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
    hdata = elem.text.strip().split("\n")
    uflow, *bdata, oflow = [float(line.split(" ")[0]) for line in hdata]

    return uflow, bdata, oflow

# blocks of dependent params from MA5 and MG5
sigma = 0.0
nevents = 50000

if __name__ == "__main__":

    for enum, hist_elem in enumerate([elem for elem in root if elem.tag == "Histo"]):
        histogram = histo(hist_elem)
        print(histogram)
        with open("histogram-" + histogram.obs, "w") as f:
            f.write(str(histogram))

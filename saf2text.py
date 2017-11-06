#! /usr/bin/env python3

from collections import namedtuple
import xml.etree.ElementTree as ET

# saf are xml from MA5
# convert to text for gnuplot

# hardcode test
test = "/home/luke/Documents/Physics/Research/saf2text/test/histos.saf"

with open(test) as f:
    xml = f.read()
root = ET.fromstring("<root>" + xml + "</root>")

# test for a single histogram
hist1 = root[2]

# structure of histogram

# histogram
#   description
#   statistics
#   data

# simple container for data we want
safhisto = namedtuple("saf_histogram", "obs, bins, xsec")

def histo(histtree):
    """"""
    des, stat, dat = histtree

    obs_name, bins = description(des)
    stat = statistics(stat)
    uflow, bdata, oflow = data(dat)

    # return obs_name, bins, bdata
    return safhisto(obs_name, bins, bdata)

def description(elem):
    """"""
    obs, _, binning, *other = elem.text.strip().split("\n")

    nbins, xmin, xmax = binning.split()
    # cast to appropriate form
    nbins = int(nbins)
    xmin = float(xmin)
    xmax = float(xmax)

    binsize = (xmax - xmin)/nbins
    bins = [xmin + binsize*(n+1/2) for n in range(int(nbins))]

    return obs, bins

def statistics(elem):
    pass

def data(elem):
    """"""
    hdata = elem.text.strip().split("\n")
    uflow, *bdata, oflow = [float(line.split(" ")[0]) for line in hdata]

    return uflow, bdata, oflow

# for elem in [elem for elem in root if elem.tag == "Histo"]:
    # print(histo(elem))

h = histo(root[2])

with open("hist", "w") as f:
    # write header
    f.write("# " + h.obs)
    f.write("\n")
    f.write("# bins  xsec")
    f.write("\n")

    # write data
    for bin, dat in zip(h.bins, h.xsec):
        f.write(str(bin) + "  " + str(dat))
        f.write("\n")

#! /usr/bin/env python3

import xml.etree.ElementTree as ET

# saf are xml from MA5
# convert to text for gnuplot

# hardcode test
test = "/home/luke/Documents/Physics/Research/saf2text/test/histos.saf"

with open(test) as f:
    xml = f.read()
root = ET.fromstring("<root>" + xml + "</root>")

for child in root:
    print(child.tag, child.attrib)

hist1 = root[1]
print(hist1.tag)

# structure of histogram

# histogram
#   description
#   statistics
#   data

# this is data
print(hist1[2].text)

# lines
# print(hist1[2].text.strip().split("\n"))
ldata = hist1[2].text.strip().split("\n")

# data from hist 1
hist = []
for line in ldata:
    ld = line.split(" ")[0]
    hist.append(ld)

# now get bins

bdata = hist1[0].text.strip().split("\n")

# line 3 is the important one...
bdata = bdata[2].split()

nbins, xmin, xmax = map(float, bdata)
binsize  = (xmax - xmin)/nbins
bins = [xmin + binsize*n for n in range(int(nbins))]

for bin, dat in zip(bins, hist):
    print(bin, dat)

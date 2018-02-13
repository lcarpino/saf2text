#! /usr/bin/env python3

from collections import namedtuple
from operator import add, sub
from functools import partial, reduce
import itertools
import numbers
from pathlib import Path
import xml.etree.ElementTree as ET
import argparse

# saf are xml from MA5
# convert to text for gnuplot and/or other generic plotting programs

# parser
parser = argparse.ArgumentParser(description="""Covert a Histo.saf file from
MadAnalysis into a set of plain text Histograms for use in an external plotting
program such as gnuplot.""")

parser.add_argument('-e', '--extended',
                    dest='extended',
                    action='store_true',
                    help="""Use extended options. Pass the saf folders rather
                    than the files i.e. ANALYSIS_NAME. This is needed if the
                    showering program does not support IDWTUP = -4 or this
                    option has not been used in the sample passed to
                    MadAnalysis.""")
parser.add_argument('input',
                    metavar='input(s)',
                    type=str,
                    action='store',
                    nargs='+',
                    help=".saf file(s) to be processed")
parser.add_argument('-o', '--output',
                    dest='output',
                    type=str,
                    action='store',
                    default='histogram',
                    help="""The prefix that is attached at the start of the
                    output filename""")
parser.add_argument('-x', '--xsec',
                    type=float,
                    action='store',
                    default=1.0,
                    help="""Parton level cross section to be used to normalise
                    an unweighted histogram""")
parser.add_argument('-n', '--nevents',
                    type=int,
                    action='store',
                    default=1,
                    help="""Number of events to be used to normalise an
                    unweighted histogram""")
parser.add_argument('--avg',
                    action='store_true',
                    help="""Average over the input files, i.e. report
                    histograms as sum over inputs divide the number of inputs.
                    One should use this if runs were performed in parallel, but
                    with the same setup""")
parser.add_argument('--rebin',
                    action='store_true',
                    help="""Rescale the bins so that they integrate to the
                    total cross section rather than the default sum. This is
                    equivalent to converting from nevents/bin to pb/bin""")
parser.add_argument('--fb',
                    action='store_true',
                    help="""Convert final cross section into fb from pb""")


# simple container for the histogram data that we want,
# supports addition of histograms
class safhisto(namedtuple("saf_histogram", "obs bins xsec")):
    def __str__(self):
        pretty_histo = ""
        for bin, xsec in zip(self.bins, self.xsec):
            pretty_histo += "{:08.3f}\t{:014.10f}\n".format(bin, xsec)

        return pretty_histo

    def __add__(self, other):
        # it doesn't make sense to add histograms of different observables
        # or that have different numbers of bins, what is the consistent way
        # to handle these extra bins?
        assert self.obs == other.obs
        assert self.bins == other.bins
        wts = tuple([(x + y) for x, y in zip(self.xsec, other.xsec)])
        return safhisto(self.obs, self.bins, wts)


# structure of histogram output from MadAnalysis
# histogram
#   description
#   statistics
#   data

# routines to process saf files and return safhisto objects
# containing the data in plaintext
def histo(hist_tree, event_wt, rebin, fb):
    """"""
    des, stat, dat = hist_tree

    obs_name, bins, binsize = description(des)
    stat = statistics(stat)
    uflow, bdata, oflow = data(dat)

    # rescale bdata
    rescale = event_wt
    if rebin:
        rescale = rescale / binsize
    if fb:
        # conversion factor for pb -> fb
        rescale = rescale * 1000.

    bdata = [dat * rescale for dat in bdata]

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

    try:
    # handle negative event weights in histogram
        hdata = [reduce(sub, map(float, line.split(" ")[:2])) for line in hline]
    except:
    # usual case, just positive weights
        hdata = [float(line.split(" ")[:1]) for line in hline]

    uflow, *bdata, oflow = hdata

    return uflow, bdata, oflow


# functions to process folder

def handle_ma5_out(path):

    # structure of folder

    # ma5_input_name
    #  ma5_input_name.saf
    #  analysis_name_0
    #    analysis_name.saf
    #    Cutflows
    #      cut_name.saf  # names in analysis_name.saf with "-" -> "_"
    #       .
    #       .
    #      cut_name.saf
    #    Histograms
    #      histos.saf
    #  analysis_name_1
    #    .
    #    .
    #  analysis_name_n

    # ma5_input_name
    ma5_in = path.name
    global_info = (path / ma5_in).with_suffix(".saf")

    # analysis_names
    analyses = [child for child in path.iterdir() if child.is_dir()]
    analysis_names = ["_".join(str(analysis.name).split('_')[:-1]) for analysis in analyses]

    # analysis names should all be the same, but check anyway
    if len(set(analysis_names)) != 1:
        exit()
    else:
        analysis = max(analyses)
        analysis_name = set(analysis_names).pop()

    # analysis_name.saf
    region_selection = (analysis / analysis_name).with_suffix(".saf")

    # Cutflows
    cutflows = [child for child in (analysis / "Cutflows").iterdir()]

    # Histograms
    histograms = (analysis / "Histograms" / "histos").with_suffix(".saf")

    return global_info, region_selection, cutflows, histograms

def global_info(path):
    global_info_xml = readsaf(path)

    saf_header, sample, saf_footer = global_info_xml
    sample_info = list(map(lambda s : s.strip("#").split(), sample.text.strip().split("\n")))
    sample_titles = sample_info[0]
    sample_data = map(lambda x, y : x(y), [float, float, int, float, float], sample_info[1])

    xsec, xsec_err, nevents, sum_poswt, sum_negwt = sample_data

    return xsec, xsec_err, nevents, sum_poswt, sum_negwt

def region_selection(path):
    region_selection_xml = readsaf(path)

    saf_header, regions, saf_footer = region_selection_xml
    region_names = map(lambda s : s.strip("\"").replace("-", "_"), regions.text.split())

    return region_names

def cutflow(path):
    cutflow = readsaf(path)

    saf_header, initial_counter, *counters, saf_footer = cutflow 
    title, *initial_data = initial_counter.text.strip().split("\n")

    title = title.split('#')[0].strip().strip("\"")

    nentries, wt, wt2 = map(lambda x : x.split("#")[0].split(), initial_data)
    nentries = map(int, nentries)
    wt = map(float, wt)
    wt2 = map(float, wt2)

    # num_poswt, num_negwt = nentries

    return nentries, wt, wt2


# helper functions
def readsaf(input):
    """Helper function to read and process saf files
    (saf's are badly formatted xml so need this hack)"""
    # fix to make path objects work in Python3.5, not needed in 3.6+
    with open(str(input), 'r') as f:
        saf = f.read()
    return ET.fromstring("<root>" + saf + "</root>")

def saf2hist(saf, event_wt, rebin, fb):
    """takes a well formatted saf from readsaf and converts it into the
    appropriate histograms"""
    histlist = []
    for enum, hist_elem in enumerate([elem for elem in saf if elem.tag == "Histo"]):
        histlist.append(histo(hist_elem, event_wt, rebin, fb))
    return histlist

def event_weight(args, path):

    regions = region_selection(path)
    nentries, wt, wt2 = cutflow


# main methods
def main_extended(inputs, event_wts, rebin, fb):
    safs = [readsaf(insaf) for insaf in inputs]

    if isinstance(event_wts, numbers.Real):
        # stuff for the old way
        fns = itertools.repeat(partial(saf2hist, event_wt=event_wts, rebin=rebin, fb=fb))
    elif isinstance(event_wts, list):
        # stuff for the extended version
        fns = map(lambda it : partial(saf2hist, event_wt=it, rebin=rebin, fb=fb), event_wts)
    else:
        return NotImplemented

    try:
        # get nested list of all histograms
        all_histos = map(lambda x, y: x(y) , fns, safs)
    except:
        exit()

    return all_histos

if __name__ == "__main__":

    args = parser.parse_args()

    # process inputs, it will let you choose paths that don't exist and
    # continue quietly
    if args.extended:
        inputs = [Path(path) for path in args.input if Path(path).exists()]
        ninputs = len(inputs)

        # process to check
        gi, rs, cfs, hist = zip(*map(handle_ma5_out, inputs))

        xsec, xsec_err, nevents, sum_poswt, sum_negwt = zip(*map(global_info, gi))
        region_names = map(region_selection, rs)
        first_region = itertools.starmap(lambda first, *args : first, region_names)
        cutflow_sample = map(lambda r, cf : next(filter(lambda p : p.stem == r, cf)), first_region, cfs) 
        nentries, wt, wt2 = zip(*map(cutflow, cutflow_sample))

        net_nevents = list(map(lambda it : reduce(sub, it), nentries))

        # perform rescaling
        rescale_nevents = [0.0 if den == 0.0 else num/den for num, den in zip(nevents, net_nevents)]
        # determine if nevents and xsec are needed for rescaling
        wt_test = list(map(lambda x : int(sum(x)), wt))
        if net_nevents == wt_test:
            rescale = [0.0 if nevents == 0.0 else xsec/nevents*scale for xsec, nevents, scale in zip(xsec, nevents, rescale_nevents)]
        else:
            rescale = rescale_nevents

        histos = main_extended(hist, rescale, args.rebin, args.fb)

        if args.avg:
            # weighted average
            tot_evts = sum(nevents)
            avg_wts = [nevts/tot_evts for nevts in nevents]
            histos = [[safhisto(hist.obs, hist.bins, [xsec*avg_wt for xsec in hist.xsec]) for hist in hist_dir]
                      for hist_dir, avg_wt in zip(histos, avg_wts)]

        histos_tp = map(lambda *sl : list(sl), *histos)
        histos = [reduce(add, hist) for hist in histos_tp]

    else:
        inputs = [Path(path) for path in args.input if Path(path).exists()]
        ninputs = len(inputs)
        histos = main_extended(inputs, args.xsec/args.nevents, args.rebin, args.fb)

        histos_tp = list(map(lambda *sl : list(sl), *histos))
        histos = [reduce(add, hist) for hist in histos_tp]

        if args.avg:
            histos = [safhisto(obs, bins, [x/ninputs for x in xsec]) for obs, bins, xsec in histos]

    for hist in histos:
        with open(str(args.output) + "-" + hist.obs, "w") as f:
            f.write(str(hist))

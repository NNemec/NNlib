import tables as pytables
from numarray import array, NewAxis
from units import G_0

def read_legend_values(table):
    legend_by = getattr(table.attrs,'legend_by',None)
    if legend_by:
        return sorted(set([ x[legend_by] for x in table.iterrows() ]))
    else:
        return [ None ]

def read_scan_value(table,legend_value=None):
    legend_by = getattr(table.attrs,'legend_by',None)
    def rows():
        if legend_by and legend_value:
            return table.where(getattr(table.cols,legend_by) == legend_value)
        else:
            return table.iterrows()
    scan_by = table.attrs.scan_by
    energy = array([ x[scan_by] for x in rows() ])
    return energy

def read_energy(table,legend_value):
    legend_by = getattr(table.attrs,'legend_by',None)
    def rows():
        if legend_by:
            return table.where(getattr(table.cols,legend_by) == legend_value)
        else:
            return table.iterrows()

    energy = array([ x['energy'] for x in rows() ])
    return energy

def read_transmission(table,legend_value,read_err=True):
    legend_by = getattr(table.attrs,'legend_by',None)
    def rows():
        if legend_by:
            return table.where(getattr(table.cols,legend_by) == legend_value)
        else:
            return table.iterrows()

    if hasattr(table.cols,'samples'):
        samples = array([ abs(x['samples']) for x in rows() ])

        transmission_sum = array([ x['transmission_sum'] for x in rows() ])
        if len(transmission_sum.shape) > len(samples.shape):
            samples = samples[:,NewAxis]
        transmission_avg = transmission_sum / samples

        if read_err and hasattr(table.cols,'transmission_sqsum'):
            transmission_sqsum = array([ x['transmission_sqsum'] for x in rows() ])
            transmission_var = abs(transmission_sqsum/samples - transmission_avg**2)
            # abs is needed for cases where var==0 mathematically, and therefore possibly
            # var<0 numerically
            transmission_stdev = transmission_var**0.5
            transmission_err = transmission_stdev/(samples**0.5)
        else:
            transmission_err = None

    else:
        transmission_avg = array([ x['transmission'] for x in rows() ])
        if read_err and hasattr(intable.cols,'transmission_err'):
            transmission_err = array([ x['transmission_err'] for x in rows() ])
        else:
            transmission_err = None
    return transmission_avg, transmission_err

def read_samples(table):
    samples = [ x['samples'] for x in table.iterrows() if x['samples'] > 1 ]
    return min(samples),max(samples)

def read_conductance(table,legend_value,read_err=True):
    legend_by = getattr(table.attrs,'legend_by',None)
    def rows():
        if legend_by:
            return table.where(getattr(table.cols,legend_by) == legend_value)
        else:
            return table.iterrows()

    conductance_avg = array([ x['conductance'] for x in rows() ]) * G_0
    if read_err and hasattr(table.cols,'conductance_err'):
        conductance_err = array([ x['conductance_err'] for x in rows() ]) * G_0
    else:
        conductance_err = None
    return conductance_avg, conductance_err



def get_timestamp():
    import time
    return time.strftime("%Y%m%d-%H%M%S",time.localtime())


def get_directory():
    import os
    return os.getenv('DDIR','.')

def get_pathnames(prefix = '.*',suffix = '.h5'):
    import sys, os, re

    namepatterns = []
    pathnames = []

    ddir = get_directory()

    for arg in sys.argv[1:]:
        if arg.endswith('.py'):
            pass
        elif arg.startswith('-'):
            pass
        elif arg.endswith(suffix):
            pathnames += [arg]
        else:
            namepatterns += [arg]

    if pathnames == [] and namepatterns == []:
        namepatterns = [r'.*']

    for namepattern in namepatterns:
        regexp = re.compile(prefix+"-........-......-"+namepattern+suffix.replace('.',r'\.')+"$")
#        regexp = re.compile(prefix+r"-........-......-"+namefragment+suffix.replace('.',r'\.')+"$")
        names = [ n for n in os.listdir(ddir) if regexp.match(n) ]
        if len(names) == 0:
            raise LookupError, "No "+prefix+"-...-" + namepattern + suffix + " file found in directory '" + ddir + "'"
        names.sort()
        pathnames = [ ddir + "/" + name for name in names ]

    return pathnames

def get_pathname(*args,**kwargs):
    names = get_pathnames(*args,**kwargs)
    if len(names) != 1:
        print names
        raise "Exactly one file should be selected!"
    return names[0]

def create_new_file(prefix = None,suffix = '.h5'):
    import sys
    if len(sys.argv) == 1:
        pathname = 'tryout.h5'
        h5file = pytables.openFile(pathname,mode='w')
    else:
        datapath = get_directory()
        filename = get_timestamp() + "-" + sys.argv[1] + ".h5"
        if prefix:
            filename = prefix + '-' + filename
        pathname = datapath + "/" + filename
        h5file = pytables.openFile(pathname,mode='a')
    return h5file

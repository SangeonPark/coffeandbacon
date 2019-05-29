# coding: utf-8
from __future__ import print_function, division
from collections import defaultdict
import gzip
import lz4.frame as lz4f
import cloudpickle as cpkl
import json
import re
import os

import uproot
import numpy as np

from fnal_column_analysis_tools import hist
from fnal_column_analysis_tools.hist import export
from collections import OrderedDict

with lz4f.open("hists_wtag_3.cpkl.lz4") as fin:
    hists_unmapped = cpkl.load(fin)

process = hist.Cat("process", "Process", sorting='placement')
process_cats = ("dataset", "AK8Puppijet0_isHadronicV")
process_map = OrderedDict()

# data
process_map["data_obs"] = ("SingleMuon*",slice(None))
# this unmatched
process_map["catp1"] = (re.compile("(?!SingleMuonRun2017B_17Nov2017_v1)(?!SingleMuonRun2017C_17Nov2017_v1)(?!SingleMuonRun2017D_17Nov2017_v1)(?!SingleMuonRun2017E_17Nov2017_v1)(?!SingleMuonRun2017F_17Nov2017_v1)"), 0)
# this is matched
process_map["catp2"] = ("*", slice(1,None))

hists = {}
for key, val in hists_unmapped.items():
    if isinstance(val, hist.Hist):
        hists[key] = val.group(process, process_cats, process_map)

# this should be custom
name_out = "T_n2_1504_hbb_incl_cat_pre.root"
name_out_pass = name_out.replace('cat','pass')
name_out_fail = name_out.replace('cat','fail')

if os.path.exists(name_out_pass):
    os.remove(name_out_pass)
fout_pass = uproot.create(name_out_pass)
if os.path.exists(name_out_fail):
    os.remove(name_out_fail)
fout_fail = uproot.create(name_out_fail)

h = hists['templates_Wtagregion']
lumi = 41.5
nodata = re.compile("(?!data_obs)")
h.scale({p: lumi for p in h[nodata].identifiers('process')}, axis="process")

rename = {
    'trigweight': 'trigger',
    'pileupweight': 'pu',
    'mutrigweight': 'mutrigger',
    'muidweight': 'muid',
    'muisoweight': 'muiso',
    'matchedUp': 'matched',
    'matchedDown': 'unmatched',
}

for proc in h.identifiers('process'):
    for syst in h.identifiers('systematic'):
        mproj = (slice(None), 'all')
        systreal = syst
        fail_template = (h.project('process', proc)
                          .project('systematic', systreal)
                          .project('AK8Puppijet0_pt', overflow='all')
                          .project('ak8jet_n2ddt', slice(None,0.), overflow='under')[50:100]
                        )
        pass_template = (h.project('process', proc)
                            .project('systematic', systreal)
                            .project('AK8Puppijet0_pt', overflow='all')
                            .project('ak8jet_n2ddt', slice(0.,None))[50:100]
                         )
        content = fail_template.sum('AK8Puppijet0_msd').values()
        if content == {} or content[()] == 0.:
            print(proc, syst)
            continue
        sname = "_%s" % syst if syst != '' else ''
        for k,v in rename.items():
            sname = sname.replace(k, v)
        #name = "%s_pass%s" % (proc, sname)
        name = "%s%s" % (proc, sname) 
        fout_pass[name] = export.export1d(pass_template)
        #name = "%s_fail%s" % (proc, sname)
        fout_fail[name] = export.export1d(fail_template)

fout_pass.close()
fout_fail.close()

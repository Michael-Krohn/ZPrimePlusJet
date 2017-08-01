import ROOT
from ROOT import TFile, TTree, TChain, gPad, gDirectory
from multiprocessing import Process
from optparse import OptionParser
from operator import add
import math
import sys
import time
import array

from plotHelpers import *
from sampleContainer_Heppy import *
DBTAGMIN=-99
#

##############################################################################
def main(options,args):
    #idir = "/eos/uscms/store/user/lpchbb/ggHsample_V11/sklim-v0-28Oct/"
    #odir = "plots_2016_10_31/"
    idir = 'root://cmseos.fnal.gov//store/user/mkrohn/DASZLE/Heppy/v2'
#    idir = 'root://cmseos.fnal.gov//store/user/lpchbb/zprimebits-v12.04/norm2/cvernier'
#    idir = options.idir
    idirData = 'root://cmseos.fnal.gov//store/user/mkrohn/DASZLE/Heppy/v2'
#    idirData = 'root://cmseos.fnal.gov//store/user/lpchbb/zprimebits-v12.04/norm2/cvernier'
#    idirData = 'root://cmseos.fnal.gov//eos/uscms/store/user/lpchbb/zprimebits-v12.05/'
    odir = options.odir
    lumi = options.lumi
    
    legname = {'ggH': 'GF H(b#bar{b})',
	       'VBF': 'VBF H(b#bar{b})',
               'VH':'VH(b#bar{b})',
	       'ttH': 't#bar{t}H(b#bar{b})',
	       'QCD': 'QCD',
	       'ZZ': 'Z',
	       'WW': 'W',	   	
	       'Phibb': ' Phi(125)(b#bar{b})'}

        
    tfiles = {'ggH': [idirData + '/GluGluHToBB_0.root'],
	       'VBF': [idir+'/VBFHToBB_0.root'],
               'VH': [idir+'/ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8_ext_1000pb_weighted.root',
			idir+'/ZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8_1000pb_weighted.root',
			idir+'/WminusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8_1000pb_weighted.root',
			idir+'/WplusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8_1000pb_weighted.root',
			idir+'/ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8_1000pb_weighted.root',
			idir+'/ggZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8_1000pb_weighted.root'],
               'ttH': [idir+'/ttHTobb_M125_13TeV_powheg_pythia8_1000pb_weighted.root'],
               'DYqq' : [idir+'/DYJetsToQQ_HT180_13TeV_0_1000pb_weighted.root'],
#                        idir+'/WplusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8_1000pb_weighted.root'],
#               'tthbb' : [idir+'/ttHTobb_M125_13TeV_powheg_pythia8_1000pb_weighted.root'],
	       'ZZ': [idir+'/ZZ_13TeV_pythia8_1000pb_weighted.root'],
#			idir+'/WWTo4Q_13TeV_powheg_0_1000pb_weighted.root'],
	       'WW':[idir+'/WW_13TeV_pythia8_1000pb_weighted.root']
               }

    color = {'ggH': ROOT.kPink+5,
	     'VBF': ROOT.kBlue+2,
             'VH': ROOT.kAzure+6,
	     'ttH': ROOT.kOrange+1,
	     'QCD': ROOT.kBlack,
	     'ZZ': ROOT.kAzure+1,	
	     'WW': ROOT.kOrange+1,
	     'Phibb':ROOT.kRed-2
               }

    style = {'ggH': 1,
	     'VBF': 2,
             'VH': 2,
	     'ttH': 2,
	     'WHbb':1,
             'QCD': 1,
	     'ZZ':1,		
	     'WW':1
               }
        
    print "Signals... "
    sigSamples = {}
    sigSamples['ggH']  = sampleContainer_Heppy('ggH',tfiles['ggH']  , 1, DBTAGMIN,lumi ) 
#    sigSamples['Zqq']  = sampleContainer('Zqq',tfiles['Zqq']  , 1, DBTAGMIN,lumi)
#    sigSamples['Wqq'] = sampleContainer('Wqq',tfiles['Wqq'], 1, DBTAGMIN,lumi) 
#    sigSamples['DYqq'] = sampleContainer('DYqq',tfiles['DYqq'], 1, DBTAGMIN,lumi ) 	
    sigSamples['VBF'] = sampleContainer_Heppy('VBF',tfiles['VBF'], 1, DBTAGMIN,lumi )
#    sigSamples['VH'] = sampleContainer('VH',tfiles['VH'], 1, DBTAGMIN,lumi )	
#    sigSamples['ttH'] = sampleContainer('ttH',tfiles['ttH'], 1, DBTAGMIN,lumi )
    #sigSamples['Phibb'] = sampleContainer('Phibb',tfiles['Phibb'], 1, lumi*0.035)   


    ofile = ROOT.TFile.Open(odir+'/Plots_1000pb_weighted.root','recreate')



    plots = [
#'h_n_ak4'           ,
'h_met_dbtagCut'     ,
'h_DeltaR_0'     ,
'h_DeltaR_1'     ,
#'h_DeltaR_Higgs_AK8' ,
#'h_pt_ak8'          ,
'h_pt_ak8_dbtagCut' ,
#'h_msd_ak8'         ,
'h_msd_ak8_dbtagCut',
'h_msd_ak8_dbtagCut_Cuts',
'h_dEta_ak4_dbtagCut',
'h_dEta_ak4',
'h_dEta_ak4_nonH_dbtagCut',
'h_dEta_ak4_nonH_dbtagCut_pT30',
'h_dEta_ak4_nonH_dbtagCut_pT50',
'h_dEta_ak4_nonH_dbtagCut_pT70',
'h_dEta_ak4_nonH_dbtagCut_pT100',
'h_mjj_ak4'          ,
'h_dEta_ak4_nonH'    ,
'h_msd_ak8_4ak4'     ,
#'h_mjj_ak4_qq'       ,
'h_mjj_ak4_qq_pT30'       ,
'h_mjj_ak4_qq_pT50'       ,
'h_mjj_ak4_qq_pT70'       ,
'h_mjj_ak4_qq_pT100'       ,
'h_ak4_multiplicity_pT30' ,
'h_ak4_multiplicity_pT50' ,
'h_ak4_multiplicity_pT70' ,
'h_ak4_multiplicity_pT100' ,
#'h_msd_ak8_t21ddtCut'  ,
#'h_msd_ak8_N2Cut'   ,
#'h_dbtag_ak8'       ,
#'h_t21ddt_ak8'      ,
#'h_t32_ak8'         ,
#'h_t32_ak8_t21ddtCut'  ,
#'h_n2b1sd_ak8'      ,
#'h_n2b1sdddt_ak8'   ,
#'h_pt_bbleading' ,
#'h_bb_bbleading' ,
#'h_msd_bbleading',
#'h_msd_ak8_inc',
#'h_msd_ak8_raw',
#'h_msd_ak8_topR6_N2_pass',
#'h_msd_ak8_topR6_N2_fail'
]
#'h_pt_ak8','h_msd_ak8','h_dbtag_ak8','h_n_ak4','h_n_ak4_dR0p8','h_pt_ak8_dbtagCut','h_msd_ak8_dbtagCut','h_t21_ak8','h_t32_ak8','h_msd_ak8_t21ddtCut','h_msd_ak8_N2Cut','h_met']
    for plot in plots:
        hs = {}
        for process, s in sigSamples.iteritems():
            hs[process] = getattr(s,plot)
        c = makeCanvasComparison(hs,legname,color,style,plot.replace('h_','signalcomparison_'),odir,lumi)

        ofile.cd()
        for process, h in hs.iteritems():
            h.Write()

        
        c.Write()
	


##----##----##----##----##----##----##
if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
    parser.add_option("--lumi", dest="lumi", type=float,default = 35.9,help="luminosity", metavar="lumi")
    parser.add_option('-i','--idir', dest='idir', default = 'data/',help='directory with data', metavar='idir')
    parser.add_option('-o','--odir', dest='odir', default = 'plots/',help='directory to write plots', metavar='odir')

    (options, args) = parser.parse_args()

     
    import tdrstyle
    tdrstyle.setTDRStyle()
    ROOT.gStyle.SetPadTopMargin(0.10)
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadRightMargin(0.10)
    ROOT.gStyle.SetPalette(1)
    ROOT.gStyle.SetPaintTextFormat("1.1f")
    ROOT.gStyle.SetOptFit(0000)
    ROOT.gROOT.SetBatch()
    
    main(options,args)
##----##----##----##----##----##----##





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
from sampleContainer import *
DBTAGMIN=-99
#

##############################################################################
def main(options,args):
    #idir = "/eos/uscms/store/user/lpchbb/ggHsample_V11/sklim-v0-28Oct/"
    #odir = "plots_2016_10_31/"
    idir = options.idir
    idirData = 'root://cmseos.fnal.gov//store/user/lpchbb/zprimebits-v12.04/cvernier'
#    idirData = 'root://cmseos.fnal.gov//eos/uscms/store/user/lpchbb/zprimebits-v12.05/'
    odir = options.odir
    lumi = options.lumi
    
    legname = {'ggH': 'H',
	       'Zqq': 'Z(q#bar{q})+Jets',
               'Wqq':'W(q#bar{q})+Jets',
	       'DYqq': ' Drell-Yan(q#bar{q})+Jets',
	       'QCD': 'QCD',
	       'ZZ': 'Z',
	       'WW': 'W',	   	
	       'Phibb': ' Phi(125)(b#bar{b})'}

        
    tfiles = {'ggH': [idirData + '/GluGluHToBB_M125_13TeV_powheg_pythia8_all_1000pb_weighted.root'],
	       'Zqq': [idir+'/ZJetsToQQ_HT_600ToInf_13TeV_1000pb_weighted.root'],
               'Wqq': [idir+'/WJetsToQQ_HT180_13TeV_0_1000pb_weighted.root'],
               'QCD': [idir+'/QCD_HT300to500_13TeV_1000pb_weighted.root',
			idir+'/QCD_HT500to700_13TeV_1000pb_weighted.root',
			idir+'/QCD_HT700to1000_13TeV_1000pb_weighted.root',
			idir+'/QCD_HT1500to2000_13TeV_1000pb_weighted.root',
			idir+'/QCD_HT2000toInf_13TeV_1000pb_weighted.root'],
               'DYqq' : [idir+'/DYJetsToQQ_HT180_13TeV_0_1000pb_weighted.root'],
#                        idir+'/WplusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8_1000pb_weighted.root'],
#               'tthbb' : [idir+'/ttHTobb_M125_13TeV_powheg_pythia8_1000pb_weighted.root'],
	       'ZZ': [idir+'/ZZ_13TeV_pythia8_1000pb_weighted.root'],
#			idir+'/WWTo4Q_13TeV_powheg_0_1000pb_weighted.root'],
	       'WW':[idir+'/WW_13TeV_pythia8_1000pb_weighted.root']
               }

    color = {'ggH': ROOT.kPink+5,
	     'Zqq': ROOT.kBlue+2,
             'Wqq': ROOT.kAzure+3,
	     'DYqq': ROOT.kPink+5,
	     'QCD': ROOT.kBlack,
	     'ZZ': ROOT.kAzure+1,	
	     'WW': ROOT.kOrange+1,
	     'Phibb':ROOT.kRed-2
               }

    style = {'ggH': 1,
	     'Zqq': 1,
             'Wqq': 1,
	     'DYqq': 1,
	     'WHbb':1,
             'QCD': 1,
	     'ZZ':1,		
	     'WW':1
               }
        
    print "Signals... "
    sigSamples = {}
    sigSamples['ggH']  = sampleContainer('ggH',tfiles['ggH']  , 1, DBTAGMIN,lumi ) 
#    sigSamples['Zqq']  = sampleContainer('Zqq',tfiles['Zqq']  , 1, DBTAGMIN,lumi)
#    sigSamples['Wqq'] = sampleContainer('Wqq',tfiles['Wqq'], 1, DBTAGMIN,lumi) 
#    sigSamples['DYqq'] = sampleContainer('DYqq',tfiles['DYqq'], 1, DBTAGMIN,lumi ) 	
    sigSamples['QCD'] = sampleContainer('QCD',tfiles['QCD'], 1, DBTAGMIN,lumi )
    sigSamples['ZZ'] = sampleContainer('ZZ',tfiles['ZZ'], 1, DBTAGMIN,lumi )	
    sigSamples['WW'] = sampleContainer('WW',tfiles['WW'], 1, DBTAGMIN,lumi )
    #sigSamples['Phibb'] = sampleContainer('Phibb',tfiles['Phibb'], 1, lumi*0.035)   


    ofile = ROOT.TFile.Open(odir+'/Plots_1000pb_weighted.root','recreate')



    plots = [
#'h_n_ak4'           ,
#'h_met'             ,
#'h_pt_ak8'          ,
#'h_pt_ak8_dbtagCut' ,
'h_msd_ak8'         ,
#'h_msd_ak8_dbtagCut',
#'h_msd_ak8_t21ddtCut'  ,
#'h_msd_ak8_N2Cut'   ,
'h_dbtag_ak8'       ,
'h_t21ddt_ak8'      ,
#'h_t32_ak8'         ,
#'h_t32_ak8_t21ddtCut'  ,
'h_n2b1sd_ak8'      ,
'h_n2b1sdddt_ak8'   ,
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





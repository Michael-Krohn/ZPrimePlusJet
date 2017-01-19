#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
from optparse import OptionParser

import ROOT

# ROOT.gROOT.ProcessLine(".L ~/tdrstyle.C")
# ROOT.setTDRStyle()
# ROOT.gStyle.SetPadTopMargin(0.06)
# ROOT.gStyle.SetPadLeftMargin(0.16)
# ROOT.gStyle.SetPadRightMargin(0.10)
# ROOT.gStyle.SetPalette(1)
# ROOT.gStyle.SetPaintTextFormat("1.1f")


############################################################

# observableTraining takes in 2 root files, list of observables, spectator observables ... launches a CONDOR job
# TMVAhelper.py is used by observableTraining
# analysis.py defines the list of trainings, the signal and background process

########################################################################################################################
########################################################################################################################

def main(options,args):

    DataDir = options.idir
    OutDir = options.odir

    tags = []
    #tags.append([ 'JetHTRun2016B_23Sep2016_v3/', 0] )
    #tags.append([ 'JetHTRun2016B_23Sep2016_v1/', 0] )
    #tags.append([ 'JetHTRun2016C_23Sep2016_v1/', 0] )
    tags.append([ 'JetHTRun2016D_23Sep2016_v1/', 0] )
    #tags.append([ 'JetHTRun2016E_23Sep2016_v1/', 0] )
    tags.append([ 'JetHTRun2016F_23Sep2016_v1/', 0] )
    #tags.append([ 'JetHTRun2016G_23Sep2016_v1/', 0] )
    tags.append([ 'JetHTRun2016H_PromptReco_v1/', 0] )
    tags.append([ 'JetHTRun2016H_PromptReco_v2/', 0] )
    #tags.append([ 'JetHTRun2016H_PromptReco_v3/', 0] )
    tags.append([ 'SingleMuonRun2016B_23Sep2016_v1/', 0] )
    tags.append([ 'SingleMuonRun2016B_23Sep2016_v3/', 0] )
    tags.append([ 'SingleMuonRun2016C_23Sep2016_v1/', 0] )
    tags.append([ 'SingleMuonRun2016D_23Sep2016_v1/', 0] )
    tags.append([ 'SingleMuonRun2016E_23Sep2016_v1/', 0] )
    tags.append([ 'SingleMuonRun2016F_23Sep2016_v1/', 0] )
    tags.append([ 'SingleMuonRun2016G_23Sep2016_v1/', 0] )
    tags.append([ 'SingleMuonRun2016H_PromptReco_v1/', 0] )
    #tags.append([ 'SingleMuonRun2016H_PromptReco_v2/', 0] )
    #tags.append([ 'SingleMuonRun2016H_PromptReco_v3/', 0] )

    postfix = ''
    for i in range(len(tags)):
        basename = tags[i][0].replace('/','') + '.root'
        filesToConvert, badFiles = getFilesRecursively(DataDir,tags[i][0],None,'sklim')
        print "files To Convert = ",filesToConvert
        print "bad files = ", badFiles
        haddCommand = 'hadd -f %s/%s %s'%(OutDir,basename,' '.join(filesToConvert))
        print haddCommand
        with open('hadd_command_%s.sh'%basename,'w') as f:
            f.write(haddCommand)
        with open('bad_files_%s.txt'%basename,'w') as f:
            for badFile in badFiles:
                f.write(badFile+'\n')

        os.system('source hadd_command_%s.sh'%basename)




def getFilesRecursively(dir,searchstring,additionalstring = None, skipString = None):
    
    # thesearchstring = "_"+searchstring+"_"
    thesearchstring = searchstring

    theadditionalstring = None
    if not additionalstring == None: 
        theadditionalstring = additionalstring

    cfiles = []
    badfiles = []
    for root, dirs, files in os.walk(dir+'/'+thesearchstring):
        nfiles = len(files)
        for ifile, file in enumerate(files):
            
            if ifile%100==0:
                print '%i/%i files checked in %s'%(ifile,nfiles,dir+'/'+thesearchstring)
            try:
                f = ROOT.TFile.Open(os.path.join(root, file))
                if f.IsZombie():
                    print 'file is zombie'
                    f.Close()
                    badfiles.append(os.path.join(root, file))                    
                    continue
                elif not f.Get('Events'):
                    print 'tree is false'
                    f.Close()
                    badfiles.append(os.path.join(root, file))                    
                    continue
                elif not f.Get('Events').InheritsFrom('TTree'):
                    print 'tree is not a tree'
                    f.Close()
                    badfiles.append(os.path.join(root, file))                    
                    continue
                else:
                    f.Close()
                    cfiles.append(os.path.join(root, file))                    
            except:
                print 'could not open file or tree'
                badfiles.append(os.path.join(root, file))                    
                continue
                
    return cfiles, badfiles



if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
    parser.add_option('--train', action='store_true', dest='train', default=False, help='train')
    parser.add_option("--lumi", dest="lumi", default = 30,type=float,help="luminosity", metavar="lumi")
    parser.add_option('-i','--idir', dest='idir', default = 'data/',help='directory with bacon bits', metavar='idir')
    parser.add_option('-o','--odir', dest='odir', default = 'skim/',help='directory to write hadded bits', metavar='odir')

    (options, args) = parser.parse_args()

    main(options,args)

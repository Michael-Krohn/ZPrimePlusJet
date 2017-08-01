import ROOT as rt
from multiprocessing import Process
from optparse import OptionParser
from operator import add
import math
import sys
import time
import array
import warnings
import os


class sampleContainer:

    def __init__( self, name, tt, tf, cutFormula='1'):
        self._name = name
        self._tf = tf
        self._tt = tt
        warnings.filterwarnings( action='ignore', category=RuntimeWarning, message='creating converter.*' )
        self._cutFormula = rt.TTreeFormula("cutFormula",cutFormula,self._tt)

        self._branches = [('l1Jets_pt','f',-999),('l1Jets_eta','f',-999),('l1Jets_phi','f',-999),('l1MET_pt','f',-999),('l1MET_phi','f',-999),
                ('Jet_pt','f',-999),('Jet_eta','f',-999),('Jet_phi','f',-999),('Jet_mass','f',-999),
                ('HCSV_reg_mass','f',-999),('met_pt','f',-999),('met_phi','f',-999), ('nl1Jets', 'i', -999),('nJet', 'i', -999),
                ('Vtype','i',-999),('tkMet_phi','f',-999),('Jet_puId','i',-999),('Jet_id','i',-999),('hJCidx','i',-999),('nselLeptons','i',-999)]

        
        self._tt.SetBranchStatus("*",0)
        for branch in self._branches:
            self._tt.SetBranchStatus(branch[0],1)
        for branch in self._branches:
            setattr(self, branch[0].replace(' ', ''), array.array(branch[1],[branch[2]]*50))
            self._tt.SetBranchAddress( branch[0], getattr(self, branch[0].replace(' ', '')) )

def convertPhi(iphi):
    return float(iphi)* 2. * rt.TMath.Pi() / 144. - rt.TMath.Pi()

def convertEta(ieta):
    return float(ieta) * 2. * rt.TMath.Pi() / 144.
    
def setStyle():
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetOptFit(0000)
    rt.gStyle.SetOptTitle(0)
    rt.gStyle.SetPaintTextFormat("1.2g")
    #rt.gStyle.SetPalette()
    rt.gStyle.SetNumberContours(999)
    rt.gROOT.SetBatch()
    rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.FATAL)    
    rt.gStyle.SetStatY(1.9)
    rt.gStyle.SetStatX(1.9)
    rt.gStyle.SetStatW(0.1)
    rt.gStyle.SetStatH(0.1)    

def main(options, args):
    tree = rt.TChain("tree")
    count = 0
    tfiles = []
    for arg in args:
        tree.Add(arg)
        tfile = rt.TFile.Open(arg)
        tfiles.append(tfile)
        count += tfile.Get("Count").Integral()


    
    VtypeZnn = 'Vtype==4'
    AntiQCD = 'Sum$(abs(TVector2::Phi_mpi_pi(Jet_phi-met_phi))<0.5 && Jet_pt>30 && Jet_puId>=4)==0'
    Jets = 'Jet_pt[hJCidx[0]]>60 && Jet_pt[hJCidx[1]]>35 && Min$(Jet_id[hJCidx])>=4 && Min$(Jet_puId[hJCidx])>=4 && Jet_id[0]>=4 && Jet_puId[0]>=4'
    NoLeptons = 'nselLeptons==0'
    NotFourJets  = 'Sum$(Jet_pt>30 && abs(Jet_eta)<5.2 && Jet_puId>=4)<4'
    MetTkMet = 'abs(TVector2::Phi_mpi_pi(met_phi-tkMet_phi))>2.64'

    if options.presel:
        cuts = '(nl1Jets >= 2) && (%s) && (%s) && (%s) && (%s) && (%s) && (%s)'%(VtypeZnn, AntiQCD, Jets, MetTkMet, NoLeptons, NotFourJets)
    else:
        #cuts = '(nl1Jets >= 2)'
        cuts = '1'
        
    sample = sampleContainer('ZeroBias', tree, tfile, cuts)
    print 'cuts = %s'%cuts
        
    #CSVM         = 'Jet_btagCSV[hJCidx[0]]>0.800'
    #NotCSVM      = 'Jet_btagCSV[hJCidx[0]]<0.800'

    l1Cuts = {'L1_ETM120': 'sample.l1MET_pt[0] > 120',
              'L1_ETM100': 'sample.l1MET_pt[0] > 100',
              'L1_ETM90': 'sample.l1MET_pt[0] > 90',
              'L1_ETM75': 'sample.l1MET_pt[0] > 75',
              'L1_Jet60': 'sample.l1Jets_pt[0] > 60',
              'L1_JetC60': 'sample.l1Jets_pt[0] > 60',
              'L1_DoubleJet60': 'sample.l1Jets_pt[0] > 60 and sample.l1Jets_pt[1] > 60',
              'L1_DoubleJet60_dR_Min1p5': 'sample.l1Jets_pt[0] > 60 and sample.l1Jets_pt[1] > 60 and dRMin > 1.5',
              'L1_DoubleJet60_dR_Max1p5': 'sample.l1Jets_pt[0] > 60 and sample.l1Jets_pt[1] > 60 and dRMax < 1.5',
              'L1_DoubleJetC60_dR_Min1p5': 'sample.l1Jets_pt[0] > 60 and abs(convertEta(sample.l1Jets_eta[0]))<3 and sample.l1Jets_pt[1] > 60 and abs(convertEta(sample.l1Jets_eta[1]))<3 and dRMinC > 1.5',
              'L1_DoubleJetC60_dR_Max1p5': 'sample.l1Jets_pt[0] > 60 and abs(convertEta(sample.l1Jets_eta[0]))<3 and sample.l1Jets_pt[1] > 60 and abs(convertEta(sample.l1Jets_eta[1]))<3 and dRMaxC < 1.5',
              'L1_DoubleJetC60': 'sample.l1Jets_pt[0] > 60 and abs(convertEta(sample.l1Jets_eta[0]))<3 and sample.l1Jets_pt[1] > 60 and abs(convertEta(sample.l1Jets_eta[1]))<3',
              'L1_DoubleJet_60_30': 'sample.l1Jets_pt[0] > 60 and sample.l1Jets_pt[1] > 30',
              'L1_DoubleJetC_60_30': 'sample.l1Jets_pt[0] > 60 and abs(convertEta(sample.l1Jets_eta[0]))<3 and sample.l1Jets_pt[1] > 30 and abs(convertEta(sample.l1Jets_eta[1]))<3',              
              'L1_DoubleJet60_30_Mj60j30_80': 'sample.l1Jets_pt[0] > 60 and sample.l1Jets_pt[1] > 30 and l1JetMass > 80',
              'L1_ETM70_DoubleJet30_dPhi_MinLead20p5': 'sample.l1Jets_pt[0] > 30 and sample.l1Jets_pt[1] > 30 and sample.l1MET_pt[0] > 70 and dPhiMinLead2 > 0.5',
              'L1_ETM75_Jet60': 'sample.l1Jets_pt[0] > 60 and sample.l1MET_pt[0] > 75',
              'L1_ETM75_Jet60_dPhi_Min0p4_Mj60j30_80': 'sample.l1Jets_pt[0] > 60 and sample.l1Jets_pt[1] > 30 and sample.l1MET_pt[0] > 75 and dPhiMin > 0.4 and l1JetMass > 80',
              'L1_ETM75_Jet60_dPhi_Min0p4': 'sample.l1Jets_pt[0] > 60 and sample.l1MET_pt[0] > 75 and dPhiMin > 0.4',
              'L1_ETM75_DoubleJet30_dPhi_MinLead20p4': 'sample.l1Jets_pt[0] > 30 and sample.l1Jets_pt[1] > 30 and sample.l1MET_pt[0] > 75 and dPhiMinLead2 > 0.4',
              'L1_ETM75_Jet60_dPhi_Min0p6': 'sample.l1Jets_pt[0] > 60 and sample.l1MET_pt[0] > 75 and dPhiMin > 0.6',
              'L1_ETM75_DoubleJet_60_30': 'sample.l1Jets_pt[0] > 60 and sample.l1Jets_pt[1] > 30 and sample.l1MET_pt[0] > 75',
              'L1_ETM85_Jet60': 'sample.l1Jets_pt[0] > 60 and sample.l1MET_pt[0] > 85',
              'L1_ETM85_Jet60_dPhi_Min0p4_Mj60j30_80': 'sample.l1Jets_pt[0] > 60 and sample.l1Jets_pt[1] > 30 and sample.l1MET_pt[0] > 85 and dPhiMin > 0.4 and l1JetMass > 80',
              'L1_ETM85_Jet60_dPhi_Min0p4': 'sample.l1Jets_pt[0] > 60 and sample.l1MET_pt[0] > 85 and dPhiMin > 0.4',
              'L1_ETM85_Jet60_dPhi_Min0p6': 'sample.l1Jets_pt[0] > 60 and sample.l1MET_pt[0] > 85 and dPhiMin > 0.6',
              'L1_ETM85_DoubleJet_60_30': 'sample.l1Jets_pt[0] > 60 and sample.l1Jets_pt[1] > 30 and sample.l1MET_pt[0] > 85',
              'L1_ETM85_DoubleJet30_dPhi_MinLead20p4': 'sample.l1Jets_pt[0] > 30 and sample.l1Jets_pt[1] > 30 and sample.l1MET_pt[0] > 85 and dPhiMinLead2 > 0.4',
              'L1_ETM85_DoubleJet30_dPhi_MinLead20p4_Mjj_60': 'sample.l1Jets_pt[0] > 30 and sample.l1Jets_pt[1] > 30 and sample.l1MET_pt[0] > 85 and dPhiMinLead2 > 0.4 and l1JetMass > 60',
              'L1_ETM90_DoubleJet30_dPhi_MinLead20p4_Mjj_60': 'sample.l1Jets_pt[0] > 30 and sample.l1Jets_pt[1] > 30 and sample.l1MET_pt[0] > 90 and dPhiMinLead2 > 0.4 and l1JetMass > 60',
              'L1_ETM90_DoubleJet30_dPhi_MinLead20p4': 'sample.l1Jets_pt[0] > 30 and sample.l1Jets_pt[1] > 30 and sample.l1MET_pt[0] > 90 and dPhiMinLead2 > 0.4',
              'L1_ETM95_DoubleJet30_dPhi_MinLead20p4': 'sample.l1Jets_pt[0] > 30 and sample.l1Jets_pt[1] > 30 and sample.l1MET_pt[0] > 95 and dPhiMinLead2 > 0.4',
              'L1_ETM95_DoubleJet30_dPhi_METLead2Jets0p4': 'sample.l1Jets_pt[0] > 30 and sample.l1Jets_pt[1] > 30 and sample.l1MET_pt[0] > 95 and dPhiMETLead2Jets > 0.4',
              'L1_ETM95_DoubleJet30_dPhi_METLead2Jets0p2': 'sample.l1Jets_pt[0] > 30 and sample.l1Jets_pt[1] > 30 and sample.l1MET_pt[0] > 95 and dPhiMETLead2Jets > 0.2',
              'L1_ETM95_Jet60': 'sample.l1Jets_pt[0] > 60 and sample.l1MET_pt[0] > 95',
              'L1_ETM95_Jet60_dPhi_Min0p4_Mj60j30_80': 'sample.l1Jets_pt[0] > 60 and sample.l1Jets_pt[1] > 30 and sample.l1MET_pt[0] > 95 and dPhiMin > 0.4 and l1JetMass > 80',
              'L1_ETM95_Jet60_dPhi_Min0p4': 'sample.l1Jets_pt[0] > 60 and sample.l1MET_pt[0] > 95 and dPhiMin > 0.4',
              'L1_ETM95_Jet60_dPhi_Min0p6': 'sample.l1Jets_pt[0] > 60 and sample.l1MET_pt[0] > 95 and dPhiMin > 0.6',
              'L1_ETM95_DoubleJet_60_30': 'sample.l1Jets_pt[0] > 60 and sample.l1Jets_pt[1] > 30 and sample.l1MET_pt[0] > 95',
              'L1_DoubleJetC60_ETM60': 'sample.l1Jets_pt[0] > 60 and abs(convertEta(sample.l1Jets_eta[0]))<3 and sample.l1Jets_pt[1] > 60 and abs(convertEta(sample.l1Jets_eta[1]))<3 and sample.l1MET_pt[0] > 60',
              'L1_DoubleJetC60_ETM60': 'sample.l1Jets_pt[0] > 60 and abs(convertEta(sample.l1Jets_eta[0]))<3 and sample.l1Jets_pt[1] > 60 and abs(convertEta(sample.l1Jets_eta[1]))<3 and sample.l1MET_pt[0] > 60',
              }
    
    if '||' in options.l1CutName:
        l1CutNames = options.l1CutName.split("||")
        l1OrCuts = ['('+l1Cuts[myCut]+')' for myCut in l1CutNames]
        l1Cut = ' or '.join(l1OrCuts)        
    else:
        l1Cut = l1Cuts[options.l1CutName]

    CSVM         = 'Jet_btagCSV[hJCidx[0]]>0.800'
    Hmass        = '(HCSV_reg_mass>100&&HCSV_reg_mass<140)'
    FourJets     = 'Sum$(Jet_pt>30 && abs(Jet_eta)<5.2 && Jet_puId>=4)>=4'
    NotFourJets  = 'Sum$(Jet_pt>30 && abs(Jet_eta)<5.2 && Jet_puId>=4)<4'
    AntiQCD      = 'Sum$(abs(TVector2::Phi_mpi_pi(Jet_phi-met_phi))<0.5 && Jet_pt>30 && Jet_puId>=4)==0'
    Signal = '<!Cuts|Jets!> && HCSV_mass<500 && <!Cuts|VtypeZnn!> && <!Cuts|CSVM!> && nselLeptons==0 && <!Cuts|AntiQCD!> && <!Cuts|NotFourJets!> && abs(TVector2::Phi_mpi_pi(met_phi-tkMet_phi))>2.64'
    SignalTight = '<!Cuts|Jets!> && HCSV_mass<500 && <!Cuts|VtypeZnn!> && <!Cuts|Hmass!> && <!Cuts|CSVM!> && nselLeptons==0 && <!Cuts|AntiQCD!> && <!Cuts|NotFourJets!> && abs(TVector2::Phi_mpi_pi(met_phi-tkMet_phi))>2.64'
    
    histDef = {'l1Jets_pt0':[100,0,1000],
               'l1Jets_phi0':[100,-3.142,3.142],
               'l1Jets_eta0':[100,-5,5],
               'l1Jets_pt1':[100,0,1000],
               'l1Jets_phi1':[100,-3.142,3.142],
               'l1Jets_eta1':[100,-5,5],
               'l1Jets_mass':[100,0,1000],
               'l1Jets_dphi':[100,0,3.142],
               'l1Jets_dr':[100,0,6],
               'l1Jets_drMin':[100,0,6],
               'l1Jets_drMax':[100,0,6],
               'l1Jets_drMinC':[100,0,6],
               'l1Jets_drMaxC':[100,0,6],
               'l1MET_pt':[100,0,600],
               'l1MET_phi':[100,-3.142,3.142],
               'l1MET_Jet_dphiMin':[100,0,3.142],
               'l1MET_Jet_dphiMinLead2':[100,0,3.142],
               'l1MET_Jet_dphiMETLead2Jets':[100,0,3.142],
               'Jet_pt0':[100,60,1000],
               'Jet_pt1':[100,30,1000],
               'Jet_eta0':[100,-5,5],
               'Jet_eta1':[100,-5,5],
               'Jet_phi0':[100,-3.142,3.142],
               'Jet_phi1':[100,-3.142,3.142],
               'Jet_mass':[100,0,1000],
               'met_pt':[100,80,600],
               'HCSV_reg_mass':[100,0,1000]                         
               }

    effDef = {'effJet_pt0': ['effJet_pt0; Leading jet p_{T} [GeV];efficiency', 50, 20, 500],
              'effJet_mass': ['effJet_mass; Leading jets mass [GeV];efficiency', 50, 20, 500],
              'effHCSV_reg_mass': ['effHCSV_reg_mass; Higgs candidate reg. mass [GeV];efficiency', 50, 20, 500],
              'effmet_pt': ['effmet_pt; E_{T}^{miss} [GeV];efficiency', 50, 80, 500]
              }
    

        
    hist = {}
    eff = {}
               
    for histName, histSpec in histDef.iteritems():
        hist[histName] = rt.TH1D(histName,histName,histSpec[0],histSpec[1],histSpec[2])
        
    for effName, effSpec in effDef.iteritems():
        eff[effName] = rt.TEfficiency(effName,effSpec[0],effSpec[1],effSpec[2],effSpec[3])        

        
    #sample._tt.Draw('>>elist',cuts,'entrylist')
    #elist = rt.gDirectory.Get('elist')
    #entry = -1

    presel_events = 0
    l1_events = 0
    l1etm120_events = 0
    l1_exclusive_events = 0
    l1_or_events = 0
    
    nent = sample._tt.GetEntries()
    print nent
    sample._tt.SetNotify(sample._cutFormula)
    for entry in xrange(nent):
        sample._tt.LoadTree(entry)
        selected = False
        for j in range(sample._cutFormula.GetNdata()):
            if (sample._cutFormula.EvalInstance(j)):
                selected = True
                break
        if not selected: continue

        sample._tt.GetEntry(entry)
        presel_events += 1.0           
        if presel_events>1000000: break
        if presel_events%10000==0:
            print "%i preselected events processsed"%presel_events
        if sample.nl1Jets[0]>0:
            hist['l1Jets_pt0'].Fill(sample.l1Jets_pt[0])
            hist['l1Jets_eta0'].Fill(convertEta(sample.l1Jets_eta[0]))
            hist['l1Jets_phi0'].Fill(convertPhi(sample.l1Jets_phi[0]))
        if sample.nl1Jets[0]>1:
            hist['l1Jets_pt1'].Fill(sample.l1Jets_pt[1])
            hist['l1Jets_eta1'].Fill(convertEta(sample.l1Jets_eta[1]))
            hist['l1Jets_phi1'].Fill(convertPhi(sample.l1Jets_phi[1]))
        hist['l1MET_pt'].Fill(sample.l1MET_pt[0])
        hist['l1MET_phi'].Fill(sample.l1MET_phi[0])

        l1MET = rt.TLorentzVector()
        l1MET.SetPtEtaPhiM(sample.l1MET_pt[0],0.,sample.l1MET_phi[0],0.)

        l1JetMass = -1
        if sample.nl1Jets[0]>1:
            l1Jet0 = rt.TLorentzVector()
            l1Jet1 = rt.TLorentzVector()
            l1Jet0.SetPtEtaPhiM(sample.l1Jets_pt[0],convertEta(sample.l1Jets_eta[0]),convertPhi(sample.l1Jets_phi[0]),0.)
            l1Jet1.SetPtEtaPhiM(sample.l1Jets_pt[1],convertEta(sample.l1Jets_eta[1]),convertPhi(sample.l1Jets_phi[1]),0.)
            l1JetMass = (l1Jet0+l1Jet1).M()
            hist['l1Jets_mass'].Fill(l1JetMass)
            hist['l1Jets_dphi'].Fill(l1Jet0.DeltaPhi(l1Jet1))
            hist['l1Jets_dr'].Fill(l1Jet0.DeltaR(l1Jet1))

        dPhiList = []
        dPhiMin = -1
        for i in range(0,sample.nl1Jets[0]):
            l1Jet = rt.TLorentzVector()
            l1Jet.SetPtEtaPhiM(sample.l1Jets_pt[i],convertEta(sample.l1Jets_eta[i]),convertPhi(sample.l1Jets_phi[i]),0.)
            if l1Jet.Pt()>60:
                dPhiList.append(abs(l1Jet.DeltaPhi(l1MET)))
        if len(dPhiList)>0:
            dPhiMin = min(dPhiList)
            hist['l1MET_Jet_dphiMin'].Fill(dPhiMin)
        dPhiListLead2 = []
        dPhiMinLead2 = -1
        for i in range(0,min(2,sample.nl1Jets[0])):
            l1Jet = rt.TLorentzVector()
            l1Jet.SetPtEtaPhiM(sample.l1Jets_pt[i],convertEta(sample.l1Jets_eta[i]),convertPhi(sample.l1Jets_phi[i]),0.)
            if l1Jet.Pt()>30:
                dPhiListLead2.append(abs(l1Jet.DeltaPhi(l1MET)))
        if len(dPhiListLead2)>1:
            dPhiMinLead2 = min(dPhiListLead2)
            hist['l1MET_Jet_dphiMinLead2'].Fill(dPhiMinLead2)
            
        dPhiMETLead2Jets = -1
        if sample.nl1Jets[0] > 1:
            l1Jet0 = rt.TLorentzVector()
            l1Jet0.SetPtEtaPhiM(sample.l1Jets_pt[0],convertEta(sample.l1Jets_eta[0]),convertPhi(sample.l1Jets_phi[0]),0.)
            l1Jet1 = rt.TLorentzVector()
            l1Jet1.SetPtEtaPhiM(sample.l1Jets_pt[1],convertEta(sample.l1Jets_eta[1]),convertPhi(sample.l1Jets_phi[1]),0.)
            if l1Jet0.Pt()>30 and l1Jet1.Pt()>30:
                dPhiMETLead2Jets = abs((l1Jet0+l1Jet1).DeltaPhi(l1MET))
                hist['l1MET_Jet_dphiMETLead2Jets'].Fill(dPhiMETLead2Jets)

        dRList = []
        dRListC = []
        dRMin = -1
        dRMinC = -1
        dRMax = 9999
        dRMaxC = 9999
        for i in range(0,sample.nl1Jets[0]):
            l1Jeta = rt.TLorentzVector()
            l1Jeta.SetPtEtaPhiM(sample.l1Jets_pt[i],convertEta(sample.l1Jets_eta[i]),convertPhi(sample.l1Jets_phi[i]),0.)
            for j in range(i+1,sample.nl1Jets[0]):
                l1Jetb = rt.TLorentzVector()
                l1Jetb.SetPtEtaPhiM(sample.l1Jets_pt[j],convertEta(sample.l1Jets_eta[j]),convertPhi(sample.l1Jets_phi[j]),0.)
                if l1Jeta.Pt()>60 and l1Jetb.Pt()>60:
                    dRList.append(abs(l1Jeta.DeltaR(l1Jetb)))
                    if abs(l1Jeta.Eta())<3 and abs(l1Jetb.Eta())<3:
                        dRListC.append(abs(l1Jeta.DeltaR(l1Jetb)))
        if len(dRListC)>0:
            dRMinC = min(dRListC)
            dRMaxC = max(dRListC)
            hist['l1Jets_drMinC'].Fill(dRMinC)
            hist['l1Jets_drMaxC'].Fill(dRMaxC)
        if len(dRList)>0:
            dRMin = min(dRList)
            dRMax = max(dRList)
            hist['l1Jets_drMin'].Fill(dRMin)
            hist['l1Jets_drMax'].Fill(dRMax)

        try:
            bNum = eval(l1Cut)
        except IndexError:
            bNum = False

        try:
            bNumETM120 = eval(l1Cuts['L1_ETM120'])
        except IndexError:
            bNumETM120 = False
            
        if bNum:
            l1_events += 1.0

        if bNumETM120:
            l1etm120_events += 1.0
            
        if bNum and not (bNumETM120):
            l1_exclusive_events += 1.0
            
        if bNum or bNumETM120:
            l1_or_events += 1.0
            
        jetMass = -1
        if sample.nJet[0]>1:            
            jet0 = rt.TLorentzVector()
            jet1 = rt.TLorentzVector()
            jet0.SetPtEtaPhiM(sample.Jet_pt[0],sample.Jet_eta[0],sample.Jet_phi[0],sample.Jet_mass[0])
            jet1.SetPtEtaPhiM(sample.Jet_pt[1],sample.Jet_eta[1],sample.Jet_phi[1],sample.Jet_mass[1])

            jetMass = (jet0+jet1).M()
            hist['Jet_mass'].Fill(jetMass)
            
        if sample.nJet[0]>0:            
            hist['Jet_pt0'].Fill(sample.Jet_pt[0])
            hist['Jet_eta0'].Fill(sample.Jet_eta[0])
            hist['Jet_phi0'].Fill(sample.Jet_phi[0])
        if sample.nJet[0]>1:    
            hist['Jet_pt1'].Fill(sample.Jet_pt[1])
            hist['Jet_eta1'].Fill(sample.Jet_eta[1])
            hist['Jet_phi1'].Fill(sample.Jet_phi[1])
            
        hist['HCSV_reg_mass'].Fill(sample.HCSV_reg_mass[0])
        hist['met_pt'].Fill(sample.met_pt[0])

            
        if sample.nJet[0]>0:            
            #eff['effJet_pt0'].Fill(sample.l1Jets_pt[0] > 60,sample.Jet_pt[0])
            eff['effJet_pt0'].Fill(bNum,sample.Jet_pt[0])
        eff['effJet_mass'].Fill(bNum,jetMass)
        eff['effHCSV_reg_mass'].Fill(bNum,sample.HCSV_reg_mass[0])
        eff['effmet_pt'].Fill(bNum,sample.met_pt[0])
        
    tfileOut = rt.TFile.Open(options.odir+'/l1Study_%s.root'%options.l1CutName.replace('||','Or'),'recreate')
    for histName, h in hist.iteritems():
        h.Write()
    for effName, pEff in eff.iteritems():
        pEff.Write()
    tfileOut.Close()

    
    setStyle()
    c = rt.TCanvas("c","c",500,400)
    c.SetRightMargin(0.1)

    sigmoids = []
    for effName, pEff in eff.iteritems():
        
        sigmoid = rt.TF1("sigmoid","[0]/(1.0+exp(-(x-[1])/[2]))",
                         pEff.GetTotalHistogram().GetXaxis().GetBinLowEdge(1),
                         pEff.GetTotalHistogram().GetXaxis().GetBinUpEdge(pEff.GetTotalHistogram().GetNbinsX()))
        sigmoid.SetParameter(0,1)
        sigmoid.SetParLimits(0,0,1)
        sigmoid.SetParameter(1,60)
        sigmoid.SetParameter(2,10)
        sigmoids.append(sigmoid)
        
        #pEff.SetLineColor(color)
        pEff.SetMarkerSize(0.8)
        pEff.SetMarkerStyle(20)
        pEff.Draw("apez")
        pEff.Fit(sigmoid,"I")
        rt.gPad.Update()        
        #pEff.GetPaintedHistogram().GetXaxis().SetRangeUser(x[0],x[-1])
        pEff.GetPaintedGraph().SetMarkerStyle(8)
        pEff.GetPaintedGraph().SetMarkerSize(20)        
        pEff.GetPaintedGraph().SetMinimum(0)
        pEff.GetPaintedGraph().SetMaximum(1.3)
        rt.gPad.Update()
        l = rt.TLatex()
        l.SetTextAlign(11)
        l.SetTextSize(0.045)
        l.SetTextFont(42)
        l.SetNDC()
        l.DrawLatex(0.12,0.92,"CMS Simulation")
        l.DrawLatex(0.8,0.92,"13 TeV")
        l.SetTextFont(52)
        l.SetTextSize(0.03)
        l.DrawLatex(0.2,0.8,options.l1CutName.replace('||',' || '))
        l.SetTextSize(0.02)
        l.SetTextFont(42)        
        c.Print(options.odir+"/"+pEff.GetName()+"_"+options.l1CutName.replace('||','Or')+".pdf")
        c.Print(options.odir+"/"+pEff.GetName()+"_"+options.l1CutName.replace('||','Or')+".C")


    #Nbunches = 2592.
    #LHCFreq = 11.2455
    Nbunches = 2800.
    LHCFreq = 11.
    print ''
    print "total events:                                         %i"%(count)
    print ''
    print "total events passing preselection:                    %i"%(presel_events)
    print "total events passing l1 selection & preselection:     %i"%(l1_events)
    print ''
    print "total efficiency passing preselection:                %f"%(float(presel_events/count))
    print "total efficiency passing l1 selection & preselection: %f"%(float(l1_events/count))
    print ''
    print "relative efficiency:                                  %f"%(float(l1_events)/float(presel_events))
    print "rate at 55 PU:                                        %f +/- %f kHz"%(float(l1_events)/float(presel_events)*Nbunches*LHCFreq, getErr(l1_events,presel_events)*Nbunches*LHCFreq)
    print "pure rate at 55 PU:                                   %f +/- %f kHz"%(float(l1_exclusive_events)/float(presel_events)*Nbunches*LHCFreq, getErr(l1_exclusive_events,presel_events)*Nbunches*LHCFreq)
    print "L1_ETM120 rate at 55 PU:                              %f +/- %f kHz"%(float(l1etm120_events)/float(presel_events)*Nbunches*LHCFreq, getErr(l1etm120_events,presel_events)*Nbunches*LHCFreq)
    print "OR(L1_ETM120) rate at 55 PU:                          %f +/- %f kHz"%(float(l1_or_events)/float(presel_events)*Nbunches*LHCFreq, getErr(l1_or_events,presel_events)*Nbunches*LHCFreq)
    

        
def getErr(num, denom):
    if num > 0 and denom > 0: 
        return ( float(num)/float(denom) ) * math.sqrt( 1./float(num) + 1./float(denom) )
    else:
        return 0
    
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
    parser.add_option('--presel', action='store_true', dest='presel', default=False, help='turn on preselection')
    parser.add_option('--l1', dest='l1CutName', default = 'L1_ETM100',help='l1 cut', metavar='L1CutName')
    parser.add_option('-o','--odir', dest='odir', default = 'eff',help='directory to write plots', metavar='odir')
        
    (options, args) = parser.parse_args()

    main(options, args)


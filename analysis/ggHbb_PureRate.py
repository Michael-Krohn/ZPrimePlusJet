import ROOT
from VBFutils import Sort,GetVariablesToFill

PU_max = 60.
PU_min = 50
lumi = 2.E34

class Jet:
    def __init__(self, pt, eta, phi, mass, csv):
        self.pt = pt
        self.eta = eta
        self.phi = phi
        self.mass = mass
        self.csv = csv

def getNum(tree):
    num = 0.
    den = 0.
    tree.SetBranchStatus("*",0)
#    tree.SetBranchStatus("pu",1)
#    tree.SetBranchStatus("ptHat",1)
#    tree.SetBranchStatus("maxPUptHat",1)
    tree.SetBranchStatus("l1Jets_num",1)
    tree.SetBranchStatus("l1Jets_pt[l1Jets_num]",1)
    tree.SetBranchStatus("l1HT",1)
    tree.SetBranchStatus("caloJets_num",1)
    tree.SetBranchStatus("caloJets_pt[caloJets_num]",1)
    tree.SetBranchStatus("caloJets_eta[caloJets_num]",1)
    tree.SetBranchStatus("caloJets_phi[caloJets_num]",1)
    tree.SetBranchStatus("caloJets_mass[caloJets_num]",1)
    tree.SetBranchStatus("caloJets_csv[caloJets_num]",1)
    tree.SetBranchStatus("AK8CaloJetsIDPassed_pt[AK8CaloJetsIDPassed_num]",1)
    tree.SetBranchStatus("AK8CaloJetsIDPassed_eta[AK8CaloJetsIDPassed_num]",1)
    tree.SetBranchStatus("AK8PFJetsMatchedToCaloJets300_pt[AK8PFJetsMatchedToCaloJets300_num]",1)
    tree.SetBranchStatus("AK8PFJetsMatchedToCaloJets300_eta[AK8PFJetsMatchedToCaloJets300_num]",1)
    tree.SetBranchStatus("AK8PFJetsMatchedToCaloJets280_pt[AK8PFJetsMatchedToCaloJets280_num]",1)
    tree.SetBranchStatus("AK8PFJetsMatchedToCaloJets280_eta[AK8PFJetsMatchedToCaloJets280_num]",1)
    tree.SetBranchStatus("AK8PFJetsMatchedToCaloJets270_pt[AK8PFJetsMatchedToCaloJets270_num]",1)
    tree.SetBranchStatus("AK8PFJetsMatchedToCaloJets270_eta[AK8PFJetsMatchedToCaloJets270_num]",1)
    tree.SetBranchStatus("AK8PFJetsMatchedToCaloJets450_pt[AK8PFJetsMatchedToCaloJets450_num]",1)
    tree.SetBranchStatus("AK8PFJetsMatchedToCaloJets450_eta[AK8PFJetsMatchedToCaloJets450_num]",1)
    tree.SetBranchStatus("pfJets_pt[pfJets_num]",1)
    tree.SetBranchStatus("pfJets_eta[pfJets_num]",1)
    tree.SetBranchStatus("pfJets_phi[pfJets_num]",1)
    tree.SetBranchStatus("pfJets_mass[pfJets_num]",1)
    tree.SetBranchStatus("pfJets_csv[pfJets_num]",1)
    tree.SetBranchStatus("pfJets_num",1)
    tree.SetBranchStatus("HighestAK8PFJetForBTag_csv",1)
    tree.SetBranchStatus("run",1)
    tree.SetBranchStatus("HLT_AK8PFJet360_TrimMass30_v9",1)
#    tree.SetBranchStatus("HLT_AK8PFJet380_TrimMass30_v2",1)
#    tree.SetBranchStatus("HLT_AK8PFJet400_TrimMass30_v3",1)
#    tree.SetBranchStatus("HLT_AK8PFJet420_TrimMass30_v2",1)
    tree.SetBranchStatus("HLT_AK8PFHT750_TrimMass50_v3",1)
    tree.SetBranchStatus("HLT_AK8PFHT800_TrimMass50_v3",1)
    tree.SetBranchStatus("HLT_AK8PFHT850_TrimMass50_v2",1)
    tree.SetBranchStatus("HLT_AK8PFHT900_TrimMass50_v2",1)
    tree.SetBranchStatus("HLT_AK8PFJet40_v7",1)
    tree.SetBranchStatus("HLT_AK8PFJet60_v6",1)
    tree.SetBranchStatus("HLT_AK8PFJet80_v6",1)
    tree.SetBranchStatus("HLT_AK8PFJet140_v6",1)
    tree.SetBranchStatus("HLT_AK8PFJet200_v6",1)
    tree.SetBranchStatus("HLT_AK8PFJet260_v7",1)
    tree.SetBranchStatus("HLT_AK8PFJet320_v7",1)
    tree.SetBranchStatus("HLT_AK8PFJet400_v7",1)
    tree.SetBranchStatus("HLT_AK8PFJet450_v7",1)
    tree.SetBranchStatus("HLT_AK8PFJet500_v7",1)
    tree.SetBranchStatus("HLT_AK8PFJet550_v2",1)
    tree.SetBranchStatus("HLT_AK8PFJetFwd40_v6",1)
    tree.SetBranchStatus("HLT_AK8PFJetFwd60_v5",1)
    tree.SetBranchStatus("HLT_AK8PFJetFwd80_v5",1)
    tree.SetBranchStatus("HLT_AK8PFJetFwd140_v5",1)
    tree.SetBranchStatus("HLT_AK8PFJetFwd200_v5",1)
    tree.SetBranchStatus("HLT_AK8PFJetFwd260_v6",1)
    tree.SetBranchStatus("HLT_AK8PFJetFwd320_v6",1)
    tree.SetBranchStatus("HLT_AK8PFJetFwd400_v6",1)
    tree.SetBranchStatus("HLT_AK8PFJetFwd450_v6",1)
    tree.SetBranchStatus("HLT_AK8PFJetFwd500_v6",1)
    tree.SetBranchStatus("HLT_PFMET110_PFMHT110_IDTight_v11",1)
    tree.SetBranchStatus("HLT_PFMET120_PFMHT120_IDTight_v11",1)
    tree.SetBranchStatus("HLT_PFMET130_PFMHT130_IDTight_v11",1)
    tree.SetBranchStatus("HLT_PFMET140_PFMHT140_IDTight_v11",1)
    tree.SetBranchStatus("HLT_PFMETTypeOne110_PFMHT110_IDTight_v3",1)
    tree.SetBranchStatus("HLT_PFMETTypeOne120_PFMHT120_IDTight_v3",1)
    tree.SetBranchStatus("HLT_PFMETTypeOne130_PFMHT130_IDTight_v3",1)
    tree.SetBranchStatus("HLT_PFMETTypeOne140_PFMHT140_IDTight_v2",1)
    tree.SetBranchStatus("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v11",1)
    tree.SetBranchStatus("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v11",1)
#    tree.SetBranchStatus("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v10",1)
#    tree.SetBranchStatus("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v10",1)
#    tree.SetBranchStatus("HLT_PFMET200_NotCleaned_v1",1)
#    tree.SetBranchStatus("HLT_PFMET200_HBHECleaned_v1",1)
#    tree.SetBranchStatus("HLT_PFMET250_HBHECleaned_v1",1)
#    tree.SetBranchStatus("HLT_PFMET300_HBHECleaned_v1",1)
#    tree.SetBranchStatus("HLT_PFMET200_HBHE_BeamHaloCleaned_v1",1)
#    tree.SetBranchStatus("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v1",1)
#    tree.SetBranchStatus("HLT_PFMET120_PFMHT120_IDTight_L1ETMnoHF_v10",1)
#    tree.SetBranchStatus("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_L1ETMnoHF_v10",1)
#    tree.SetBranchStatus("HLT_PFMETTypeOne120_PFMHT120_IDTight_L1ETMnoHF_v2",1)
#    tree.SetBranchStatus("HLT_PFMET100_PFMHT100_IDTight_BTagCaloCSV_p7_v1",1)
#    tree.SetBranchStatus("HLT_PFMET110_PFMHT110_IDTight_BTagCaloCSV_p7_v1",1)
#    tree.SetBranchStatus("HLT_PFMET120_PFMHT120_IDTight_BTagCaloCSV_p7_v1",1)
#    tree.SetBranchStatus("HLT_PFMET130_PFMHT130_IDTight_BTagCaloCSV_p7_v1",1)
#    tree.SetBranchStatus("HLT_PFMET140_PFMHT140_IDTight_BTagCaloCSV_p7_v1",1)
    tree.SetBranchStatus("HLT_PFHT430_SixJet40_BTagCSV_p080_v2",1)
    tree.SetBranchStatus("HLT_AK8PFJet360_TrimMass30_v9",1)
    tree.SetBranchStatus("HLT_AK8PFJet360_TrimMass30_v9",1)
    tree.SetBranchStatus("AK8TrimModJets_mass[AK8TrimModJets_num]",1)
    for ev in tree:
        if not ev.run==297292: continue
        den+=1


#	if len(ev.l1Jets_pt) == 0 : continue
#	if not (ev.l1Jets_pt[0] > 170): continue
#	if len(ev.AK8CaloJetsIDPassed_pt) ==0: continue
#        if not (ev.caloJets_num>0 and ev.AK8CaloJetsIDPassed_pt[0]>300 and abs(ev.AK8CaloJetsIDPassed_eta[0])<5.0): continue
#	if len(ev.AK8PFJetsMatchedToCaloJets300_pt) == 0: continue
#        if not (ev.AK8PFJetsMatchedToCaloJets300_pt[0]>310 and abs(ev.AK8PFJetsMatchedToCaloJets300_eta[0])<5.0): continue
#	if not (ev.HighestAK8PFJetForBTag_csv > 0.4): continue 
#        if len(ev.AK8TrimModJets_mass)==0:continue
#        if not (ev.AK8TrimModJets_mass[0] > 30): continue

        if(ev.HLT_AK8PFJet360_TrimMass30_v9==0): continue
#	if(ev.HLT_AK8PFJet360_TrimMass30_v9==1 or ev.HLT_AK8PFJet450_v7==1 or ev.HLT_AK8PFHT900_TrimMass50_v2==1): continue
#        if (ev.HLT_AK8PFHT750_TrimMass50_v3==1 or ev.HLT_AK8PFHT800_TrimMass50_v3==1 or ev.HLT_AK8PFHT850_TrimMass50_v2==1 or ev.HLT_AK8PFHT900_TrimMass50_v2==1 or ev.HLT_AK8PFJet40_v7==1 or ev.HLT_AK8PFJet60_v6==1 or ev.HLT_AK8PFJet80_v6==1 or ev.HLT_AK8PFJet140_v6==1 or ev.HLT_AK8PFJet200_v6==1 or ev.HLT_AK8PFJet260_v7==1 or ev.HLT_AK8PFJet320_v7==1 or ev.HLT_AK8PFJet400_v7==1 or ev.HLT_AK8PFJet450_v7==1 or ev.HLT_AK8PFJet500_v7==1 or ev.HLT_AK8PFJet550_v2 or ev.HLT_AK8PFJetFwd40_v6==1 or HLT_AK8PFJetFwd60_v5==1 or ev.HLT_AK8PFJetFwd80_v5==1 or ev.HLT_AK8PFJetFwd140_v5 or ev.HLT_AK8PFJetFwd200_v5==1 or ev.HLT_AK8PFJetFwd260_v6==1 or ev.HLT_AK8PFJetFwd320_v6==1 or ev.HLT_AK8PFJetFwd400_v6==1 or ev.HLT_AK8PFJetFwd450_v6==1 or ev.HLT_AK8PFJetFwd500_v6==1 or ev.HLT_PFMET110_PFMHT110_IDTight_v11==1 or ev.HLT_PFMET120_PFMHT120_IDTight_v11==1 or ev.HLT_PFMET130_PFMHT130_IDTight_v11==1 or ev.HLT_PFMET140_PFMHT140_IDTight_v11==1 or ev.HLT_PFMETTypeOne110_PFMHT110_IDTight_v3==1 or ev.HLT_PFMETTypeOne120_PFMHT120_IDTight_v3==1 or ev.HLT_PFMETTypeOne130_PFMHT130_IDTight_v3==1 or ev.HLT_PFMETTypeOne140_PFMHT140_IDTight_v2==1 or ev.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v11==1 or ev.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v11==1): continue
#	if (ev.HLT_AK8PFJet360_TrimMass30_v9==1 or ev.HLT_AK8PFJet380_TrimMass30_v2==1 or ev.HLT_AK8PFJet400_TrimMass30_v3==1 or ev.HLT_AK8PFJet420_TrimMass30_v2 or ev.HLT_AK8PFHT750_TrimMass50_v3==1 or ev.HLT_AK8PFHT800_TrimMass50_v3==1 or ev.HLT_AK8PFHT850_TrimMass50_v2==1 or ev.HLT_AK8PFHT900_TrimMass50_v2==1 or ev.HLT_AK8PFJet40_v7==1 or ev.HLT_AK8PFJet60_v6==1 or ev.HLT_AK8PFJet80_v6==1 or ev.HLT_AK8PFJet140_v6==1 or ev.HLT_AK8PFJet200_v6==1 or ev.HLT_AK8PFJet260_v7==1 or ev.HLT_AK8PFJet320_v7==1 or ev.HLT_AK8PFJet400_v7==1 or ev.HLT_AK8PFJet450_v7==1 or ev.HLT_AK8PFJet500_v7==1 or ev.HLT_AK8PFJet550_v2 or ev.HLT_AK8PFJetFwd40_v6==1 or HLT_AK8PFJetFwd60_v5==1 or ev.HLT_AK8PFJetFwd80_v5==1 or ev.HLT_AK8PFJetFwd140_v5 or ev.HLT_AK8PFJetFwd200_v5==1 or ev.HLT_AK8PFJetFwd260_v6==1 or ev.HLT_AK8PFJetFwd320_v6==1 or ev.HLT_AK8PFJetFwd400_v6==1 or ev.HLT_AK8PFJetFwd450_v6==1 or ev.HLT_AK8PFJetFwd500_v6==1 or ev.HLT_PFMET110_PFMHT110_IDTight_v11==1 or ev.HLT_PFMET120_PFMHT120_IDTight_v11==1 or ev.HLT_PFMET130_PFMHT130_IDTight_v11==1 or ev.HLT_PFMET140_PFMHT140_IDTight_v11==1 or ev.HLT_PFMETTypeOne110_PFMHT110_IDTight_v3==1 or ev.HLT_PFMETTypeOne120_PFMHT120_IDTight_v3==1 or ev.HLT_PFMETTypeOne130_PFMHT130_IDTight_v3==1 or ev.HLT_PFMETTypeOne140_PFMHT140_IDTight_v2==1 or ev.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v11==1 or ev.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v11==1 or ev.HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v10==1 or ev.HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v10==1 or ev.HLT_PFMET200_NotCleaned_v1==1 or ev.HLT_PFMET200_HBHECleaned_v1==1 or ev.HLT_PFMET250_HBHECleaned_v1==1 or ev.HLT_PFMET300_HBHECleaned_v1==1 or ev.HLT_PFMET200_HBHE_BeamHaloCleaned_v1==1 or ev.HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v1==1 or ev.HLT_PFMET120_PFMHT120_IDTight_L1ETMnoHF_v10==1 or ev.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_L1ETMnoHF_v10==1 or ev.HLT_PFMETTypeOne120_PFMHT120_IDTight_L1ETMnoHF_v2==1 or ev.HLT_PFMET100_PFMHT100_IDTight_BTagCaloCSV_p7_v1==1 or ev.HLT_PFMET110_PFMHT110_IDTight_BTagCaloCSV_p7_v1 or ev.HLT_PFMET120_PFMHT120_IDTight_BTagCaloCSV_p7_v1==1 or ev.HLT_PFMET130_PFMHT130_IDTight_BTagCaloCSV_p7_v1==1 or ev.HLT_PFMET140_PFMHT140_IDTight_BTagCaloCSV_p7_v1==1 or ev.HLT_PFHT430_SixJet40_BTagCSV_p080_v2==1): continue
#        calojetswithcsv = []
#        for (pt,eta,phi,mass,csv) in zip(ev.caloJets_pt,ev.caloJets_eta,ev.caloJets_phi,ev.caloJets_mass,ev.caloJets_csv):
#            jet = ROOT.TLorentzVector()
#            jet.SetPtEtaPhiM(pt,eta,phi,mass)
#            jet.csv = csv
#            calojetswithcsv.append(jet)
        
#        pfjetswithcsv = []
#        for (pt,eta,phi,mass,csv) in zip(ev.pfJets_pt,ev.pfJets_eta,ev.pfJets_phi,ev.pfJets_mass,ev.pfJets_csv):
#            jet = ROOT.TLorentzVector()
#            jet.SetPtEtaPhiM(pt,eta,phi,mass)
#            jet.csv = csv
#            pfjetswithcsv.append(jet)

#        (b1,b2,q1,q2) = Sort(calojetswithcsv,'Eta')
#        (Detaqq_eta,Dphibb_eta,Mqq_eta,Mbb_eta) = GetVariablesToFill(b1,b2,q1,q2)
        
#        (b1,b2,q1,q2) = Sort(pfjetswithcsv,'1BTagAndEta')
#        (Detaqq_1b,Dphibb_1b,Mqq_1b,Mbb_1b) = GetVariablesToFill(b1,b2,q1,q2)
        
#        (b1,b2,q1,q2) = Sort(pfjetswithcsv,'2BTagAndPt')
#        (Detaqq_2b,Dphibb_2b,Mqq_2b,Mbb_2b) = GetVariablesToFill(b1,b2,q1,q2)
        
#        if not(Detaqq_eta>1.5 and Mqq_eta>150): continue
        
        
#        pfCSV = [csv for csv in ev.pfJets_csv]
#        pfCSV.sort(reverse=True)
#        if not (pfCSV[0]>0.82): continue
#        if not(Detaqq_1b>4.1 and Mqq_1b>500 and Dphibb_1b<1.6): continue
        
#        if not (pfCSV[0]>0.82 and pfCSV[1]>0.47): continue
#        if not(Detaqq_2b>2.3 and Mqq_2b>240 and Dphibb_2b<2.1): continue
        num += 1

    print "Numerator: ", num   
    print "Denominator: ", den 
    return num,den

def getFraction(fileName):
    file_ = ROOT.TFile(fileName)
    tree = file_.Get("tree")
#    den = file_.Get("Count").GetBinContent(1)
#    den *= (PU_max-PU_min)/(63-28)
    num,den = getNum(tree)
    if num<=0: return (0.,0.)
    file_.Close()
    fraction = num/den
    return (fraction , (fraction*(1-fraction)/den)**0.5)

def getRate(fileName, xsection, lumi):
    sampleRate = 1.E-36*xsection*lumi
    (fraction, errFraction) = getFraction(fileName)
    return (sampleRate * fraction , sampleRate * errFraction)

## https://github.com/cms-steam/RateEstimate/blob/master/datasetCrossSections/datasetCrossSectionsSummer16.py#L66
bkgCrossSection = [
#    ("QCD15to30.root"   , 1837410000.),
#    ("QCD30to50.root"   , 140932000.),
    ("QCD50to80.root"   , 19204300.),
    ("QCD80to120.root"  , 2762530.),
    ("QCD120to170.root" , 471100.),
    ("QCD170to300.root" , 117276.),
    ("QCD300to470.root" , 7823.),
#    ("QCD470to600.root" , 648.2),
]

totalRate = 0
totalRateError = 0

#for (fileName,xsection) in bkgCrossSection:
#    rate , rateError = getRate(fileName, xsection, lumi)
#    print fileName,"\t",rate," +/- ",rateError
#    totalRate += rate
#    totalRateError += rateError
#
#print "totalRate\t",totalRate," +/- ",totalRateError

for signal in ["ParkedData_Sklims/HLTPhysics_2017B_Run297292to297296.root"]:
    (fraction, errFraction) = getFraction(signal)
    print signal,"\t",fraction," +/- ",errFraction
    print signal, "Rate(90kHz L1 Rate)\t",fraction*75000.,"+/-",errFraction*75000.


import ROOT
from VBFutils import Sort,GetVariablesToFill

PU_max = 60.
PU_min = 50.
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
    num2 = 0.
    num3 = 0.
    tree.SetBranchStatus("*",0)
    tree.SetBranchStatus("maxPUptHat",1)
    tree.SetBranchStatus("l1Jets_num",1)
    tree.SetBranchStatus("l1Jets_pt[l1Jets_num]",1)
    tree.SetBranchStatus("l1HT",1)
    tree.SetBranchStatus("caloJets_num",1)
    tree.SetBranchStatus("caloJets_pt[caloJets_num]",1)
    tree.SetBranchStatus("caloJets_eta[caloJets_num]",1)
    tree.SetBranchStatus("caloJets_phi[caloJets_num]",1)
    tree.SetBranchStatus("caloJets_mass[caloJets_num]",1)
    tree.SetBranchStatus("caloJets_csv[caloJets_num]",1)
    tree.SetBranchStatus("pfJets_pt[pfJets_num]",1)
    tree.SetBranchStatus("pfJets_eta[pfJets_num]",1)
    tree.SetBranchStatus("pfJets_phi[pfJets_num]",1)
    tree.SetBranchStatus("pfJets_mass[pfJets_num]",1)
    tree.SetBranchStatus("pfJets_csv[pfJets_num]",1)
    tree.SetBranchStatus("pfJets_num",1)
    tree.SetBranchStatus("run",1)
    for ev in tree:
#        if not ev.run==297292: continue
        den+=1

        if not (ev.caloJets_num>3 and ev.caloJets_pt[0]>80 and ev.caloJets_pt[1]>65 and ev.caloJets_pt[2]>50 and ev.caloJets_pt[3]>15): continue
        if not (max(ev.caloJets_csv)>0.8): continue
        if not (ev.pfJets_num>3 and ev.pfJets_pt[0]>92 and ev.pfJets_pt[1]>76 and ev.pfJets_pt[2]>64 and ev.pfJets_pt[3]>15): continue
        
        calojetswithcsv = []
	num2 += 1
	if not (ev.pfJets_pt[0] < 9999 and ev.pfJets_pt[1] < 9999 and ev.pfJets_pt[2] < 9999 and ev.pfJets_pt[3] < 9999): continue
	num3 +=1
        for (pt,eta,phi,mass,csv) in zip(ev.caloJets_pt,ev.caloJets_eta,ev.caloJets_phi,ev.caloJets_mass,ev.caloJets_csv):
            jet = ROOT.TLorentzVector()
            jet.SetPtEtaPhiM(pt,eta,phi,mass)
            jet.csv = csv
            calojetswithcsv.append(jet)
        
        pfjetswithcsv = []
        for (pt,eta,phi,mass,csv) in zip(ev.pfJets_pt,ev.pfJets_eta,ev.pfJets_phi,ev.pfJets_mass,ev.pfJets_csv):
            jet = ROOT.TLorentzVector()
            jet.SetPtEtaPhiM(pt,eta,phi,mass)
            jet.csv = csv
            pfjetswithcsv.append(jet)

        (b1,b2,q1,q2) = Sort(calojetswithcsv,'Eta')
        (Detaqq_eta,Dphibb_eta,Mqq_eta,Mbb_eta) = GetVariablesToFill(b1,b2,q1,q2)
        
        (b1,b2,q1,q2) = Sort(pfjetswithcsv,'1BTagAndEta')
        (Detaqq_1b,Dphibb_1b,Mqq_1b,Mbb_1b) = GetVariablesToFill(b1,b2,q1,q2)
        
        (b1,b2,q1,q2) = Sort(pfjetswithcsv,'2BTagAndPt')
        (Detaqq_2b,Dphibb_2b,Mqq_2b,Mbb_2b) = GetVariablesToFill(b1,b2,q1,q2)
        
        if not(Detaqq_eta>1.5 and Mqq_eta>150): continue
        
        
        pfCSV = [csv for csv in ev.pfJets_csv]
        pfCSV.sort(reverse=True)
        if not (pfCSV[0]>0.82): continue
        if not(Detaqq_1b>4.1 and Mqq_1b>500 and Dphibb_1b<1.6): continue
        
#        if not (pfCSV[0]>0.82 and pfCSV[1]>0.47): continue
#        if not(Detaqq_2b>2.3 and Mqq_2b>240 and Dphibb_2b<2.1): continue
        num += 1

    print "Numerator: ", num
    print "Numerator2: ", num2
    print "Numerator3: ", num3
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
#    ("QCD50to80.root"   , 19204300.),
    ("../QCD80to120.root"  , 2762530.),
    ("../QCD120to170.root" , 471100.),
    ("../QCD170to300.root" , 117276.),
    ("../QCD300to470.root" , 7823.),
    ("../QCD470to600.root" , 648.2),
]

totalRate = 0
totalRateError = 0

#for (fileName,xsection) in bkgCrossSection:
#    rate , rateError = getRate(fileName, xsection, lumi)
#    print fileName,"\t",rate," +/- ",rateError
#    totalRate += rate
#    totalRateError += rateError

#print "totalRate\t",totalRate," +/- ",totalRateError

for signal in ["ParkedData_Sklims/HLTPhysics_2017C.root"]:
#for signal in ["ParkedData_Sklims/HLTPhysics_2017C.root","ParkedData_Sklims/HLTPhysics_2017B_Run297292to297296.root","ParkedData_Sklims/HLTPhysics_2017B_Run297424to297435.root"]:
#for signal in ["ParkedData_Sklims/HLTPhysics_2017B_Run297292to297296.root"]:
    (fraction, errFraction) = getFraction(signal)
    print signal,"\t",fraction," +/- ",errFraction
    print signal, "Rate(54kHz L1 Rate)\t",fraction*54000.,"+/-",errFraction*54000.

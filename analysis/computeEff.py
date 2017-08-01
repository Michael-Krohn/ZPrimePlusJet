import ROOT as rt, sys, math, os
##-------------------------------------------------------------------------------------



lumi = 2.E34
pu_min = 0.
pu_max = 100.

def getFraction(num, den):
    fraction = 1.*num/den
    if num<=0: return (0.,0.)
    print fraction
    print den
    print fraction*(1.-fraction)/den
    return (fraction , (fraction*(1.-fraction)/den)**0.5)

def getRate(xsection, num, den):    
    sampleRate = 1.E-36*xsection*lumi
    (fraction, errFraction) = getFraction(num, den)
    return (sampleRate * fraction , sampleRate * errFraction)
    
if __name__ == '__main__':

    rt.gROOT.SetBatch()

    rt.gROOT.ProcessLine(".L DeltaPhi.C+")
    rt.gSystem.Load("DeltaPhi_C.so")
    #bkg_names = ["QCD15to30","QCD30to50","QCD50to80", "QCD80to120", "QCD120to170", "QCD170to300", "QCD300to470", "QCD470to600"]
    #sig_names = ["gluglu260_Trigger","gluglu300_Trigger", "gluglu450_Trigger","ZH_Trigger"]
    
    bkg_names = ['QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8',
                 'QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8',
                 'QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8',
                 'QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8',
                 'QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8',
                 'QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8',
                 'QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8',                 
                 'QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8']
    sig_names = ['tree_11']
#    sig_names = ['ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8']
    names = sig_names
#    names = bkg_names+sig_names
#    files = ['root://cmseos.fnal.gov//store/user/mkrohn/HLT_Ntuple_Hbb_Signal_v8_AODSIM_hltScript/GluGluHToBB_M125_13TeV_powheg_pythia8/HLT_Ntuple_Hbb_Signal_v8_AODSIM_hltScript/170626_172957/0000/'+name+'.root' for name in bkg_names]
    files = ['.root' for name in sig_names]
    files = ['root://cmseos.fnal.gov//eos/uscms/store/user/mkrohn/HLT_Ntuple_Hbb_Signal_v8_AODSIM_FINAL_6/GluGluHToBB_M125_13TeV_powheg_pythia8/HLT_Ntuple_Hbb_Signal_v8_AODSIM_FINAL_6/170717_191655/0000/'+name+'.root' for name in sig_names]
#    files.extend(['root://cmseos.fnal.gov//store/user/mkrohn/HLT_Ntuple_Hbb_Signal_v8_AODSIM_hltScript/GluGluHToBB_M125_13TeV_powheg_pythia8/HLT_Ntuple_Hbb_Signal_v8_AODSIM_hltScript/170626_172957/0000/'+name+'.root' for name in sig_names])
    xsec ={}
    for name in sig_names: xsec[name] = 0.
    xsec["QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8"] = 1837410000.
    xsec["QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8"] = 140932000.
    xsec["QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8"] = 19204300.
    xsec["QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8"] = 2762530.
    xsec["QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8"] = 471100.
    xsec["QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8"] = 117276.
    xsec["QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8"] = 7823.
    xsec["QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8"] = 648.2
     
    mysum=0
    eff = {}
    eff_error = {}
    ratesL1 = {}
    ratesHLT = {}
    errorL1 = {}
    errorHLT = {}
    total_rate_HLT = {}
    total_rate_L1 = {}
    total_error_HLT = {}
    total_error_L1 = {}

    triggers = {}
    #l1MetString = "l1l1Met>120. || ( Sum$(l1Jets_pt>60 && abs(l1Jets_eta)<2.5)>=2 && l1Met > 90) || ( Alt$(l1Jets_pt[0],0)>=60 && l1Met > 100 && abs(mymath::deltaPhi(l1Met_phi,Alt$(l1Jets_phi[0],l1Met_phi))) < 0.4)"

    l1MetString = "1"
    l1PFJet360String = "l1Jets_pt>180"
    l1PFAK8BTagString = "l1Jets_pt>170"
    l1HTString  = "l1HT > 380"
    ## triggers['HLT_PFMET100_PFMHT100_BTagCaloCSV_p7'] = [l1MetString,
    ##                                      "pfMet_pt > 100 && pfMht_pt > 100. && Sum$(caloJets_csv>0.7 && caloJets_pt>30 && abs(caloJets_eta)<2.4)>=1"]
    ## triggers['HLT_PFMET110_PFMHT110_BTagCaloCSV_p7'] = [l1MetString,
    ##                                      "pfMet_pt > 110 && pfMht_pt > 110. && Sum$(caloJets_csv>0.7 && caloJets_pt>30 && abs(caloJets_eta)<2.4)>=1"]
    ## triggers['HLT_PFMET120_PFMHT120_BTagCaloCSV_p7'] = [l1MetString,
    ##                                      "pfMet_pt > 120 && pfMht_pt > 120. && Sum$(caloJets_csv>0.7 && caloJets_pt>30 && abs(caloJets_eta)<2.4)>=1"]
    ## triggers['HLT_PFMET130_PFMHT130_BTagCaloCSV_p7'] = [l1MetString,
    ##                                      "pfMet_pt > 130 && pfMht_pt > 130. && Sum$(caloJets_csv>0.7 && caloJets_pt>30 && abs(caloJets_eta)<2.4)>=1"]
    ## triggers['HLT_PFMET140_PFMHT140_BTagCaloCSV_p7'] = [l1MetString,
    ##                                      "pfMet_pt > 140 && pfMht_pt > 140. && Sum$(caloJets_csv>0.7 && caloJets_pt>30 && abs(caloJets_eta)<2.4)>=1"]
    ## triggers['HLT_PFMET100_PFMHT100'] = [l1MetString,
    ##                                      "pfMet_pt > 100 && pfMht_pt > 100."]
    ## triggers['HLT_PFMET110_PFMHT110'] = [l1MetString,
    ##                                      "pfMet_pt > 110 && pfMht_pt > 110."]
    triggers['HLT_AK8PFJet360_TrimMass30'] = [l1PFJet360String,
                                        "AK8CaloJetsIDPassed_pt>260 && abs(AK8CaloJetsIDPassed_eta) < 2.5 && AK8PFJetsCHS_pt > 360 && AK8TrimModJets_mass > 30. && abs(AK8PFJetsCHS_eta) < 2.5"]
    triggers['HLT_AK8PFHT750_TrimMass50'] = [l1HTString,
                                        "AK8PFJetsCHS_pt > 360 && AK8TrimModJets_mass > 50. && abs(AK8PFJetsCHS_eta) < 2.5"]
    triggers['HLT_AK8PFJet550_PFAK8BTagCSV'] = [l1PFAK8BTagString,
                                        "AK8CaloJetsIDPassed_pt>500 && abs(AK8CaloJetsIDPassed_eta) < 2.5 && AK8PFJetsMatchedToCaloJets500_pt > 550 && abs(AK8PFJetsMatchedToCaloJets500_eta) < 5.0 && AK8PFJetsCHS_csv > 0.63"]
#    triggers['HLT_PFMET120_PFMHT120'] = [l1MetString,
#                                         "pfMet_pt > 120 && pfMht_pt > 120."]
    triggers['HLT_PFMET120_PFMHT120_HT60'] = [l1MetString,
                                         "pfMet_pt > 120 && pfMht_pt > 120. && Sum$(pfJets_pt*(pfJets_eta<3)*(pfJets_pt>30))> 60."]
    ## triggers['HLT_PFMET130_PFMHT130'] = [l1MetString,
    ##                                      "pfMet_pt > 130 && pfMht_pt > 130."]
    ## triggers['HLT_PFMET140_PFMHT140'] = [l1MetString,
    ##                                      "pfMet_pt > 140 && pfMht_pt > 140."]
    ## triggers['OR'] = [l1MetString,
    ##                                      "(pfMet_pt > 120 && pfMht_pt > 120.) || ( pfMet_pt > 110 && pfMht_pt > 110. && Sum$(caloJets_csv>0.7 && caloJets_pt>30 && abs(caloJets_eta)<2.4)>=1 )"]



    for trigger in triggers:
        total_rate_HLT[trigger] = 0
        total_rate_L1[trigger] = 0
        total_error_HLT[trigger] = 0
        total_error_L1[trigger] = 0
        
    for name, ifile, in zip(names,files):
        print "reading %s"%(ifile)
        for trigger, (sel_L1, sel_HLT) in triggers.iteritems():
            f = rt.TFile(ifile)
            t = f.Get("tree")

            sel_pu = 'pu >= %f && pu <= %f'%(pu_min, pu_max)
            if name not in sig_names:
                #sel_pu = 'pu >= %f && pu <= %f'%(pu_min, pu_max)
                sel_pu = 'ptHat>maxPUptHat && pu >= %f && pu <= %f'%(pu_min, pu_max)

                
            if name == 'QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8' or name== 'QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8':
                total_events = f.Get('Count').GetBinContent(1)
                print total_events
            else:
                total_events = t.Draw("",sel_pu)
            presel = '1==1'
#            presel = 'offMet_pt>130 && Sum$(offJets_csv>0.5426 && abs(offJets_eta)<2.5 && offJets_pt>25)>=2 && Sum$(offJets_csv>0.8484 && abs(offJets_eta)<2.5 && offJets_pt>25)>=1'
            if 'QCD' in name:
                total_presel_events = total_events
            else:
                total_presel_events = t.Draw("","(" + presel + ") && (" + sel_pu + ")")
            #h = f.Get("Count")
            #total_events = h.GetBinContent(1)

            sel_L1_tot = "(" + sel_L1 + ") && (" + sel_pu + ")"
            sel_HLT_tot = "(" + sel_L1 + ") && (" + sel_HLT + ") && (" + sel_pu + ")"
            sel_HLT_presel =  "(" + sel_L1 + ") && (" + sel_HLT + ") && (" + sel_pu + ") && ("  + presel + ")"
                
            pass_L1 = t.Draw("l1HT", sel_L1_tot)
            pass_HLT = t.Draw("l1HT",sel_HLT_tot)
            if 'QCD' in name:
                pass_HLT_presel = pass_HLT
            else:
                pass_HLT_presel = t.Draw("l1HT",sel_HLT_presel)

	    print "presel:"
	    print presel 
            print trigger, total_events, total_presel_events, pass_L1, pass_HLT
            ratesL1[name,trigger], errorL1[name,trigger]  = getRate(xsec[name],pass_L1,total_events)
            ratesHLT[name,trigger], errorHLT[name, trigger] =  getRate(xsec[name],pass_HLT,total_events)
            total_rate_L1[trigger] += ratesL1[name,trigger] 
            total_rate_HLT[trigger] += ratesHLT[name,trigger]
            total_error_L1[trigger] += errorL1[name,trigger] 
            total_error_HLT[trigger] += errorHLT[name,trigger]
            eff[name,trigger], eff_error[name, trigger]  = getFraction(pass_HLT_presel,total_presel_events)


    #print "eff", eff
    #print "rates @ HLT", ratesHLT
    #signal = 'ZH_Trigger'
    signal = 'tree_11'
#    signal = 'ggHToBB'
    print "trigger | total rate @ L1 | total rate @ HLT | eff for signal %s"%signal
    for trigger in sorted(triggers):
        print "%s | %.2f +- %.2f kHz | %.2f +- %.2f Hz | %f +- %f%%"%(trigger,
                                                              total_rate_L1[trigger]/1000.,
                                                              total_error_L1[trigger]/1000.,
                                                              total_rate_HLT[trigger],
                                                              total_error_HLT[trigger],
                                                              eff[signal,trigger]*100.,
                                                              eff_error[signal,trigger]*100.)



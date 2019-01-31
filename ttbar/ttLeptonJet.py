import os,sys
basepath = os.path.abspath(__file__).rsplit('/ttCMSDAS/',1)[0]+'/ttCMSDAS/'
sys.path.append(basepath)

from framework.analysis import analysis
from framework.functions import DeltaPhi, DiPt, InvMass, lepton, jet
from ROOT.TMath import Sqrt as sqrt
from ROOT import *

################ Analysis
class ttLeptonJet(analysis):
  def init(self):
    # Load SF files
    if not self.isData:
      self.LoadHisto('MuonIsoSF', basepath+'./inputs/MuonISO.root', 'NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta') # pt, abseta
      self.LoadHisto('MuonIdSF',  basepath+'./inputs/MuonID.root',  'NUM_TightID_DEN_genTracks_pt_abseta') # pt, abseta
      self.LoadHisto('ElecSF',    basepath+'./inputs/ElecTightCBid94X.root',  'EGamma_SF2D') # eta, pt

    # Objects for the analysis
    self.selLeptons = []
    self.selJets = []
    self.pmet = TLorentzVector()
    self.selLeptons_QCD = []
    self.selJets_QCD = []
    self.pmet_QCD = TLorentzVector()
    self.selLeptons_noMet = []
    self.selJets_noMet = []
    self.pmet_noMet = TLorentzVector()
    self.selLeptons_QCD_noMet = []
    self.selJets_QCD_noMet = []
    self.pmet_QCD_noMet = TLorentzVector()

    # Create output histograms
    for selN in ("signal", "qcd", "signal_nomet","qcd_nomet"):
        self.myCreateTH1F(selN, "Lep0Pt", "", 24, 0, 120)
        ## etc...
	self.CreateTH1F("Lep0Pt",   "", 24, 0, 120)
    	self.CreateTH1F("Lep0Eta",  "", 50, -2.5, 2.5)

  def myCreateTH1F(self, sel, var, title, nBin, xMin, xMax):
      self.CreateTH1F("{0}_{1}".format(sel, var), title, nBin, xMin, xMax)
0
.
  def resetObjects(self):
    ''' Reset the list where the objects are stored '''
    self.selLeptons = []
    self.selJets = []
    self.pmet = TLorentzVector()

  def FillHistograms(self, leptons, jets, pmet, sel="signal"):
    ''' Fill all the histograms. Take the inputs from lepton list, jet list, pmet '''
    self.obj['DeltaPhi'].Fill(dphi/3.141592, self.weight)
    self.weight = self.EventWeight * self.SFmuon * self.SFelec * self.PUSF

    # Re-calculate the observables
    if not len(leptons) == 1: return # Just in case
    lep0  = leptons[0]
    l0pt  = lep0.Pt()
    l0eta = lep0.Eta()
    
    ### Fill the histograms
    self.obj['{0}_Lep0Pt'.format(sel)].Fill(l0pt, self.weight)
    self.obj['Lep0Eta'].Fill(l0eta, self.weight)

  def insideLoop(self, t):
    self.resetObjects()
    signalPass1 = False
    signalPass2 = False

    ### Lepton selection
    ###########################################
    if not self.isData: nGenLep = t.nGenDressedLepton 
    
    ##### Muons
    QCDPass = False
    NoMetPass = False	
    for i in range(t.nMuon):
      p = TLorentzVector()
      p.SetPtEtaPhiM(t.Muon_pt[i], t.Muon_eta[i], t.Muon_phi[i], t.Muon_mass[i])
      charge = t.Muon_charge[i]

      # Tight ID, tight ISO, RelIso04 < 0.15, tight IP
      if not t.Muon_tightId[i]: continue # Tight ID
	dxy = abs(t.Muon_dxy[i])
      	dz  = abs(t.Muon_dz[i] )
     	if dxy > 0.05 or dz > 0.1: continue

      	# pT > 30 GeV, |eta| < 2.1
      	if p.Pt() < 30 or abs(p.Eta()) > 2.1: continue

      if t.Muon_pfRelIso04_all[i] > 0.15 and t.MET_pt < 20 :
	self.selLeptons_QCD_noMet.append(lepton(p, charge, 13))
      elif t.Muon_pfRelIso04_all[i] < 0.15 and t.MET_pt < 20 :
	self.selLeptons__noMet.append(lepton(p, charge, 13))
      elif t.Muon_pfRelIso04_all[i] < 0.15 and t.MET_pt > 20 :
      	self.selLeptons.append(lepton(p, charge, 13))
      else:
	self.selLeptons_QCD.append(lepton(p,charge,13))


    ##### Electrons
    for i in range(t.nElectron):
      p = TLorentzVector()
      p.SetPtEtaPhiM(t.Electron_pt[i], t.Electron_eta[i], t.Electron_phi[i], t.Electron_mass[i])
      charge = t.Electron_charge[i]
      etaSC    = abs(p.Eta());
      convVeto = t.Electron_convVeto[i]

      # Tight cut-based Id, convVeto, RelIso03 tight, tight IP
      if not t.Electron_cutBased[i] >= 4: continue
      if not convVeto: continue
      relIso03 = t.Electron_pfRelIso03_all[i]
       
      dxy = abs(t.Electron_dxy[i])
      dz  = abs(t.Electron_dz[i] )
      if dxy > 0.05 or dz > 0.1: continue

      # pT > 30 GeV, |eta| < 2.1
      if p.Pt() < 30 or abs(p.Eta()) > 2.1: continue
      if   etaSC <= 1.479 and relIso03 > 0.0361: signalPass1 = True
      elif etaSC >  1.479 and relIso03 > 0.094:  signalPass2 = True
      if signalPass1 and signalPass2 and t.MET_pt < 20 :
      	self.selLeptons_noMet.append(lepton(p, charge, 11))
      elif signalPass1 is False and signalPass2 is False and t.MET_pt < 20 :
      	self.selLeptons_QCD_noMet.append(lepton(p, charge, 11))
      elif signalPass1 is False and signalPass2 is False and t.MET_pt > 20 :
      	self.selLeptons_QCD.append(lepton(p, charge, 11))
      else:self.selLeptons.append(lepton(p, charge, 11))





    leps = self.selLeptons
    pts  = [lep.Pt() for lep in leps]
    self.selLeptons = [lep for _,lep in sorted(zip(pts,leps))]
	
    #### Jets
    for i in range(t.nJets):
    	p = TLorentzVector()
	p.SetPtEtaPhiM(t.Jet_pt[i], t.Jet_eta[i], t.Jet_phi[i], t.Jet_mass[i])
	### loose PFJetID is applied in order to reject jets induced by pure instrumental noise in the calorimeters

	if t.jetId & (0x1<<1)

        ### pT>30GeV and |eta|< 2.5
	if p.Pt() < 30 or abs(p.Eta()) > 2.5: continue
        self.selJets.append(p)

	
    ### Calculate the weights
    self.SFelec = 1; self.SFmuon = 1; self.SFelecErr = 0; self. SFmuonErr = 0
    if not self.isData:
      for lep in self.selLeptons:
        if lep.IsMuon():
          sf, err = self.GetSFandErr('MuonIsoSF, MuonIdSF', lep.Pt(), TMath.Abs(lep.Eta()))
          self.SFmuon*=sf
          self.SFmuonErr+=err*err
        else:
          sf, err = self.GetSFandErr('ElecSF', lep.Eta(), lep.Pt())
          self.SFelec*=sf
          self.SFelecErr+=err*err
      self.SFelecErr = sqrt(self.SFelecErr)
      self.SFmuonErr = sqrt(self.SFmuonErr)

    # PU SF --> PLEASE CHECK THAT THE WEIGHTS ARE IN THE TREES THAT YOU'RE USING!
    if not self.isData:
      self.PUSF   = t.puWeight
      self.PUUpSF = t.puWeightUp
      self.PUDoSF = t.puWeightDown
    else:
      self.PUSF   = 1; self.PUUpSF = 1; self.PUDoSF = 1

    ### Event selection
    ###########################################

    ### SIGNAL
 
    ### One electron or muon and at least 4 jets
    ### veto events with more than one lepton.
    ### veto events containing leptons passing relaxed ID/ISO cuts. 
    if not len(leps) == 1; 	return
    if not len(selJets)>=4; 	return
     
    l0 = leps[0]

    ### Fill the histograms
    self.FillHistograms(self.selLeptons, self.selJets, self.pmet)


    ### QCD selection : non-prompt and less isolated leptons
	### the control region : isolation cut is reversed 
    
	

   self.FillHistograms(self.selLeptons_QCD, self.selJets, self.pmet)

   self.FillHistograms(self.selLeptons_QCD_noMet, self.selJets, self.pmet)
	
   self.FillHistograms(self.selLeptons_noMet, self.selJets, self.pmet)

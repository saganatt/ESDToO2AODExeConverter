/*
 * Taken from macro of Peter Hristov
*/

#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDMuonTrack.h"

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TStopwatch.h>
#include <TRandom.h>
#include <TChain.h>

#include <iostream>
using namespace std;

UInt_t numberOfSetBits(UInt_t i) {
  // Function that calculates the Hemming weight of an UInt_t
  // see
  // https://en.wikipedia.org/wiki/Hamming_weight
  // https://stackoverflow.com/questions/109023/how-to-count-the-number-of-set-bits-in-a-32-bit-integer
  i = i - ((i >> 1) & 0x55555555);
  i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
  return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}
void printBits(UInt_t ii) {
  // Prints the bits of "ii" in increasing order
  for (UInt_t j=0; j<32; ++j)
    putchar(ii & (1 << j) ? '1' : '0');
  putchar('\n');
}
void hemw(UInt_t iii){
  // Test the function that calculates the number of set bits
  printf("Value %i\n", iii);
  printf("Bits ");
  printBits(iii);
  printf("Number of set bits %i\n", numberOfSetBits(iii));
}

Float_t float2half(Float_t x, UInt_t mask = 0xFFFFFFFF){
  // Mask the less significant bits in the float fraction (1 bit sign, 8 bits exponent, 23 bits fraction)
  // mask 0xFFFF0000 means only 23 - 16 = 7 bits in the fraction
  union { Float_t y; UInt_t iy;} myu;
  myu.y = x;
  myu.iy &= mask;
  return myu.y;
}

void convertO2(TTree * tEsd) {
  //---------------------------------------------------------------------------

  // ESD Initialisation
  AliESDEvent * esd = new AliESDEvent();
  esd->ReadFromTree(tEsd);

  // File to store the results
  TFile * fO2 = TFile::Open("O2aod.root","recreate");
  //---------------------------------------------------------------------------
  // Collision vertex
  //---------------------------------------------------------------------------
  TTree * tO2vtx = new TTree("O2vtx","Collision vertices");
  struct {
    // Vertex position
    Float_t           fX;         // X coordinate of the primary vertex
    Float_t           fY;         // X coordinate of the primary vertex
    Float_t           fZ;         // X coordinate of the primary vertex
    // Covariance matrix
    Float_t           fCovXX;     // vertex covariance matrix
    Float_t           fCovXY;
    Float_t           fCovXZ;
    Float_t           fCovYY;
    Float_t           fCovYZ;
    Float_t           fCovZZ;
    // Quality parameters
    Float_t           fChi2;      // Chi2 of the vertex
    UInt_t            fN;         // Number of contributors
    // Time in LHC representation
    UInt_t            fOrbitNumber;       // Orbit Number
    UInt_t            fPeriodNumber;      // Period Number
    UShort_t          fBCNumber;  // Bunch Crossing Number
  } vtx;
  // Branches
  tO2vtx->Branch("fX", &vtx.fX, "fX/F");
  tO2vtx->Branch("fY", &vtx.fY, "fY/F");
  tO2vtx->Branch("fZ", &vtx.fZ, "fZ/F");
  //
  tO2vtx->Branch("fCovXX", &vtx.fCovXX, "fCovXX/F");
  tO2vtx->Branch("fCovXY", &vtx.fCovXY, "fCovXY/F");
  tO2vtx->Branch("fCovXZ", &vtx.fCovXZ, "fCovXZ/F");
  tO2vtx->Branch("fCovYY", &vtx.fCovYY, "fCovYY/F");
  tO2vtx->Branch("fCovYZ", &vtx.fCovYZ, "fCovYZ/F");
  tO2vtx->Branch("fCovZZ", &vtx.fCovZZ, "fCovZZ/F");
  //
  tO2vtx->Branch("fChi2", &vtx.fChi2, "fChi2/F");
  tO2vtx->Branch("fN", &vtx.fN, "fN/i");
  //
  tO2vtx->Branch("fOrbitNumber", &vtx.fOrbitNumber, "fOrbitNumber/i");
  tO2vtx->Branch("fPeriodNumber", &vtx.fPeriodNumber, "fPeriodNumber/i");
  tO2vtx->Branch("fBCNumber", &vtx.fBCNumber, "fBCNumber/s");
  
  //---------------------------------------------------------------------------
  // Barrel tracks
  //---------------------------------------------------------------------------
  TTree * tO2tracks = new TTree("O2tracks","Barrel tracks");

  // Data members according to the table of Ruben

  // Identifier to associate tracks to collisions.
  // The inversed association (collisions to tracks) is generated on the fly
  struct {
    Int_t             fID4Tracks;  // The index of the collision vertex, to which the track is attached
    UChar_t           fDIDbefore;  // The difference of indexes for the first preceeding collision the track can be attached to
    UChar_t           fDIDafter;   // The difference of indexes for the last succeeding collision the track can be attached to
    // Coordinate system parameters
    Float_t           fX;         // X coordinate for the point of parametrisation
    Float_t           fAlpha;     // Local <--> global coor.system rotation angle
    // Track parameters
    Float_t           fY;         // fP[0] local Y-coordinate of a track (cm) 
    Float_t           fZ;         // fP[1] local Z-coordinate of a track (cm)
    Float_t           fSnp;       // fP[2] local sine of the track momentum azimuthal angle
    Float_t           fTgl;       // fP[3] tangent of the track momentum dip angle
    Float_t           fSigned1Pt; // fP[4] 1/pt (1/(GeV/c))
    Int_t             fPID;
    // Covariance matrix
    Float_t           fCYY;       // fC[0]
    Float_t           fCZY;       // fC[1]
    Float_t           fCZZ;       // fC[2]
    Float_t           fCSnpY;     // fC[3]
    Float_t           fCSnpZ;     // fC[4]
    Float_t           fCSnpSnp;   // fC[5]
    Float_t           fCTglY;     // fC[6]
    Float_t           fCTglZ;     // fC[7]
    Float_t           fCTglSnp;   // fC[8]
    Float_t           fCTglTgl;   // fC[9]
    Float_t           fC1PtY;     // fC[10]
    Float_t           fC1PtZ;     // fC[11]
    Float_t           fC1PtSnp;   // fC[12]
    Float_t           fC1PtTgl;   // fC[13]
    Float_t           fC1Pt21Pt2; // fC[14]
    // Additional track parameters
    Float_t           fTPCinnerP; // Full momentum at the inner wall of TPC for dE/dx PID
    // Track quality parameters
    ULong64_t         fFlags;     // Reconstruction status flags
    // Clusters
    UChar_t           fITSClusterMap;  // ITS map of clusters, one bit per a layer
    UShort_t          fTPCncls;        // number of clusters assigned in the TPC
    UChar_t           fTRDntracklets;  // number of TRD tracklets used for tracking/PID (TRD/TOF pattern)	
    // Chi2
    Float_t           fITSchi2Ncl;     // chi2/Ncl ITS	
    Float_t           fTPCchi2Ncl;     // chi2/Ncl TPC		
    Float_t           fTRDchi2;  // chi2 TRD match (?)	
    Float_t           fTOFchi2;  // chi2 TOF match (?)	
    // PID
    Float_t           fTPCsignal;// dE/dX TPC 		
    Float_t           fTRDsignal;// dE/dX TRD		
    Float_t           fTOFsignal;// TOFsignal		
    Float_t           fLength;   // Int.Lenght @ TOF	
  } tracks;
  // Branches
  tO2tracks->Branch("fID4Tracks", &tracks.fID4Tracks, "fID4Tracks/I");
  tO2tracks->Branch("fDIDbefore", &tracks.fDIDbefore, "fDIDbefore/b");
  tO2tracks->Branch("fDIDafter", &tracks.fDIDafter, "fDIDafter/b");
  //
  tO2tracks->Branch("fX",&tracks.fX,"fX/F");
  tO2tracks->Branch("fAlpha", &tracks.fAlpha, "fAlpha/F");
  tO2tracks->Branch("fY", &tracks.fY, "fY/F");
  tO2tracks->Branch("fZ", &tracks.fZ, "fZ/F");
  tO2tracks->Branch("fSnp", &tracks.fSnp, "fSnp/F");
  tO2tracks->Branch("fTgl", &tracks.fTgl, "fTgl/F");
  tO2tracks->Branch("fSigned1Pt", &tracks.fSigned1Pt, "fSigned1Pt/F");
  //
  tO2tracks->Branch("fPID", &tracks.fPID, "fPID/I");
  //
  tO2tracks->Branch("fCYY", &tracks.fCYY, "fCYY/F");
  tO2tracks->Branch("fCZY", &tracks.fCZY, "fCZY/F");
  tO2tracks->Branch("fCZZ", &tracks.fCZZ, "fCZZ/F");
  tO2tracks->Branch("fCSnpY", &tracks.fCSnpY, "fCSnpY/F");
  tO2tracks->Branch("fCSnpZ", &tracks.fCSnpZ, "fCSnpZ/F");
  tO2tracks->Branch("fCSnpSnp", &tracks.fCSnpSnp, "fCSnpSnp/F");
  tO2tracks->Branch("fCTglY", &tracks.fCTglY, "fCTglY/F");
  tO2tracks->Branch("fCTglZ", &tracks.fCTglZ, "fCTglZ/F");
  tO2tracks->Branch("fCTglSnp", &tracks.fCTglSnp, "fCTglSnp/F");
  tO2tracks->Branch("fCTglTgl", &tracks.fCTglTgl, "fCTglTgl/F");
  tO2tracks->Branch("fC1PtY", &tracks.fC1PtY, "fC1PtY/F");
  tO2tracks->Branch("fC1PtZ", &tracks.fC1PtZ, "fC1PtZ/F");
  tO2tracks->Branch("fC1PtSnp", &tracks.fC1PtSnp, "fC1PtSnp/F");
  tO2tracks->Branch("fC1PtTgl", &tracks.fC1PtTgl, "fC1PtTgl/F");
  tO2tracks->Branch("fC1Pt21Pt2", &tracks.fC1Pt21Pt2, "fC1Pt21Pt2/F");
  //
  tO2tracks->Branch("fTPCinnerP", &tracks.fTPCinnerP, "fTPCinnerP/F");
  //
  tO2tracks->Branch("fFlags", &tracks.fFlags, "fFlags/l");
  //
  tO2tracks->Branch("fITSClusterMap", &tracks.fITSClusterMap, "fITSClusterMap/b");
  tO2tracks->Branch("fTPCncls", &tracks.fTPCncls, "fTPCncls/s");
  tO2tracks->Branch("fTRDntracklets", &tracks.fTRDntracklets, "fTRDntracklets/b");
  //
  tO2tracks->Branch("fITSchi2Ncl", &tracks.fITSchi2Ncl, "fITSchi2Ncl/F");
  tO2tracks->Branch("fTPCchi2Ncl", &tracks.fTPCchi2Ncl, "fTPCchi2Ncl/F");
  tO2tracks->Branch("fTRDchi2", &tracks.fTRDchi2, "fTRDchi2/F");
  tO2tracks->Branch("fTOFchi2", &tracks.fTOFchi2, "fTOFchi2/F");
  //
  tO2tracks->Branch("fTPCsignal", &tracks.fTPCsignal, "fTPCsignal/F");
  tO2tracks->Branch("fTRDsignal", &tracks.fTRDsignal, "fTRDsignal/F");
  tO2tracks->Branch("fTOFsignal", &tracks.fTOFsignal, "fTOFsignal/F");
  tO2tracks->Branch("fLength", &tracks.fLength, "fLength/F");
  //---------------------------------------------------------------------------

  // Calorimeter cells
  TTree * tO2calo = new TTree("O2calo","Calorimeter cells");
  // Data members according to https://alice.its.cern.ch/jira/browse/O2-423
  // No compression is used so far
  struct {
    Int_t             fID4Calo;     // The index of the collision vertex, to which the cell is attached
    Short_t           fCellNumber; // Cell absolute Id. number
    Float_t           fAmplitude;  // Cell amplitude (= energy!)
    Float_t           fTime;       // Cell time
    Char_t            fType;       // Cell type (-1 is undefined, 0 is PHOS, 1 is EMCAL)
  } calo;
  // Branches
  tO2calo->Branch("fID4Calo", &calo.fID4Calo, "fID4Calo/I");
  //
  tO2calo->Branch("fCellNumber", &calo.fCellNumber, "fCellNumber/S");
  tO2calo->Branch("fAmplitude", &calo.fAmplitude, "fAmplitude/F");
  tO2calo->Branch("fTime", &calo.fTime, "fTime/F");
  tO2calo->Branch("fType", &calo.fType, "fType/B");
  //---------------------------------------------------------------------------
  // MUON tracks
  //---------------------------------------------------------------------------
  TTree * tO2mu = new TTree("O2mu","Muon tracks");
  struct {
    Int_t             fID4mu;     // The index of the collision vertex, to which the muon is attached
    // parameters at vertex
    Float_t fInverseBendingMomentum; ///< Inverse bending momentum (GeV/c ** -1) times the charge 
    Float_t fThetaX;                 ///< Angle of track at vertex in X direction (rad)
    Float_t fThetaY;                 ///< Angle of track at vertex in Y direction (rad)
    Float_t fZ;                      ///< Z coordinate (cm)
    Float_t fBendingCoor;            ///< bending coordinate (cm)
    Float_t fNonBendingCoor;         ///< non bending coordinate (cm)
    /// reduced covariance matrix of UNCORRECTED track parameters, ordered as follow:      <pre>
    /// [0] =  <X,X>
    /// [1] =<X,ThetaX>  [2] =<ThetaX,ThetaX>
    /// [3] =  <X,Y>     [4] =  <Y,ThetaX>     [5] =  <Y,Y>
    /// [6] =<X,ThetaY>  [7] =<ThetaX,ThetaY>  [8] =<Y,ThetaY>  [9] =<ThetaY,ThetaY>
    /// [10]=<X,InvP_yz> [11]=<ThetaX,InvP_yz> [12]=<Y,InvP_yz> [13]=<ThetaY,InvP_yz> [14]=<InvP_yz,InvP_yz>  </pre>
    Float_t fCovariances[15]; ///< \brief reduced covariance matrix of parameters AT FIRST CHAMBER
    // global tracking info
    Float_t fChi2;                ///< chi2 in the MUON track fit
    Float_t fChi2MatchTrigger;    ///< chi2 of trigger/track matching
    //
  } muons;
  tO2mu->Branch("fID4mu", &muons.fID4mu, "fID4mu/I");
  tO2mu->Branch("fInverseBendingMomentum", &muons.fInverseBendingMomentum, "fInverseBendingMomentum/F");
  tO2mu->Branch("fThetaX", &muons.fThetaX, "fThetaX/F");
  tO2mu->Branch("fThetaY", &muons.fThetaY, "fThetaY/F");
  tO2mu->Branch("fZ", &muons.fZ, "fZ/F");
  tO2mu->Branch("fBendingCoor", &muons.fBendingCoor, "fBendingCoor/F");
  tO2mu->Branch("fNonBendingCoor", &muons.fNonBendingCoor, "fNonBendingCoor/F");
  //
  tO2mu->Branch("fCovariances", muons.fCovariances, "fCovariances[15]/F");
  //
  tO2mu->Branch("fChi2", &muons.fChi2, "fChi2/F");
  tO2mu->Branch("fChi2MatchTrigger", &muons.fChi2MatchTrigger, "fChi2MatchTrigger/F");
  //---------------------------------------------------------------------------
  // Muon clisters
  //---------------------------------------------------------------------------
  // VZERO as proxy for FIT
  TTree * tO2vz = new TTree("O2vz","VZERO");
  struct {
    Int_t             fID4vz;     // The index of the collision vertex, to which the muon is attached
    //
    Float_t fAdc[64];          //  adc for each channel
    Float_t fTime[64];         //  time for each channel
    Float_t fWidth[64];        //  time width for each channel
  } vzero;
  //
  tO2mu->Branch("fID4vz", &vzero.fID4vz, "fID4vz/I");
  tO2vz->Branch("fAdc", vzero.fAdc, "fAdc[64]/F");
  tO2vz->Branch("fTime", vzero.fTime, "fTime[64]/F");
  tO2vz->Branch("fWidth", vzero.fWidth, "fWidth[64]/F");
  //---------------------------------------------------------------------------
  UInt_t compressionMask = 0xFFD0;
  //---------------------------------------------------------------------------

  // Loop on events
  Long64_t nev = tEsd->GetEntries();

  cout << nev << " ESD events" << endl;

  for (Long64_t iev=0; iev<nev; iev++) {
	cout << "Event: " << iev << endl;
	
    // Reset the current ESD event
    esd->Reset();
    
    // Get the next ESD
    tEsd->GetEntry(iev);

    // Needed by TOF
    esd->ConnectTracks();
    
    //---------------------------------------------------------------------------
    // Store the vertex

    const AliESDVertex * vert = esd->GetPrimaryVertex();
    //
    vtx.fX            = float2half(vert->GetX());
    vtx.fY            = float2half(vert->GetY());
    vtx.fZ            = float2half(vert->GetZ());
    //
    Double_t cov[6] = {0};
    vert->GetCovarianceMatrix(cov);
    vtx.fCovXX        = float2half(cov[0]);
    vtx.fCovXY        = float2half(cov[1]);
    vtx.fCovXZ        = float2half(cov[3]);
    vtx.fCovYY        = float2half(cov[2]);
    vtx.fCovYZ        = float2half(cov[4]);
    vtx.fCovZZ        = float2half(cov[5]);
    //
    vtx.fChi2         = float2half(vert->GetChi2());
    vtx.fN            = vert->GetNContributors();
    //
    vtx.fOrbitNumber  = esd->GetOrbitNumber();
    vtx.fPeriodNumber = esd->GetPeriodNumber();
    vtx.fBCNumber     = esd->GetBunchCrossNumber();

    tO2vtx->Fill();
    
    //---------------------------------------------------------------------------
    // Loop on tracks
    Int_t ntrk = esd->GetNumberOfTracks();
    // Note: some events have zero reconstructed tracks
    
    tracks.fID4Tracks  = iev; // The ID of the current event to which the track belong
    tracks.fDIDbefore  = gRandom->Integer(3);
    tracks.fDIDafter   = gRandom->Integer(3);
    for(Int_t itrk=0; itrk<ntrk; itrk++) {
      AliESDtrack * track = esd->GetTrack(itrk);
      track->SetESDEvent(esd); // Needed by TOF (?)

      //
      tracks.fX         = float2half(track->GetX());
      tracks.fAlpha     = float2half(track->GetAlpha());
      //
      tracks.fY         = float2half(track->GetY());
      tracks.fZ         = float2half(track->GetZ());
      tracks.fSnp       = float2half(track->GetSnp());
      tracks.fTgl       = float2half(track->GetTgl());
      tracks.fSigned1Pt = float2half(track->GetSigned1Pt());
      //
      tracks.fPID = track->GetPID();
      //
      tracks.fCYY       = float2half(track->GetSigmaY2());
      tracks.fCZY       = float2half(track->GetSigmaZY());
      tracks.fCZZ       = float2half(track->GetSigmaZ2());
      tracks.fCSnpY     = float2half(track->GetSigmaSnpY());
      tracks.fCSnpZ     = float2half(track->GetSigmaSnpZ());
      tracks.fCSnpSnp   = float2half(track->GetSigmaSnp2());
      tracks.fCTglY     = float2half(track->GetSigmaTglY());
      tracks.fCTglZ     = float2half(track->GetSigmaTglZ());
      tracks.fCTglSnp   = float2half(track->GetSigmaTglSnp());
      tracks.fCTglTgl   = float2half(track->GetSigmaTgl2());
      tracks.fC1PtY     = float2half(track->GetSigma1PtY());
      tracks.fC1PtZ     = float2half(track->GetSigma1PtZ());
      tracks.fC1PtSnp   = float2half(track->GetSigma1PtSnp());
      tracks.fC1PtTgl   = float2half(track->GetSigma1PtTgl());
      tracks.fC1Pt21Pt2 = float2half(track->GetSigma1Pt2());
      //
      const AliExternalTrackParam * intp = track->GetTPCInnerParam();
      tracks.fTPCinnerP = float2half((intp?intp->GetP():0)); // Set the momentum to 0 if the track did not reach TPC
      //
      tracks.fFlags     = track->GetStatus();
      //
      tracks.fITSClusterMap = track->GetITSClusterMap();
      tracks.fTPCncls = track->GetTPCNcls();
      tracks.fTRDntracklets = track->GetTRDntracklets();
      //
      tracks.fITSchi2Ncl = float2half((track->GetITSNcls()?track->GetITSchi2()/track->GetITSNcls():0));
      tracks.fTPCchi2Ncl = float2half((track->GetTPCNcls()?track->GetTPCchi2()/track->GetTPCNcls():0));
      tracks.fTRDchi2 = float2half(track->GetTRDchi2());
      tracks.fTOFchi2 = float2half(track->GetTOFchi2());
      //
      tracks.fTPCsignal = float2half(track->GetTPCsignal());
      tracks.fTRDsignal = float2half(track->GetTRDsignal());
      tracks.fTOFsignal = float2half(track->GetTOFsignal());
      tracks.fLength = float2half(track->GetIntegratedLength());
      
      // Store the exctracted information
      tO2tracks->Fill(); // Fill the current track
    } // End loop on tracks
    //---------------------------------------------------------------------------
    // Calorimeters
    // Loop on EMCAL cells
    AliESDCaloCells * cells = esd->GetEMCALCells();
    Short_t nCells = cells->GetNumberOfCells();
    calo.fID4Calo = iev;
    for (Short_t ice = 0; ice<nCells; ++ice) {
      Short_t cellNumber;
      Double_t amplitude; 
      Double_t time;
      Int_t mclabel;
      Double_t efrac;
      
      cells->GetCell(ice, cellNumber, amplitude, time, mclabel, efrac);
      calo.fCellNumber = cellNumber;
      calo.fAmplitude = float2half(amplitude);
      calo.fTime = float2half(time);
      calo.fType = cells->GetType(); // common for all cells
      // Store the exctracted information
      tO2calo->Fill(); // Fill the current calorimeter cell
    } // End loop on EMCAL cells

    // Loop on PHOS cells
    cells = esd->GetPHOSCells();
    nCells = cells->GetNumberOfCells();
    for (Short_t icp = 0; icp<nCells; ++icp) {
      Short_t cellNumber;
      Double_t amplitude; 
      Double_t time;
      Int_t mclabel;
      Double_t efrac;
      
      cells->GetCell(icp, cellNumber, amplitude, time, mclabel, efrac);
      calo.fCellNumber = cellNumber;
      calo.fAmplitude = float2half(amplitude);
      calo.fTime = float2half(time);
      calo.fType = cells->GetType(); // common for all cells
      // Store the exctracted information
      tO2calo->Fill(); // Fill the current calorimeter cell
    } // End loop on PHOS cells
    //---------------------------------------------------------------------------
    // Muon tracks
    muons.fID4mu  = iev; // The ID of the current event to which the track belong
    Int_t nmu = esd->GetNumberOfMuonTracks();
    for (Int_t imu=0; imu<nmu; ++imu) {
      AliESDMuonTrack* mutrk = esd->GetMuonTrack(imu);
      //
      muons.fInverseBendingMomentum = float2half(mutrk->GetInverseBendingMomentum());
      muons.fThetaX = float2half(mutrk->GetThetaX());
      muons.fThetaY = float2half(mutrk->GetThetaY());
      muons.fZ = float2half(mutrk->GetZ());
      muons.fBendingCoor = float2half(mutrk->GetBendingCoor());
      muons.fNonBendingCoor = float2half(mutrk->GetNonBendingCoor());
      //
      TMatrixD covm;
      mutrk->GetCovariances(covm);
      for (Int_t i = 0; i < 5; i++)
	    for (Int_t j = 0; j <= i; j++)
	      muons.fCovariances[i*(i+1)/2 + j] = float2half(covm(i,j));
      //
      muons.fChi2 = float2half(mutrk->GetChi2());
      muons.fChi2MatchTrigger = float2half(mutrk->GetChi2MatchTrigger());
      //
      tO2mu->Fill();
    } // End loop on muon tracks
    //---------------------------------------------------------------------------
    // VZERO
    AliESDVZERO * vz = esd->GetVZEROData();
    vzero.fID4vz  = iev; // The ID of the current event to which the VZERO data belong
    for (Int_t ich=0; ich<64; ++ich) {
      vzero.fAdc[ich] = float2half(vz->GetAdc(ich));
      vzero.fTime[ich] = float2half(vz->GetTime(ich));
      vzero.fWidth[ich] = float2half(vz->GetWidth(ich));
    }
    tO2vz->Fill();
    //---------------------------------------------------------------------------
  }

  // Write the resulting O2 trees and close the file
  tO2vtx->Write();
  tO2tracks->Write();
  tO2calo->Write();
  tO2mu->Write();
  tO2vz->Write();
  fO2->Close();

}

int main(int argc, char** argv) {

  TStopwatch timer;
  timer.Start();

  // Make a chain and convert the original ESD to the Run3 format

  TChain * chain = new TChain("esdTree");
  chain->Add("AliESDs.root");
//  chain->Add("15000244918034.508.AliESDs.root");
//  chain->Add("15000244918034.509.AliESDs.root");
//  chain->Add("15000244918034.510.AliESDs.root");
  
  convertO2(chain);

  delete chain;
  
  timer.Stop();
  timer.Print();

  return 0;
}

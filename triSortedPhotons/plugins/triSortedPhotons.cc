// -*- C++ -*-
//
// Package:    triphoton/triSortedPhotons
// Class:      triSortedPhotons
// 
/**\class triSortedPhotons triSortedPhotons.cc triphoton/triSortedPhotons/plugins/triSortedPhotons.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Cilicia Uzziel Perez
//         Created:  Thu, 15 Jun 2017 23:54:47 GMT
//
//


// system include files
#include <memory>
#include <iostream>
#include <fstream>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Photon.h"

#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class triSortedPhotons : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit triSortedPhotons(const edm::ParameterSet&);
      ~triSortedPhotons();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      
      edm::EDGetToken photonsMiniAODToken_;
      TTree *fTree;
      struct eventInfo_t {
        Long64_t run;
        Long64_t LS;
        Long64_t evnum;
      };
      eventInfo_t fEventInfo;

      struct photonInfo_t {
        double pt;
        double eta; 
        double phi;
        double sceta;
        double scphi;
      };
      //Instantiate the different photon structs to store photoninfo
      photonInfo_t photon1;
      photonInfo_t photon2; 
      photonInfo_t photon3; 

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
triSortedPhotons::triSortedPhotons(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs; 
      fTree = fs->make<TTree>("fTree","TriphotonTree");
      fTree->Branch("Event",&fEventInfo, "run/L:LS:evnum");
      fTree->Branch("Photon1", &photon1, "pt/F:eta:phi:sceta:scphi");
      fTree->Branch("Photon2", &photon2, "pt/F:eta:phi:sceta:scphi");
      fTree->Branch("Photon3", &photon3, "pt/F:eta:phi:sceta:scphi");
 
   photonsMiniAODToken_ = mayConsume<edm::View<pat::Photon>>(edm::InputTag("slimmedPhotons"));
}


triSortedPhotons::~triSortedPhotons()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
triSortedPhotons::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std; 
 
   //Initialize leaf values 
   fEventInfo.run = -999;
   fEventInfo.LS = -999;
   fEventInfo.evnum = -999;
   
  photon1.pt = -999;
  photon1.eta = -999;
  photon1.phi = -999;
  photon2.pt = -999;
  photon2.eta = -999;
  photon2.phi = -999;
  photon3.pt = -999;
  photon3.eta = -999;
  photon3.phi = -999;
  photon1.sceta = -999;
  photon1.scphi = -999;
  photon2.sceta = -999;
  photon2.scphi = -999;
  photon3.sceta = -999;
  photon3.scphi = -999;  
  
  //Initialize variables for boolean sorting 
  double pholead1 = 0; 
  double pholead2 = 0; 
  double pholead3 = 0;

  //Store Print out in a file: 
  ofstream cout("photon.txt", ios::app);

  //Print out run information 
  cout << "Run: " << iEvent.id().run() << ", LS: " << iEvent.id().luminosityBlock() << ", Event: " << iEvent.id().event() << endl;
  
  //Store Run info
  fEventInfo.run = iEvent.id().run();
  fEventInfo.LS = iEvent.id().luminosityBlock();
  fEventInfo.evnum = iEvent.id().event();

  edm::Handle<edm::View<pat::Photon>> photons;   //To get photon collection
  iEvent.getByToken(photonsMiniAODToken_,photons);

  //Loop over each photon in each event
  //for (edm::View<pat::Photon>::const_iterator pho = photons->begin(); pho != photons->end(); ++pho) {
  for (size_t i = 0; i<photons->size(); ++i) {
  
    const auto pho = photons->ptrAt(i);
  
  //print out all Photon information
  cout << "Photon: " << "pt = " << pho->pt() << "; eta = " << pho->eta() <<"; phi = " <<pho->phi() << "; sceta = " 
    << pho->superCluster()->eta() << "; scphi = " << pho->superCluster()->phi() << endl;
     
    //LOGICAL SORTING EXERCISE: Although the Photon pts are already sorted, we're going to make sure that they are. 
    //Here we write a boolean sorting algorithm and we'll print out the results for comparison. We'll also store the photon
    //info herein. 

   if (pho->pt()>= pholead1){
     pholead3 = pholead2; pholead2=pholead1; pholead1= pho->pt();
     photon1.pt = pholead1;  photon1.eta = pho->eta(); photon1.phi = pho->phi();  
     photon1.sceta = pho->superCluster()->eta(); photon1.scphi = pho->superCluster()->phi();  
   }
  else if (pho->pt() > pholead2){
    pholead3 = pholead2; pholead2 = pho->pt();
    photon2.pt = pholead2; photon2.eta = pho->eta(); photon2.phi = pho->phi();  
    photon2.sceta = pho->superCluster()->eta(); photon2.scphi = pho->superCluster()->phi();  
  }
   else if (pho->pt() > pholead3){
    pholead3 = pho->pt();
    photon3.pt = pholead3; photon3.eta = pho->eta(); photon3.phi = pho->phi();  
    photon3.sceta = pho->superCluster()->eta(); photon3.scphi = pho->superCluster()->phi();  
   }
   
  }//end of photon loop 

  //****** CHECKING IF WE ARE STORING THE RIGHT INFO **********
  int phosize = photons->size();
  cout <<"SIZE: "<< phosize  << "; lead_pt: " << pholead1 << "; sub_leading_photon: " << pholead2 << "; sub_sub_leading: " << pholead3 << endl;
  cout <<"CHECK INFO:" << endl;
  //Photon1
  cout << "pt: " << photon1.pt << "; eta: " << photon1.eta << "; phi: " << photon1.phi 
    << "; sceta: " << photon1.sceta<< "; scphi: " << photon1.scphi <<endl;  
  //Photon2
  cout << "pt: " << photon2.pt << "; eta: " << photon2.eta << "; phi: " << photon2.phi 
    << "; sceta: " << photon2.sceta<< "; scphi: " << photon2.scphi <<endl;    
  //Photon3
  cout << "pt: " << photon3.pt << "; eta: " << photon3.eta << "; phi: " << photon3.phi 
    << "; sceta: " << photon3.sceta<< "; scphi: " << photon3.scphi <<endl;  
  cout << "--------------------------------------------- RUN ENDS -----------------------------------------"<<endl;
//We only fill tree for events with at least three photons: 
   if (phosize >2) fTree->Fill();



#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
triSortedPhotons::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
triSortedPhotons::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
triSortedPhotons::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(triSortedPhotons);

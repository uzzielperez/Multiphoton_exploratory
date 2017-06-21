// -*- C++ -*-
//
// Package:    triphoton/objtriphotons
// Class:      objtriphotons
// 
/**\class objtriphotons objtriphotons.cc triphoton/objtriphotons/plugins/objtriphotons.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Cilicia Uzziel Perez
//         Created:  Wed, 21 Jun 2017 01:32:32 GMT
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

class objtriphotons : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit objtriphotons(const edm::ParameterSet&);
      ~objtriphotons();

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

    struct photonInfo_t{
      double pt; 
      double eta; 
      double phi;
      double sceta; 
      double scphi;
    };

    //Instantiate the different photon structs to store photoninfo
    photonInfo_t photon1Info;
    photonInfo_t photon2Info;
    photonInfo_t photon3Info;
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
objtriphotons::objtriphotons(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs;
   fTree = fs->make<TTree>("fTree","TriphotonTree");
   fTree->Branch("Event",&fEventInfo, "run/L:LS:evnum");
   fTree->Branch("Photon1", &photon1Info, "pt/F:eta:phi:sceta:scphi");
   fTree->Branch("Photon2", &photon2Info, "pt/F:eta:phi:sceta:scphi");
   fTree->Branch("Photon3", &photon3Info, "pt/F:eta:phi:sceta:scphi");

   photonsMiniAODToken_ = mayConsume<edm::View<pat::Photon>>(edm::InputTag("slimmedPhotons"));

}


objtriphotons::~objtriphotons()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
objtriphotons::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std; 

   //Initialize leaf values
   fEventInfo.run = -999;
   fEventInfo.LS = -999;
   fEventInfo.evnum = -999;

   photon1Info.pt = -999;
   photon1Info.eta = -999;
   photon1Info.phi = -999;
   photon1Info.sceta = -999;
   photon1Info.scphi = -999;
   
   photon2Info.pt = -999;
   photon2Info.eta = -999;
   photon2Info.phi = -999;
   photon2Info.sceta = -999;
   photon2Info.scphi = -999;
   
   photon3Info.pt = -999;
   photon3Info.eta = -999;
   photon3Info.phi = -999;
   photon3Info.sceta = -999;
   photon3Info.scphi = -999;

   //Initialize variables for boolean sorting 
  double pholead1 = 0; 
  double pholead2 = 0; 
  double pholead3 = 0;

  //Store Print out in a file: 
  ofstream cout("photon.txt", ios::app);

  //Print out event information 
  cout << "Run: " << iEvent.id().run() << ", LS: " << iEvent.id().luminosityBlock() 
    << ", Event: " << iEvent.id().event() << endl;
  
  //Store event info
  fEventInfo.run = iEvent.id().run();
  fEventInfo.LS = iEvent.id().luminosityBlock();
  fEventInfo.evnum = iEvent.id().event();
  
  edm::Handle<edm::View<pat::Photon>> photons;   //To get photon collection
  iEvent.getByToken(photonsMiniAODToken_,photons);
  
  const pat::Photon *photon1 = NULL;
  const pat::Photon *photon2 = NULL;
  const pat::Photon *photon3 = NULL; 

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
        // store appropriate photon object
        photon1 = &(*pho);
    }
    else if (pho->pt() > pholead2){
         pholead3 = pholead2; pholead2 = pho->pt();
         //store appropriate photon object
         photon2 = &(*pho);
    }
    else if (pho->pt() > pholead3){
         pholead3 = pho->pt();
         photon3 = &(*pho);
    }

  }//End of photon loop

 //****** CHECKING IF WE ARE STORING THE RIGHT INFO **********
 
  int phosize = photons->size();
  if (phosize >2){
    cout <<"SIZE: "<< phosize  << "; lead_pt: " << pholead1 << "; sub_leading_photon: " 
      << pholead2 << "; sub_sub_leading: " << pholead3 << endl;
    cout <<"FromOBJ: "<< phosize  << "; lead_pt: " << photon1->pt() << "; sub_leading_photon: "
      << photon2->pt() << "; sub_sub_leading: " << photon3->pt() << endl;
  
    //Fill struct info
   //PHOTON1:
   photon1Info.pt = photon1->pt(); 
   photon1Info.eta = photon1->eta();
   photon1Info.phi = photon1->phi();
   photon1Info.sceta = photon1->superCluster()->eta();
   photon1Info.scphi = photon1->superCluster()->phi();
   
   //PHOTON2:
   photon2Info.pt = photon2->pt(); 
   photon2Info.eta = photon2->eta();
   photon2Info.phi = photon2->phi();
   photon2Info.sceta = photon2->superCluster()->eta();
   photon2Info.scphi = photon2->superCluster()->phi();

    //PHOTON3:
   photon3Info.pt = photon3->pt(); 
   photon3Info.eta = photon3->eta();
   photon3Info.phi = photon3->phi();
   photon3Info.sceta = photon3->superCluster()->eta();
   photon3Info.scphi = photon3->superCluster()->phi();
   
   cout <<"CHECK INFO:" << endl;
     //Photon1
    cout << "pt: " << photon1Info.pt << "; eta: " << photon1Info.eta << "; phi: " << photon1Info.phi 
         << "; sceta: " << photon1Info.sceta<< "; scphi: " << photon1Info.scphi <<endl;  
    //Photon2
    cout << "pt: " << photon2Info.pt << "; eta: " << photon2Info.eta << "; phi: " << photon2Info.phi 
         << "; sceta: " << photon2Info.sceta<< "; scphi: " << photon2Info.scphi <<endl;    
    //Photon3
    cout << "pt: " << photon3Info.pt << "; eta: " << photon3Info.eta << "; phi: " << photon3Info.phi 
         << "; sceta: " << photon3Info.sceta<< "; scphi: " << photon3Info.scphi <<endl;  
   
  }
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
objtriphotons::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
objtriphotons::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
objtriphotons::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(objtriphotons);

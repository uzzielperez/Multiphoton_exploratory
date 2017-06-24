// -*- C++ -*-
//
// Package:    triphoton/triobj
// Class:      triobj
// 
/**\class triobj triobj.cc triphoton/triobj/plugins/triobj.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Cilicia Uzziel Perez
//         Created:  Fri, 23 Jun 2017 08:07:06 GMT
//
//


// system include files
#include <memory>
#include <iostream>
#include <fstream>
#include <vector>

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

class triobj : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit triobj(const edm::ParameterSet&);
      ~triobj();

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
    photonInfo_t photonInfo;

    struct photonPtComparer{
      bool operator()(const photonInfo_t& x, const photonInfo_t& y)const{
        return x.pt> y.pt; 
      }
    };

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
triobj::triobj(const edm::ParameterSet& iConfig)

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


triobj::~triobj()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
triobj::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  
  /*const pat::Photon *photon1 = NULL;
  const pat::Photon *photon2 = NULL;
  const pat::Photon *photon3 = NULL; 
  */

  //Triphoton object container 
  vector<photonInfo_t> triphoton_obj;
  //vector<edm::Ptr<pat::Photon>> triphotons;

  //Loop over each photon in each event
  //for (edm::View<pat::Photon>::const_iterator pho = photons->begin(); pho != photons->end(); ++pho) {
  for (size_t i = 0; i<photons->size(); ++i) {
    const auto pho = photons->ptrAt(i);
    //print out all Photon information
    cout << "Photon: " << "pt = " << pho->pt() << "; eta = " << pho->eta() <<"; phi = " <<pho->phi() << "; sceta = " 
         << pho->superCluster()->eta() << "; scphi = " << pho->superCluster()->phi() << endl;
     
    //If there are triphotons, we store the objects in the vector for sorting 
    photonInfo.pt = pho->pt();
    photonInfo.eta = pho->eta();
    photonInfo.phi = pho->phi();
    photonInfo.sceta = pho->superCluster()->eta();
    photonInfo.scphi = pho->superCluster()->phi();
    if (photons->size()>2){
      triphoton_obj.push_back(photonInfo); 
    } 
   
    //Sort HERE. See photonPtComparer() struct above  
   sort(triphoton_obj.rbegin(), triphoton_obj.rend(), photonPtComparer());

}//End of photon loop

//***********CHECK if the pts are sorted in ascending order******************
if (photons->size()>2) cout << "SORTED PHOTON pt" << endl;
for (auto i = triphoton_obj.begin(); i != triphoton_obj.end(); i++) cout << i->pt <<endl;

   

//*************WORK IN PROGRESS******************
// Here I am trying to access each object in the vector and then assigning them to photon1Info, 
// photon2Info, photon3Info. IT DOES NOT WORK!!!! 
//Fill struct info
/*for (auto i = triphoton_obj.begin(); i != triphoton_obj.end(); i++){
    if (i==0) photon1Info.pt = i->pt();
    if (i==2) photon2Info.pt = i->pt();
    else if (i==3) photon3Info.pt = i->pt();   
 }*/

//cout << (*triphoton_obj)[0].pt<<endl; //pointer 
//cout << triphoton_obj[1].pt <<endl;
//cout << triphoton_obj[2].pt <<endl;

   //We only fill tree for events with at least three photons: 

   if (photons->size()>2) fTree->Fill();

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
triobj::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
triobj::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
triobj::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(triobj);

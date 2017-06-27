// -*- C++ -*-
//
// Package:    triphoton/threephoton
// Class:      threephoton
// 
/**\class threephoton threephoton.cc triphoton/threephoton/plugins/threephoton.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Cilicia Uzziel Perez
//         Created:  Mon, 26 Jun 2017 03:32:17 GMT
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

class threephoton : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit threephoton(const edm::ParameterSet&);
      ~threephoton();

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
    photonInfo_t photonInfo;
    photonInfo_t triphoton_info[3]; //one for each photon
    
    static bool comparePhotonsByPt(const edm::Ptr<pat::Photon> photon1, const edm::Ptr<pat::Photon> photon2) {
          return(photon1->pt()>=photon2->pt());
            }


    struct photonPtComparer{
     bool operator()(const photonInfo_t& x, const photonInfo_t& y)const{
        return x.pt<y.pt; //for sorting in descending order 
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
threephoton::threephoton(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   usesResource("TFileService");
 edm::Service<TFileService> fs;
   fTree = fs->make<TTree>("fTree","TriphotonTree");
   fTree->Branch("Event",   &fEventInfo, "run/L:LS:evnum");
   fTree->Branch("Photon1", &triphoton_info[0], "pt/D:eta:phi:sceta:scphi");
   fTree->Branch("Photon2", &triphoton_info[1], "pt/D:eta:phi:sceta:scphi");
   fTree->Branch("Photon3", &triphoton_info[2], "pt/D:eta:phi:sceta:scphi");

   photonsMiniAODToken_ = mayConsume<edm::View<pat::Photon>>(edm::InputTag("slimmedPhotons"));


}


threephoton::~threephoton()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
threephoton::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;  
//Initialize leaf values
   fEventInfo.run = -999;
   fEventInfo.LS = -999;
   fEventInfo.evnum = -999;

  for(int k = 0; k<3; k++){
  triphoton_info[k] = (photonInfo_t){-999,-999,-999,-999,-999};
  }
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

  //Triphoton object container 
  vector<photonInfo_t> triphoton_obj;
  vector<edm::Ptr<pat::Photon>> photons;

//  vector<pat::Photon*> photons; 

  //Loop over each photon in each event
  //for (edm::View<pat::Photon>::const_iterator pho = photons->begin(); pho != photons->end(); ++pho) {
  for (size_t i = 0; i<photons->size(); ++i) {
    const auto pho = photons->ptrAt(i);
    //print out all Photon information
    cout << "Photon: " << "pt = " << pho->pt() << "; eta = " << pho->eta() <<"; phi = " <<pho->phi() << "; sceta = " 
         << pho->superCluster()->eta() << "; scphi = " << pho->superCluster()->phi() << endl;
    
   // if(photons->size()>2){
         //triphotons.push_back(&(*pho));
          photons.push_back(pho);
    //}

    //Sort HERE. See photonPtComparer() struct above  
   //sort(photon_obj.rbegin(), photon_obj.rend(), photonPtComparer());
   sort(photons.rbegin(), photons.rend(), comparePhotonsByPt);

  }//End of photon loop
 
//***CHECK sorting from the triphotons vector 
//triphotons vector contents: for(auto i = triphotons.begin(); i != triphotons.end(); i++) cout<< i->pt() <<endl;
for (int n = 0; n < 3; n++){
 cout<< photons[n]->pt()<<endl;
}

//We only fill tree for events with at least three photons: 

   if (photons->size()>2) fTree->Fill();
cout << "======================================RUN ENDS==================================" <<endl;


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
threephoton::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
threephoton::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
threephoton::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(threephoton);

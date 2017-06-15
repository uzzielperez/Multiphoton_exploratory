// -*- C++ -*-
//
// Package:    triphoton/triphoton
// Class:      triphoton
// 
/**\class triphoton triphoton.cc triphoton/triphoton/plugins/triphoton.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Cilicia Uzziel Perez
//         Created:  Sun, 11 Jun 2017 15:21:24 GMT
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

class triphoton : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit triphoton(const edm::ParameterSet&);
      ~triphoton();

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
     
     //PHOTON STRUCT 
     
     struct photonInfo_t {
       float pt;
       float eta;
       float phi;
       float sceta;
      float scphi; 
     };
     
     //Declare Instance of Photon1, Photon2, Photon3
     photonInfo_t fPhoton1Info;
     photonInfo_t fPhoton2Info;
     photonInfo_t fPhoton3Info;
     
     photonInfo_t photon1;
     photonInfo_t photon2;
     photonInfo_t photon3;

     //Using a Container for the first time
     std::vector<photonInfo_t> photon1info;
     std::vector<photonInfo_t> photon2info;
     std::vector<photonInfo_t> photon3info;

     // ExoDiPhotons::photonInfo_t fPhoton1Info; // leading
     // ExoDiPhotons::photonInfo_t fPhoton2Info; // subleading
     //  ExoDiPhotons::photonInfo_t fPhoton3Info; // leading
        
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
triphoton::triphoton(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   usesResource("TFileService");
  
   edm::Service<TFileService> fs; 
   fTree = fs->make<TTree>("fTree","TriphotonTree");
   fTree->Branch("Event",&fEventInfo, "run/L:LS:evnum");
   fTree->Branch("Photon1", &fPhoton1Info, "pt/F:eta:phi:sceta:scphi");
   fTree->Branch("Photon2", &fPhoton2Info, "pt/F:eta:phi:sceta:scphi");
   fTree->Branch("Photon3", &fPhoton3Info, "pt/F:eta:phi:sceta:scphi");

  // fTree->Branch("Photon1",&fPhoton1Info,ExoDiPhotons::photonBranchDefString.c_str());
  // //fTree->Branch("Photon2",&fPhoton2Info,ExoDiPhotons::photonBranchDefString.c_str()); 
  // //fTree->Branch("Photon3",&fPhoton3Info,ExoDiPhotons::photonBranchDefString.c_str());  
   
   // photonsMiniAODToken_ = mayConsume<edm::View<pat::Photon>>(iConfig.getParameter<edm::InputTag>("photonsMiniAOD"));
  photonsMiniAODToken_ = mayConsume<edm::View<pat::Photon>>(edm::InputTag("slimmedPhotons"));
}


triphoton::~triphoton()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
triphoton::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   fEventInfo.run = -999;
   fEventInfo.LS = -999;
   fEventInfo.evnum = -999;
   
   ofstream cout("photon.txt", ios::app);
//   cout << endl;
   cout << "Run: " << iEvent.id().run() << ", LS: " << iEvent.id().luminosityBlock() << ", Event: " << iEvent.id().event() << endl;

   fEventInfo.run = iEvent.id().run();
   fEventInfo.LS = iEvent.id().luminosityBlock();
   fEventInfo.evnum = iEvent.id().event();

edm::Handle<edm::View<pat::Photon>> photons;   //To get photon collection
iEvent.getByToken(photonsMiniAODToken_,photons);
  //Initializing photon variables 

  fPhoton1Info.pt = -999;
  fPhoton1Info.eta = -999;
  fPhoton1Info.phi = -999;
  fPhoton2Info.pt = -999;
  fPhoton2Info.eta = -999;
  fPhoton2Info.phi = -999;
  fPhoton3Info.pt = -999;
  fPhoton3Info.eta = -999;
  fPhoton3Info.phi = -999;
  float lead_photon_pt = 0; 
  float sub_leading_photon = 0; 
  float sub_sub_photon = 0;
  float subthree = 0;
  fPhoton1Info.sceta = -999;
  fPhoton1Info.scphi = -999;
  fPhoton2Info.sceta = -999;
  fPhoton2Info.scphi = -999;
  fPhoton3Info.sceta = -999;
  fPhoton3Info.scphi = -999;
 
  //for (edm::View<pat::Photon>::const_iterator pho = photons->begin(); pho != photons->end(); ++pho) {
  //print out attributes 
  
  //const pat::Photon *photon = NULL;
  //sort(photons->pt());
 
  
 // ofstream cout("allPhotons.txt");
  //ofstream myfile;
  //myfile.open("photons.txt");

 // myfile.open("triphotons.txt");
 
  //const pat::Photon *photon = NULL; 
 // sort(photon.begin(), photon.end());
 
 //Create container for photons in each event 
  std::vector<int> photoncontainer; 
  for (size_t i = 0; i<photons->size(); ++i) {
    const auto pho = photons->ptrAt(i);
    cout << "Photon: " << "pt = " << pho->pt() << "; eta = " << pho->eta() <<"; phi = " <<pho->phi() << "; sceta = " << pho->superCluster()->eta() << "; scphi = " << pho->superCluster()->phi() << endl;
    //sort(photon->begin(), photon->end()); 
    //}
  //Sort and store info
 // for (size_t i = 0; i < photons->size(); ++i) { //loop over photon collection
   //   const auto pho = photons->ptrAt(i);
         
         //fill photon container 
         photoncontainer.push_back(pho->pt());
      

         //sort(pho->pt().begin(), pho->pt().end);
         //sort(photon, photon +photons->size());
    //   if(photons->size()>2){  
         
        //LOGICAL SORTING exercise: (Make sure they're sorted). 
        if (pho->pt() >= lead_photon_pt){
        sub_sub_photon = sub_leading_photon; sub_leading_photon = lead_photon_pt; lead_photon_pt = pho->pt(); 
        } 

        else if (pho->pt() > sub_leading_photon){
        sub_sub_photon = sub_leading_photon; sub_leading_photon = pho->pt();
        }

        else if (pho->pt() > sub_sub_photon){
        sub_sub_photon = pho->pt();
        }
        
         /*  if(pho->pt()>=lead_photon_pt){
           lead_photon_pt = pho->pt(); //sub_leading_photon = lead_photon_pt; sub_sub_photon = sub_leading_photon;
            } //photon with maxpt
         
         
         else if(pho->pt()<lead_photon_pt){
           if(pho->pt()>= sub_leading_photon){
             sub_leading_photon = pho->pt(); //sub_sub_photon = sub_leading_photon;
           }
         }
         else if(pho->pt()<sub_leading_photon){
            if(pho->pt() > sub_sub_photon){
             sub_sub_photon = pho->pt();
             }
           // else if (pho->pt()>subthree) continue;
            }*/
         
          //Photon Info
          //store Photon Info here  
          if (i== 0){
          fPhoton1Info.pt = pho->pt();
          fPhoton1Info.sceta = pho->superCluster()->eta();
          fPhoton1Info.scphi = pho->superCluster()->phi();
          fPhoton1Info.eta = pho->eta();
          fPhoton1Info.phi = pho->phi();
          photon1info.push_back(fPhoton1Info);       
          }
         if (i ==1){
          fPhoton2Info.pt = pho->pt();
          fPhoton2Info.sceta = pho->superCluster()->eta();
          fPhoton2Info.scphi = pho->superCluster()->phi();
          fPhoton2Info.eta = pho->eta();
          fPhoton2Info.phi = pho->phi();
          photon2info.push_back(fPhoton2Info);
         } 
         if (i==2){ 
           fPhoton3Info.pt = pho->pt();
           fPhoton3Info.sceta = pho->superCluster()->eta();
           fPhoton3Info.scphi = pho->superCluster()->phi();
           fPhoton3Info.eta = pho->eta();
           fPhoton3Info.phi = pho->phi();
          // photon3info.push_back(pho->pt);
         }
        // cout << "Photon: " << "pt = " << pho->pt() << "; eta = " << pho->eta() <<"; phi = " <<pho->phi() << "; sceta = " << pho->superCluster()->eta() << "; scphi = " << pho->superCluster()->phi() << endl;
        // }
        //
       // if(photons->size()>2){ 
      //  ofstream cout("minthreepho.txt", ios::app);
        //cout <<"TRI_Photon: " << "pt = " << pho->pt() << "; eta = " << pho->eta() <<"; phi = " <<pho->phi() << "; sceta = " << pho->superCluster()->eta() << "; scphi = " << pho->superCluster()->phi() << endl; 
       // } 
        //else continue;
  }

  //Print out sorted photon pt 
  cout << "lead_pt: " << lead_photon_pt << "; sub_leading_photon: " << sub_leading_photon << "; sub_sub_photon: " << sub_sub_photon << endl;  

  // Photon container now contains all the photon pts in one event 
 // sort(photoncontainer.begin(), photoncontainer.end()); 
 //sorted photons now look like 
 //for (std::vector <int>::iterator it = photoncontainer.begin() ; it != photoncontainer.end(); ++it){
 // cout << *it<<" " << photoncontainer.at(*it) <<endl; 
// }

  //cout << "lead_pt" << "pt = " << lead_photon_pt;
  //cout << "sub_lead" << "pt = " << sub_leading_photon;
  //cout << "subsub_lead" <<"pt = " << sub_sub_photon;
  // myfile.close();
  
  //std::cout << "My photon container contains:";
  //for (std::vector <int>::iterator it = photoncontainer.begin() ; it != photoncontainer.end(); ++it){
    //    std::cout << "photonpt: " << photoncontainer.at(it)<<endl;}
  //cout << "size" << ": " << photon3info.size();
 
  //print out Triphotoi containers
  //for(int j = 0 ; j< photon1info.size(); j++){cout<< photon1info.at(j) <<endl;}
  //for(int j = 0 ; j< photon2info.size(); j++){cout<< photon2info.at(j) <<endl;}
  //for(int j = 0 ; j< (photon3info.size(); j++){cout<< photon3info[j] <<endl;}
  //for(std::vector <int>::iterator it = photon1info.begin(); it !=photon1info.end(); it++)  std::cout << *it << ' '; 
   
  for (size_t i = 0; i<photons->size(); ++i) {
   const auto pho = photons->ptrAt(i);
    if(photons->size()>2){
    cout << "Photon: " << "pt = " << pho->pt() << "; eta = " << pho->eta() <<"; phi = " <<pho->phi() << "; sceta = " << pho->superCluster()->eta() << "; scphi = " << pho->superCluster()->phi() << endl;
    }
   }
  //Fill Tree for Events with at least 3 photons
  if (photons->size()>2){
  fTree->Fill();
}


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
triphoton::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
triphoton::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
triphoton::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(triphoton);

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "GlastEvent/TopLevel/Event.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "GlastEvent/TopLevel/EventModel.h"
#include "GlastEvent/MonteCarlo/McParticle.h"
#include "GlastEvent/MonteCarlo/McIntegratingHit.h"
#include "GlastEvent/MonteCarlo/McPositionHit.h"

#include "facilities/Util.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"

#include "mcRootData/McEvent.h"

#include <map>

/** @class mcRootWriterAlg
 * @brief Writes Monte Carlo TDS data to a persistent ROOT file.
 *
 * @author Heather Kelly
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/mcRootWriterAlg.cxx,v 1.3 2002/05/01 23:34:40 heather Exp $
 */

class mcRootWriterAlg : public Algorithm
{	
public:
    
    mcRootWriterAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    /// Handles setup by opening ROOT file in write mode and creating a new TTree
    StatusCode initialize();
   
    /// Orchastrates reading from TDS and writing to ROOT for each event
    StatusCode execute();
    
    /// Closes the ROOT file and cleans up
    StatusCode finalize();

private:

    /// Retrieves event Id and run Id from TDS and fills the McEvent ROOT object
    StatusCode writeMcEvent();

    /// Retrieves McParticles from TDS and write McParticles to the ROOT file
    StatusCode writeMcParticles();

    /// Retrieves McPositionHits from the TDS and writes McPostionHits 
    /// to the ROOT file
    StatusCode writeMcPositionHits();

    /// Retrieves McIntegratingHits from the TDS and 
    /// writes McIntegratingHits to the ROOT file
    StatusCode writeMcIntegratingHits();

    /// Converts idents::VolumeIdentifier into ROOT's VolumeIdentifier
    void convertVolumeId(idents::VolumeIdentifier tdsVolId, 
        VolumeIdentifier &rootVolId);

    /// Calls TTree::Fill for each event and clears m_mcEvt
    void writeEvent();

    /// Performs the final write to the ROOT file and closes
    void close();
   
    /// ROOT file pointer
    TFile *m_mcFile;
    /// ROOT tree pointer
    TTree *m_mcTree;
    /// Top-level Monte Carlo ROOT object
    McEvent *m_mcEvt;
    /// name of the output ROOT file
    std::string m_fileName;
    /// name of the TTree in the ROOT file
    std::string m_treeName;
    /// ROOT split mode
    int m_splitMode;
    /// Buffer Size for the ROOT file
    int m_bufSize;
    /// Compression level for the ROOT file
    int m_compressionLevel;

    /// Keep track of MC particles as we retrieve them from TDS
    /// Used only to aid in writing the ROOT file.
    std::map<const mc::McParticle*, McParticle*> m_particleMap;
    
};


static const AlgFactory<mcRootWriterAlg>  Factory;
const IAlgFactory& mcRootWriterAlgFactory = Factory;

mcRootWriterAlg::mcRootWriterAlg(const std::string& name, 
                                 ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator)
{
    // Input parameters available to be set via the jobOptions file
    declareProperty("mcRootFile",m_fileName="mc.root");
    declareProperty("splitMode", m_splitMode=1);
    declareProperty("bufferSize", m_bufSize=64000);
    // ROOT default compression
    declareProperty("compressionLevel", m_compressionLevel=1);
    declareProperty("treeName", m_treeName="Mc");

    m_particleMap.clear();
}

StatusCode mcRootWriterAlg::initialize()
{
    // Purpose and Method:  Called once before the run begins.  This method
    //    opens a new ROOT file and prepares for writing.

    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    // Use the Job options service to set the Algorithm's parameters
    // This will retrieve parameters set in the job options file
    setProperties();
    
    facilities::Util::expandEnvVar(&m_fileName);

    // Save the current directory for the ntuple writer service
    TDirectory *saveDir = gDirectory;   
    // Create the new ROOT file
    m_mcFile = new TFile(m_fileName.c_str(), "RECREATE");
    if (!m_mcFile->IsOpen()) sc = StatusCode::FAILURE;
    m_mcFile->cd();
    m_mcFile->SetCompressionLevel(m_compressionLevel);
    m_mcTree = new TTree(m_treeName.c_str(), "GLAST Monte Carlo Data");
    m_mcEvt = new McEvent();
    m_mcTree->Branch("McEvent","McEvent", &m_mcEvt, m_bufSize, m_splitMode);
    
    saveDir->cd();
    return sc;
    
}

StatusCode mcRootWriterAlg::execute()
{
    // Purpose and Method:  Called once per event.  This method calls
    //   the appropriate methods to read data from the TDS and write data
    //   to the ROOT file.

    MsgStream log(msgSvc(), name());

    StatusCode sc = StatusCode::SUCCESS;

    m_particleMap.clear();

    sc = writeMcEvent();
    if (sc.isFailure()) return sc;

    sc = writeMcParticles();
    if (sc.isFailure()) return sc;

    sc = writeMcPositionHits();
    if (sc.isFailure()) return sc;

    sc = writeMcIntegratingHits();
    if (sc.isFailure()) return sc;
   
    writeEvent();
    return sc;
}


StatusCode mcRootWriterAlg::writeMcEvent() {
    // Purpose and Method:  Retrieve the Event object from the TDS and set the
    //    event and run numbers in the McEvent ROOT object

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Retrieve the Event data for this event
    SmartDataPtr<Event> evt(eventSvc(), EventModel::Event);

    if (!evt) return sc;

    UInt_t evtId = evt->event();
    UInt_t runId = evt->run();

    log << MSG::DEBUG << evt->fillStream(log.stream()) << endreq;

    m_mcEvt->initialize(evtId, runId);

    return sc;
}

StatusCode mcRootWriterAlg::writeMcParticles() {
    // Purpose and Method:  Retrieve the McParticle collection from the TDS 
    //    and fill the ROOT McParticle collection
    // TDS Input:  EventModel::MC::McParticleCol
    // Output:  ROOT McParticle Collection

    MsgStream log(msgSvc(), name());

    StatusCode sc = StatusCode::SUCCESS;

    SmartDataPtr<mc::McParticleCol> particles(eventSvc(), EventModel::MC::McParticleCol);

    if (!particles) return sc;

    mc::McParticleCol::const_iterator p;

    // Create map of TDS McParticles and ROOT McParticles
    for (p = particles->begin(); p != particles->end(); p++) {
        log << MSG::DEBUG << (*p)->fillStream(log.stream()) << endreq;
        McParticle *mcPart = new McParticle();
        m_particleMap[(*p)] = mcPart;
    }

    // Now that the full map of McParticles is created 
    // we initialize and save the ROOT McParticles 
    for (p = particles->begin(); p != particles->end(); p++) {
        
        McParticle *mcPart = m_particleMap[(*p)];
        // particleProperty() returns a StdHepId which is defined as an int
        Int_t idRoot = (*p)->particleProperty();
        
        UInt_t statFlagsRoot = (*p)->statusFlags();
                
        HepPoint3D finalPosTds = (*p)->finalPosition();
        TVector3 finalPosRoot(finalPosTds.x(), finalPosTds.y(), finalPosTds.z());
        
        HepLorentzVector initMomTds = (*p)->initialFourMomentum();
        TLorentzVector initMomRoot(initMomTds.x(), initMomTds.y(), initMomTds.z(), initMomTds.t());
        
        HepLorentzVector finalMomTds = (*p)->finalFourMomentum();
        TLorentzVector finalMomRoot(finalMomTds.x(), finalMomTds.y(), finalMomTds.z(), finalMomTds.t());

        const mc::McParticle *momTds = &((*p)->mother());
        McParticle *momRoot = 0;
        if (momTds != 0) {
            // The case where this is the primary particle
            if (momTds == (*p)) {
                momRoot = mcPart;
            } else if (m_particleMap.find(momTds) != m_particleMap.end()) {
                // Otherwise we retrieve the McParticle from the map
                momRoot = m_particleMap[momTds];
            } else {
                log << MSG::WARNING << "Did not find mother McParticle in the map!!" << endreq;
            }
        } else {
            log << MSG::WARNING << "TDS McParticle Mother is null" << endreq;
        }
                
        // Setup the ROOT McParticle
        mcPart->initialize(momRoot, idRoot, statFlagsRoot, initMomRoot, finalMomRoot, finalPosRoot);
        // Add the ROOT McParticle to the ROOT collection of McParticle
        m_mcEvt->addMcParticle(mcPart);     

    }

    return sc;
}


StatusCode mcRootWriterAlg::writeMcPositionHits() {
    // Purpose and Method:  Retrieve the McPositionHit collection from the TDS 
    //  and fill the ROOT McPositionHit collection.
    // TDS Input:  EventModel::MC::McPositionHitCol
    // Output:  ROOT McPositionHit collection

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Get the McPositionHits Collection from the TDS
    SmartDataPtr<McPositionHitVector> posHits(eventSvc(), EventModel::MC::McPositionHitCol);
    if (!posHits) return sc;
    
    log << MSG::DEBUG << "Number of McPositionHits in the event = " << posHits->size() << endreq;
    
    McPositionHitVector::const_iterator hit;
    for (hit = posHits->begin(); hit != posHits->end(); hit++ ) {
    
        log << MSG::DEBUG << (*hit)->fillStream(log.stream()) << endreq;

        idents::VolumeIdentifier volIdTds = (*hit)->volumeID();
        VolumeIdentifier volIdRoot;
        convertVolumeId(volIdTds, volIdRoot);

        HepPoint3D entryTds = (*hit)->entryPoint();
        TVector3 entryRoot(entryTds.x(), entryTds.y(), entryTds.z());

        HepPoint3D exitTds = (*hit)->exitPoint();
        TVector3 exitRoot(exitTds.x(), exitTds.y(), exitTds.z());
        
        Double_t edepRoot = (*hit)->depositedEnergy();
        
        Double_t epartRoot = (*hit)->particleEnergy();
        
        //bool primaryOrigin = (*hit)->primaryOrigin();
        //log << MSG::INFO << "primaryOrigin " << primaryOrigin << endreq;
        
        //bool caloShowerOrigin = (*hit)->caloShowerOrigin();
        //log << MSG::INFO << "calShowerO " << caloShowerOrigin << endreq;
        
        //bool needDigi = (*hit)->needDigi();
        //log << MSG::INFO << "need Digi " << needDigi << endreq;
                
        Double_t tofRoot = (*hit)->timeOfFlight();
        
        const mc::McParticle *mcTds = (*hit)->mcParticle();
        McParticle *mcRoot = 0;
        if (mcTds != 0) {
            mcRoot = m_particleMap[mcTds];
        }
      
        const mc::McParticle *originTds = (*hit)->originMcParticle();
        McParticle *originRoot = 0;
        if (originTds != 0) {
            originRoot = m_particleMap[originTds];
        }

        UInt_t flagsRoot = 0;

        McPositionHit *mcPosHit = new McPositionHit();
        // Setup the ROOT McPositionHit
        mcPosHit->initialize(edepRoot, volIdRoot, entryRoot, exitRoot, mcRoot, originRoot, epartRoot, tofRoot, flagsRoot);
        // Add the ROOT McPositionHit to the ROOT collection of McPositionHits
        m_mcEvt->addMcPositionHit(mcPosHit);
    }
    
    return sc;
}

StatusCode mcRootWriterAlg::writeMcIntegratingHits() {
    // Purpose and Method:  Retrieve the McIntegratingHit collection from the
    //     TDS and fill the ROOT McIntegratingHit collection.
    // TDS Input:  EventModel::MC::McIntegratingHitCol
    // Output:  ROOT McIntegratingHit collection

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Retrieve the McIntegratingHit collection from the TDS
    SmartDataPtr<McIntegratingHitVector> intHits(eventSvc(), 
        EventModel::MC::McIntegratingHitCol);
    if (!intHits) return sc;
    
    
    log << MSG::DEBUG << "Number of McIntegratingHits in the event = " 
        << intHits->size() << endreq;
    
    McIntegratingHitVector::const_iterator hit;
    for (hit = intHits->begin(); hit != intHits->end(); hit++ ) {

        log << MSG::DEBUG << (*hit)->fillStream(log.stream()) << endreq;

        const idents::VolumeIdentifier idTds = (*hit)->volumeID();
        VolumeIdentifier idRoot;
        convertVolumeId(idTds, idRoot);

        Double_t e = (*hit)->totalEnergy();
        
        const HepPoint3D moment1Tds = (*hit)->moment1();
        TVector3 moment1Root(moment1Tds.x(), moment1Tds.y(), moment1Tds.z());

        const HepPoint3D moment2Tds = (*hit)->moment2();
        TVector3 moment2Root(moment2Tds.x(), moment2Tds.y(), moment2Tds.z());

        McIntegratingHit *mcIntHit = new McIntegratingHit();

        // Setup the ROOT McIntegratingHit
        mcIntHit->initialize(idRoot);

        mc::McIntegratingHit::energyDepositMap tdsMap = (*hit)->itemizedEnergy();
        mc::McIntegratingHit::energyDepositMap::const_iterator mapIt;
        for (mapIt = tdsMap.begin(); mapIt != tdsMap.end(); mapIt++) {
            mc::McParticle *mcPartTds = mapIt->first;
            McParticle *mcPartRoot = m_particleMap[mcPartTds];
            Double_t e = mapIt->second;
            HepPoint3D posTds = mcPartTds->finalPosition();
            TVector3 posRoot(posTds.x(), posTds.y(), posTds.z());
            mcIntHit->addEnergyItem(e, mcPartRoot, posRoot);
        }

        // Add the ROOT McIntegratingHit to the ROOT collection of McIntegratingHits
        m_mcEvt->addMcIntegratingHit(mcIntHit);
    }

    return sc;
}

void mcRootWriterAlg::convertVolumeId(idents::VolumeIdentifier tdsVolId, 
                     VolumeIdentifier& rootVolId) 
{
    // Purpose and Method:  We must store the volume ids as two 32 bit UInt_t
    //     in the ROOT class.  Hence, we must convert the 64 bit representation
    //     used in the idents::VolumeIdentifier class into two 32 bit UInt_t.
    //     To perform the conversion, we iterate over all the ids in the TDS
    //     version of the idents::VolumeIdentifier and append each to the ROOT
    //     VolumeIdentifier
    
    unsigned int index;
    rootVolId.Clear();
    for (index = 0; index < tdsVolId.size(); index++) {
        rootVolId.append(tdsVolId.operator [](index));
    }
}


void mcRootWriterAlg::writeEvent() 
{
    // Purpose and Method:  Stores the McEvent data for this event in the ROOT
    //    tree.  The m_mcEvt object is cleared for the next event.

    TDirectory *saveDir = gDirectory;
    m_mcFile->cd();
    m_mcTree->Fill();
    m_mcEvt->Clear();
    saveDir->cd();
    return;
}

void mcRootWriterAlg::close() 
{
    // Purpose and Method:  Writes the ROOT file at the end of the run.
    //    The TObject::kOverWrite parameter is used in the Write method
    //    since ROOT will periodically write to the ROOT file when the bufSize
    //    is filled.  Writing would create 2 copies of the same tree to be
    //    stored in the ROOT file, if we did not specify kOverwrite.

    TDirectory *saveDir = gDirectory;
    m_mcFile->cd();
    m_mcFile->Write(0, TObject::kOverwrite);
    m_mcFile->Close();
    saveDir->cd();
    return;
}

StatusCode mcRootWriterAlg::finalize()
{
    close();
    
    StatusCode sc = StatusCode::SUCCESS;
    return sc;
}


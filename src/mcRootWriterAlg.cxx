#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/MCEvent.h"
//#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/MonteCarlo/McPositionHit.h"

#include "facilities/Util.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TProcessID.h"

#include "mcRootData/McEvent.h"

#include "commonData.h"
#include "RootIo/IRootIoSvc.h"

#include <map>

// ADDED FOR THE FILE HEADERS DEMO
#include "RootIo/FhTool.h"
#include <cstdlib>

/** @class mcRootWriterAlg
 * @brief Writes Monte Carlo TDS data to a persistent ROOT file.
 *
 * @author Heather Kelly
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/mcRootWriterAlg.cxx,v 1.32 2004/09/24 19:14:05 heather Exp $
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
    //std::map<const Event::McParticle*, McParticle*> m_mcPartMap;

    commonData m_common;
    IRootIoSvc* m_rootIoSvc;

    // ADDED FOR THE FILE HEADERS DEMO
    IFhTool * m_headersTool ;
};


static const AlgFactory<mcRootWriterAlg>  Factory;
const IAlgFactory& mcRootWriterAlgFactory = Factory;

mcRootWriterAlg::mcRootWriterAlg(const std::string& name, 
                                 ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator)
{
    // Input parameters available to be set via the jobOptions file
    // Output ROOT file name
    declareProperty("mcRootFile",m_fileName="mc.root");
    // Split mode for writing the TTree [0,99]
    declareProperty("splitMode", m_splitMode=1);
    // Buffer size for writing ROOT data
    declareProperty("bufferSize", m_bufSize=64000);
    // ROOT default compression
    declareProperty("compressionLevel", m_compressionLevel=1);
    // ROOT TTree name
    declareProperty("treeName", m_treeName="Mc");

    
    // ADDED FOR THE FILE HEADERS DEMO
    m_headersTool = 0 ;

    m_common.m_mcPartMap.clear();
    m_common.m_mcPosHitMap.clear();
    m_common.m_mcIntHitMap.clear();
}

StatusCode mcRootWriterAlg::initialize()
{
    // Purpose and Method:  Called once before the run begins.  This method
    //    opens a new ROOT file and prepares for writing.

    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    // ADDED FOR THE FILE HEADERS DEMO
    StatusCode headersSc = toolSvc()->retrieveTool("FhTool",m_headersTool) ;
    if (headersSc.isFailure()) {
        log<<MSG::WARNING << "Failed to retreive headers tool" << endreq;
    }
    headersSc = m_headersTool->newMcHeader() ;
    if (headersSc.isFailure()) {
        log<<MSG::WARNING << "Failed to create a new Mc FileHeader" << endreq;
    }
    
    // Use the Job options service to set the Algorithm's parameters
    // This will retrieve parameters set in the job options file
    setProperties();

    if ( service("RootIoSvc", m_rootIoSvc, true).isFailure() ){
        log << MSG::INFO << "Couldn't find the RootIoSvc!" << endreq;
        log << MSG::INFO << "No Auto Saving" << endreq;
        m_rootIoSvc = 0;
    } 

    facilities::Util::expandEnvVar(&m_fileName);

    // Save the current directory for the ntuple writer service
    TDirectory *saveDir = gDirectory;   
    // Create the new ROOT file
    m_mcFile = new TFile(m_fileName.c_str(), "RECREATE");
    if (!m_mcFile->IsOpen()) {
        log << MSG::ERROR << "ROOT file " << m_fileName 
            << " could not be opened for writing." << endreq;
        return StatusCode::FAILURE;
    }

    m_mcFile->cd();
    m_mcFile->SetCompressionLevel(m_compressionLevel);
    m_mcTree = new TTree(m_treeName.c_str(), "GLAST Monte Carlo Data");
    m_mcEvt = new McEvent();
    m_mcTree->Branch("McEvent","McEvent", &m_mcEvt, m_bufSize, m_splitMode);

    m_common.m_mcEvt = m_mcEvt;

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
    
    if (!m_mcFile->IsOpen()) {
        log << MSG::ERROR << "ROOT file " << m_fileName 
            << " could not be opened for writing." << endreq;
        return StatusCode::FAILURE;
    }

    Int_t ObjectNumber = TProcessID::GetObjectCount();

    m_common.m_mcPartMap.clear();
    m_common.m_mcPosHitMap.clear();
    m_common.m_mcIntHitMap.clear();

    sc = writeMcEvent();
    if (sc.isFailure()) return sc;

    sc = writeMcParticles();
    if (sc.isFailure()) return sc;

    sc = writeMcPositionHits();
    if (sc.isFailure()) return sc;

    sc = writeMcIntegratingHits();
    if (sc.isFailure()) return sc;
   
    writeEvent();

    // reset object nr in order to avoid memleak
    TProcessID::SetObjectCount(ObjectNumber);

    return sc;
}


StatusCode mcRootWriterAlg::writeMcEvent() {
    // Purpose and Method:  Retrieve the Event object from the TDS and set the
    //    event and run numbers in the McEvent ROOT object
    
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    // Retrieve the Event data for this event
    SmartDataPtr<Event::EventHeader> evt(eventSvc(), EventModel::EventHeader);
    
    if (!evt) return sc;
    m_mcEvt->Clear();
    
    UInt_t evtId = evt->event();
    UInt_t runId = evt->run();
    
    SmartDataPtr<Event::MCEvent> mcEvt(eventSvc(), EventModel::MC::Event);
    Int_t sourceId = mcEvt->getSourceId();
    UInt_t sequence = mcEvt->getSequence();
    Double_t timestamp = mcEvt->time();  /// note: converting from timestamp to double
    
    log << MSG::DEBUG;
    if( log.isActive()) evt->fillStream(log.stream());
    log << endreq;
    
    m_mcEvt->initialize(evtId, runId, sourceId, sequence, timestamp);
    
    return sc;
}

StatusCode mcRootWriterAlg::writeMcParticles() {
    // Purpose and Method:  Retrieve the McParticle collection from the TDS 
    //    and fill the ROOT McParticle collection
    // TDS Input:  EventModel::EVENT::McParticleCol
    // Output:  ROOT McParticle Collection

    MsgStream log(msgSvc(), name());

    StatusCode sc = StatusCode::SUCCESS;

    SmartDataPtr<Event::McParticleCol> particles(eventSvc(), 
        EventModel::MC::McParticleCol);

    if (!particles) return sc;

    Event::McParticleCol::const_iterator p;

    // Create map of TDS McParticles and ROOT McParticles
    for (p = particles->begin(); p != particles->end(); p++) {
        log << MSG::DEBUG;
        if( log.isActive())(*p)->fillStream(log.stream());
        log << endreq;
        McParticle *mcPart = new McParticle();
        TRef ref = mcPart;
        m_common.m_mcPartMap[(*p)] = ref;
    }

    // Now that the full map of McParticles is created 
    // we initialize and save the ROOT McParticles 
    for (p = particles->begin(); p != particles->end(); p++) {
        
        TRef ref = m_common.m_mcPartMap[(*p)];
        McParticle *mcPart = (McParticle*)ref.GetObject();
        // particleProperty() returns a StdHepId which is defined as an int
        Int_t idRoot = (*p)->particleProperty();
        
        UInt_t statFlagsRoot = (*p)->statusFlags();
                
        HepPoint3D initPosTds = (*p)->initialPosition();
        TVector3 initPosRoot(initPosTds.x(), initPosTds.y(), initPosTds.z());

        HepPoint3D finalPosTds = (*p)->finalPosition();
        TVector3 finalPosRoot(finalPosTds.x(), finalPosTds.y(), finalPosTds.z());
        
        HepLorentzVector initMomTds = (*p)->initialFourMomentum();
        TLorentzVector initMomRoot(initMomTds.x(), initMomTds.y(), 
            initMomTds.z(), initMomTds.t());
        
        HepLorentzVector finalMomTds = (*p)->finalFourMomentum();
        TLorentzVector finalMomRoot(finalMomTds.x(), finalMomTds.y(), 
            finalMomTds.z(), finalMomTds.t());

        const Event::McParticle *momTds = &((*p)->mother());
        McParticle *momRoot = 0;
        if (momTds != 0) {
            // The case where this is the primary particle
            if (momTds == (*p)) {
                momRoot = mcPart;
            } else if (m_common.m_mcPartMap.find(momTds) != m_common.m_mcPartMap.end()) {
                // Otherwise we retrieve the McParticle from the map
                TRef ref = m_common.m_mcPartMap[momTds];
                momRoot = (McParticle*)ref.GetObject();//m_common.m_mcPartMap[momTds];
            } else {
                log << MSG::WARNING << "Did not find mother McParticle in the"
                    << "map!!" << endreq;
            }
        } else {
            log << MSG::WARNING << "TDS McParticle Mother is null" << endreq;
        }
                
        // Setup the ROOT McParticle
        mcPart->initialize(momRoot, idRoot, statFlagsRoot, initMomRoot, 
            finalMomRoot, initPosRoot, finalPosRoot, (*p)->getProcess().c_str());
        
        /* HMK This is probably no longer necessary - as this is handled in McParticle::initialize
        // Process the Daughter Particles
        const SmartRefVector<Event::McParticle> daughterCol = (*p)->daughterList();
        SmartRefVector<Event::McParticle>::const_iterator daughterIt;
        for (daughterIt = daughterCol.begin(); daughterIt != daughterCol.end(); daughterIt++) {
            if (m_mcPartMap.find((*daughterIt)) != m_mcPartMap.end()) {
                McParticle *daughter = m_mcPartMap[(*daughterIt)];
                mcPart->addDaughter(daughter);
            } else {
                log << MSG::WARNING << "Did not find daughter McParticle in the"
                    << " map!" << endreq;
            }
        }
        */

        // Add the ROOT McParticle to the ROOT collection of McParticle
        m_mcEvt->addMcParticle(mcPart);     

    }

    return sc;
}


StatusCode mcRootWriterAlg::writeMcPositionHits() {
    // Purpose and Method:  Retrieve the McPositionHit collection from the TDS 
    //  and fill the ROOT McPositionHit collection.
    // TDS Input:  EventModel::EVENT::McPositionHitCol
    // Output:  ROOT McPositionHit collection

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Get the McPositionHits Collection from the TDS
    SmartDataPtr<Event::McPositionHitVector> posHits(eventSvc(), 
        EventModel::MC::McPositionHitCol);
    if (!posHits) return sc;
    
    log << MSG::DEBUG;
    if( log.isActive()) 
        log.stream() << "Number of McPositionHits in the event = " << posHits->size();
    log << endreq;
    
    Event::McPositionHitVector::const_iterator hit;
    for (hit = posHits->begin(); hit != posHits->end(); hit++ ) {
    
        log << MSG::DEBUG;
        if( log.isActive()) (*hit)->fillStream(log.stream());
        log << endreq;

        Int_t particleId = (*hit)->getMcParticleId();

        Int_t originPartId = (*hit)->getOriginMcParticleId();

        idents::VolumeIdentifier volIdTds = (*hit)->volumeID();
        VolumeIdentifier volIdRoot;
        convertVolumeId(volIdTds, volIdRoot);

        HepPoint3D entryTds = (*hit)->entryPoint();
        TVector3 entryRoot(entryTds.x(), entryTds.y(), entryTds.z());

        HepPoint3D exitTds = (*hit)->exitPoint();
        TVector3 exitRoot(exitTds.x(), exitTds.y(), exitTds.z());

        HepPoint3D globalEntryTds = (*hit)->globalEntryPoint();
        TVector3 globalEntryRoot(globalEntryTds.x(), globalEntryTds.y(), globalEntryTds.z());

        HepPoint3D globalExitTds = (*hit)->globalExitPoint();
        TVector3 globalExitRoot(globalExitTds.x(), globalExitTds.y(), globalExitTds.z());
        
        Double_t edepRoot = (*hit)->depositedEnergy();
        
        Double_t epartRoot = (*hit)->particleEnergy();
        
        //bool primaryOrigin = (*hit)->primaryOrigin();
        //log << MSG::INFO << "primaryOrigin " << primaryOrigin << endreq;
        
        //bool caloShowerOrigin = (*hit)->caloShowerOrigin();
        //log << MSG::INFO << "calShowerO " << caloShowerOrigin << endreq;
        
        //bool needDigi = (*hit)->needDigi();
        //log << MSG::INFO << "need Digi " << needDigi << endreq;
                
        Double_t tofRoot = (*hit)->timeOfFlight();
              
        const Event::McParticle *mcPartTds = (*hit)->mcParticle();
        McParticle *mcPartRoot = 0;
        if (mcPartTds != 0) {
            TRef ref;
            ref = m_common.m_mcPartMap[mcPartTds];
            mcPartRoot = (McParticle*)ref.GetObject();
        }

        const Event::McParticle *originTds = (*hit)->originMcParticle();
        McParticle *originRoot = 0;
        if (originTds != 0) {
            TRef ref;
            ref = m_common.m_mcPartMap[originTds];
            originRoot = (McParticle*)ref.GetObject();
        }

        UInt_t flagsRoot = 0;

        // Set up the ROOT McPositionHit
        McPositionHit *mcPosHit = new McPositionHit();

        TRef ref = mcPosHit;
        m_common.m_mcPosHitMap[(*hit)] = ref;

        mcPosHit->initialize(particleId, originPartId, edepRoot, volIdRoot, 
            entryRoot, exitRoot, globalEntryRoot, globalExitRoot, 
            mcPartRoot, originRoot, epartRoot, tofRoot, flagsRoot);

        // Add the ROOT McPositionHit to the ROOT collection of McPositionHits
        m_mcEvt->addMcPositionHit(mcPosHit);
    }
    
    return sc;
}

StatusCode mcRootWriterAlg::writeMcIntegratingHits() {
    // Purpose and Method:  Retrieve the McIntegratingHit collection from the
    //     TDS and fill the ROOT McIntegratingHit collection.
    // TDS Input:  EventModel::EVENT::McIntegratingHitCol
    // Output:  ROOT McIntegratingHit collection

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Retrieve the McIntegratingHit collection from the TDS
    SmartDataPtr<Event::McIntegratingHitVector> intHits(eventSvc(), 
        EventModel::MC::McIntegratingHitCol);
    if (!intHits) return sc;
    
    log << MSG::DEBUG;
    if(log.isActive()) 
        log.stream() << "Number of McIntegratingHits in the event = " << intHits->size();
    log << endreq;
    
    Event::McIntegratingHitVector::const_iterator hit;
    for (hit = intHits->begin(); hit != intHits->end(); hit++ ) {

        log << MSG::DEBUG;
        if( log.isActive()) (*hit)->fillStream(log.stream());
        log << endreq;

        const idents::VolumeIdentifier idTds = (*hit)->volumeID();
        VolumeIdentifier idRoot;
        convertVolumeId(idTds, idRoot);

        Double_t e = (*hit)->totalEnergy();
        
        const HepPoint3D moment1Tds = (*hit)->moment1();
        TVector3 moment1Root(moment1Tds.x(), moment1Tds.y(), moment1Tds.z());

        const HepPoint3D moment2Tds = (*hit)->moment2();
        TVector3 moment2Root(moment2Tds.x(), moment2Tds.y(), moment2Tds.z());

        McIntegratingHit *mcIntHit = new McIntegratingHit();

        TRef ref = mcIntHit;
        m_common.m_mcIntHitMap[(*hit)] = ref;

        // Setup the ROOT McIntegratingHit
        mcIntHit->initialize(idRoot);
        
        Event::McIntegratingHit::energyDepositMap tdsMap = (*hit)->itemizedEnergy();
        if (tdsMap.size() <= 0) {
            // If we have not stored the full tree, assume we are using the minimal tree
            mcIntHit->setEnergyItems(e, (*hit)->energyArray(), moment1Root, moment2Root);
            
        } else {
            Event::McIntegratingHit::energyDepositMap::const_iterator mapIt;
            for (mapIt = tdsMap.begin(); mapIt != tdsMap.end(); mapIt++) {
                Event::McParticle *mcPartTds = mapIt->first;
                TRef ref = m_common.m_mcPartMap[mcPartTds];
                McParticle *mcPartRoot = (McParticle*)ref.GetObject();
                Double_t e = mapIt->second;
                HepPoint3D posTds = mcPartTds->finalPosition();
                TVector3 posRoot(posTds.x(), posTds.y(), posTds.z());
                mcIntHit->addEnergyItem(e, mcPartRoot, posRoot);
            }
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
    
    int index;
    rootVolId.Clear();
    for (index = 0; index < tdsVolId.size(); index++) {
        rootVolId.append(tdsVolId.operator [](index));
    }
}


void mcRootWriterAlg::writeEvent() 
{
    // Purpose and Method:  Stores the McEvent data for this event in the ROOT
    //    tree.  The m_mcEvt object is cleared for the next event.

    static int eventCounter = 0;
    TDirectory *saveDir = gDirectory;
    m_mcTree->GetCurrentFile()->cd();
    //m_mcFile->cd();
    m_mcTree->Fill();
    //m_mcEvt->Clear();
    saveDir->cd();
    ++eventCounter;
    if (m_rootIoSvc)
        if (eventCounter % m_rootIoSvc->getAutoSaveInterval()== 0) m_mcTree->AutoSave();
    return;
}

void mcRootWriterAlg::close() 
{
    // Purpose and Method:  Writes the ROOT file at the end of the run.
    //    The TObject::kWriteDelete parameter is used in the Write method
    //    Used rather than TObject::kOverwrite - supposed to be safer but slower
    //    since ROOT will periodically write to the ROOT file when the bufSize
    //    is filled.  Writing would create 2 copies of the same tree to be
    //    stored in the ROOT file, if we did not specify kOverwrite.

    TDirectory *saveDir = gDirectory;
    TFile *f = m_mcTree->GetCurrentFile();
    f->cd();
    //m_mcFile->cd();
    m_mcTree->BuildIndex("m_runId", "m_eventId");
    f->Write(0, TObject::kWriteDelete);
    f->Close();
    saveDir->cd();
    return;
}

StatusCode mcRootWriterAlg::finalize()
{
    // ADDED FOR THE FILE HEADERS DEMO
    m_headersTool->writeMcHeader(m_mcTree->GetCurrentFile()) ;
    
    close();
    
    StatusCode sc = StatusCode::SUCCESS;
    setFinalized();
    return sc;
}

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/MCEvent.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/MonteCarlo/McPositionHit.h"
#include "Event/MonteCarlo/McTrajectory.h"
#include "Event/MonteCarlo/McTkrStrip.h"

#include "facilities/Util.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TProcessID.h"
#include "TString.h"

#include "mcRootData/McEvent.h"

#include "commonData.h"
#include "RootIo/IRootIoSvc.h"

// ADDED FOR THE FILE HEADERS DEMO
#include "RootIo/FhTool.h"
#include <cstdlib>

// low level converters
#include <RootConvert/MonteCarlo/McPositionHitConvert.h>
#include <RootConvert/MonteCarlo/McTrajectoryConvert.h>
#include <RootConvert/MonteCarlo/McTkrStripConvert.h>
#include <RootConvert/Utilities/Toolkit.h>

#include <map>

/** @class mcRootWriterAlg
 * @brief Writes Monte Carlo TDS data to a persistent ROOT file.
 *
 * @author Heather Kelly
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/RootIo/src/mcRootWriterAlg.cxx,v 1.56 2011/01/16 03:41:21 usher Exp $
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

    /// Retrieves McTrajectory's from the TDS and writes McTrajectory's 
    /// to the ROOT file
    StatusCode writeMcTrajectories();

    /// Retrieves McTkrStrip's from the TDS and writes McTkrStrip's
    /// to the ROOT file
    StatusCode writeMcTkrStrips();

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
    /// Flag to write McTrajectory objects
    bool m_saveTrajectories;
    /// Flag to write McTkrStrips
    bool m_saveMcTkrStrips;

    commonData m_common;
    IRootIoSvc* m_rootIoSvc;
    /// Option string which will be passed to McEvent::Clear
    std::string m_clearOption;

    // ADDED FOR THE FILE HEADERS DEMO
    IFhTool * m_headersTool ;
};


//static const AlgFactory<mcRootWriterAlg>  Factory;
//const IAlgFactory& mcRootWriterAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(mcRootWriterAlg);

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
    // Save McTrajectory's?
    declareProperty("saveMcTrajectory", m_saveTrajectories=true);
    declareProperty("saveMcTkrStrips", m_saveMcTkrStrips=false);
    declareProperty("clearOption", m_clearOption="");

    
    // ADDED FOR THE FILE HEADERS DEMO
    m_headersTool = 0 ;

    m_common.m_mcPartMap.clear();
    m_common.m_mcPosHitMap.clear();
    m_common.m_mcIntHitMap.clear();
    m_common.m_mcTrajectoryMap.clear();
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
    
    // Use the Job options service to set the Algorithm's parameters
    // This will retrieve parameters set in the job options file
    setProperties();

    m_rootIoSvc = 0 ;
    if ( service("RootIoSvc", m_rootIoSvc, true).isFailure() ){
        log << MSG::WARNING << "Couldn't find the RootIoSvc!" << endreq;
        m_rootIoSvc = 0;
        // Need RootIoSvc for writing/reading, cannot continue without it
        return StatusCode::FAILURE;
    } 

    m_mcTree = m_rootIoSvc->prepareRootOutput("mc", m_fileName, m_treeName, 
        m_compressionLevel, "GLAST Monte Carlo Data");

    m_mcEvt = new McEvent();
    m_rootIoSvc->setupBranch("mc", "McEvent", "McEvent", &m_mcEvt, m_bufSize, m_splitMode);

    m_common.m_mcEvt = m_mcEvt;

    return sc;
    
}

StatusCode mcRootWriterAlg::execute()
{
    // Purpose and Method:  Called once per event.  This method calls
    //   the appropriate methods to read data from the TDS and write data
    //   to the ROOT file.

    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG;
    if( log.isActive()) log<<"execute"<<endreq ;
    log << endreq;

    StatusCode sc = StatusCode::SUCCESS;
    
    if (!m_mcTree->GetCurrentFile()->IsOpen()) {
        log << MSG::ERROR << "ROOT file " << m_fileName 
            << " could not be opened for writing." << endreq;
        return StatusCode::FAILURE;
    }

    m_common.m_mcPartMap.clear();
    m_common.m_mcPosHitMap.clear();
    m_common.m_mcIntHitMap.clear();
    m_common.m_mcTrajectoryMap.clear();
    m_common.m_mcTrajectoryPointMap.clear();

    sc = writeMcEvent();
    if (sc.isFailure()) return sc;

    sc = writeMcParticles();
    if (sc.isFailure()) return sc;

    sc = writeMcPositionHits();
    if (sc.isFailure()) return sc;

    sc = writeMcIntegratingHits();
    if (sc.isFailure()) return sc;

    if (m_saveTrajectories)
    {
        sc = writeMcTrajectories();
        if (sc.isFailure()) return sc;
    }

    if (m_saveMcTkrStrips)
    {
        sc = writeMcTkrStrips();
        if (sc.isFailure()) return sc;
    }
   
    writeEvent();


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
    m_mcEvt->Clear(m_clearOption.c_str());
    
    UInt_t evtId = evt->event();
    UInt_t runId = evt->run();
    
    SmartDataPtr<Event::MCEvent> mcEvt(eventSvc(), EventModel::MC::Event);
    Int_t sourceId = mcEvt->getSourceId();
    UInt_t sequence = mcEvt->getSequence();
    Double_t timestamp = mcEvt->time();  /// note: converting from timestamp to double

    TString sourceName = mcEvt->getSourceName();
    
    log << MSG::DEBUG;
    if( log.isActive()) evt->fillStream(log.stream());
    log << endreq;
    
    m_mcEvt->initialize(evtId, runId, sourceId, sequence, timestamp, sourceName);
    
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
        
        CLHEP::HepLorentzVector initMomTds = (*p)->initialFourMomentum();
        TLorentzVector initMomRoot(initMomTds.x(), initMomTds.y(), 
            initMomTds.z(), initMomTds.t());
        
        CLHEP::HepLorentzVector finalMomTds = (*p)->finalFourMomentum();
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
              
        McPositionHit *mcPosHit = new McPositionHit();
        TRef ref = mcPosHit;
        m_common.m_mcPosHitMap[(*hit)] = ref ;
        RootPersistence::convert(**hit,*mcPosHit) ;  
         
        const Event::McParticle *mcPartTds = (*hit)->mcParticle();
        McParticle *mcPartRoot = 0;
        if (mcPartTds != 0) {
            TRef ref;
            ref = m_common.m_mcPartMap[mcPartTds];
            mcPartRoot = (McParticle*)ref.GetObject();
        }
        mcPosHit->setMcParticle(mcPartRoot) ;

        const Event::McParticle *originTds = (*hit)->originMcParticle();
        McParticle *originRoot = 0;
        if (originTds != 0) {
            TRef ref;
            ref = m_common.m_mcPartMap[originTds];
            originRoot = (McParticle*)ref.GetObject();
        }
        mcPosHit->setOriginMcParticle(originRoot) ;

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
        idRoot = RootPersistence::convert(idTds) ;

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

        // Now copy the XtalEnergyDep objects from the map
        for(Event::McIntegratingHit::XtalEnergyDepMap::const_iterator depItr  = (*hit)->getXtalEnergyDepMap().begin();
                                                                      depItr != (*hit)->getXtalEnergyDepMap().end();
                                                                      depItr++)
        {
            std::string evtKey = depItr->first;
            TObjString* key = new TObjString(evtKey.c_str());
            McXtalEnergyDep* xtalEDep = new McXtalEnergyDep();

            TVector3 moment1(depItr->second.moment1());
            TVector3 moment2(depItr->second.moment2());

            xtalEDep->initialize(depItr->second.getTotalEnergy(),
                                 depItr->second.getDirectEnergy(),
                                 moment1,
                                 moment2                          );

            mcIntHit->getXtalEnergyDepMap().Add(key, xtalEDep);
        }

        // Add the ROOT McIntegratingHit to the ROOT collection of McIntegratingHits
        m_mcEvt->addMcIntegratingHit(mcIntHit);
    }

    return sc;
}

StatusCode mcRootWriterAlg::writeMcTrajectories() 
{
    // Purpose and Method:  Retrieve the McTrajectory collection from the TDS 
    //  and fill the ROOT McTrajectory collection.
    // TDS Input:  EventModel::EVENT::McTrajectoryCol
    // Output:  ROOT McTrajectory collection

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Get the McPositionHits Collection from the TDS
    // Interestingly... there does not seem to be a string def for McTrajectoryCol! 
    SmartDataPtr<Event::McTrajectoryCol> trajectories(eventSvc(), "/Event/MC/TrajectoryCol");

    if (!trajectories) return sc;
    
    log << MSG::DEBUG;
    if( log.isActive()) 
        log.stream() << "Number of McTrajectory's in the event = " << trajectories->size();
    log << endreq;

    Event::McTrajectoryVector::const_iterator trajectory;
    for (trajectory = trajectories->begin(); trajectory != trajectories->end(); trajectory++ ) 
    {
        log << MSG::DEBUG;
        if( log.isActive()) (*trajectory)->fillStream(log.stream());
        log << endreq;
              
        McTrajectory *mcTrajectory = new McTrajectory();
        TRef ref = mcTrajectory;
        m_common.m_mcTrajectoryMap[(*trajectory)] = ref ;
        RootPersistence::convert(**trajectory,*mcTrajectory) ;  

        // Add the ROOT McPositionHit to the ROOT collection of McPositionHits
        m_mcEvt->addMcTrajectory(mcTrajectory);

        // Now loop through and set up the McTrajectoryPoints
        std::vector<Event::McTrajectoryPoint*> points = (*trajectory)->getPoints();
        std::vector<Event::McTrajectoryPoint*>::const_iterator pointIter;

        for(pointIter = points.begin(); pointIter != points.end(); pointIter++) 
        {
            // De-reference the McTrajectoryPoint pointer
            Event::McTrajectoryPoint* mcPointTds  = *pointIter;
            McTrajectoryPoint*        mcPointRoot = new McTrajectoryPoint();

            RootPersistence::convert(*mcPointTds, *mcPointRoot);

            mcTrajectory->addMcPoint(mcPointRoot);

            TRef pointRef = mcPointRoot;
            m_common.m_mcTrajectoryPointMap[mcPointTds] = pointRef;
        }
    }
    
    return sc;
}

StatusCode mcRootWriterAlg::writeMcTkrStrips() 
{
    // Purpose and Method:  Retrieve the McTrajectory collection from the TDS 
    //  and fill the ROOT McTrajectory collection.
    // TDS Input:  EventModel::EVENT::McTrajectoryCol
    // Output:  ROOT McTrajectory collection

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Get the McPositionHits Collection from the TDS
    // Interestingly... there does not seem to be a string def for McTrajectoryCol! 
    SmartDataPtr<Event::McTkrStripCol> stripCol(eventSvc(), "/Event/MC/StripCol");

    if (!stripCol) return sc;
    
    log << MSG::DEBUG;
    if( log.isActive()) 
        log.stream() << "Number of McTkrStrip's in the event = " << stripCol->size();
    log << endreq;

    ObjectVector<Event::McTkrStrip>::const_iterator stripIter;
    Event::McTrajectoryVector::const_iterator trajectory;
    for (stripIter = stripCol->begin(); stripIter != stripCol->end(); stripIter++ ) 
    {
        log << MSG::DEBUG;
        if( log.isActive()) (*stripIter)->fillStream(log.stream());
        log << endreq;

        const Event::McTkrStrip* stripTds  = *stripIter;
        McTkrStrip*              stripRoot = new McTkrStrip();

        RootPersistence::convert(*stripTds, *stripRoot);

        const Event::McTkrStrip::hitList& posHitListTds = stripTds->getHits();
        for(Event::McTkrStrip::hitList::const_iterator tdsIter = posHitListTds.begin();
            tdsIter != posHitListTds.end(); tdsIter++)
        {
            const Event::McPositionHit *mcPosHitTds = *tdsIter;
            McPositionHit *mcPosHitRoot = 0;
            if (mcPosHitTds != 0) {
                TRef ref;
                ref = m_common.m_mcPosHitMap[mcPosHitTds];
                mcPosHitRoot = (McPositionHit*)ref.GetObject();
            }
            stripRoot->addHit(mcPosHitRoot) ;
        }

        m_mcEvt->addMcTkrStrip(stripRoot);
    }

    return sc;
}



void mcRootWriterAlg::writeEvent() 
{
    // Purpose and Method:  Stores the McEvent data for this event in the ROOT
    //    tree.  The m_mcEvt object is cleared for the next event.
    m_rootIoSvc->fillTree("mc");

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

    m_rootIoSvc->closeFile("mc");

    return;
}

StatusCode mcRootWriterAlg::finalize()
{
    // ADDED FOR THE FILE HEADERS DEMO
    m_headersTool->writeMcHeader(m_mcTree->GetCurrentFile()) ;
    
    close();
    
    StatusCode sc = StatusCode::SUCCESS;
    //setFinalized(); No longer available in Gaudi v21r7
    return sc;
}

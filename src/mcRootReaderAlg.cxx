#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/IIncidentListener.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/MCEvent.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/MonteCarlo/McPositionHit.h"
#include "Event/MonteCarlo/McTrajectory.h"
#include "Event/Utilities/TimeStamp.h"

#include "mcRootData/McEvent.h"

#include "facilities/Util.h"
#include "commonData.h"
#include "RootConvert/Utilities/RootReaderUtil.h"
#include "RootIo/IRootIoSvc.h"

// ADDED FOR THE FILE HEADERS DEMO
#include "RootIo/FhTool.h"

// low level converters
#include <RootConvert/MonteCarlo/McPositionHitConvert.h>
#include <RootConvert/MonteCarlo/McTrajectoryConvert.h>
#include <RootConvert/Utilities/Toolkit.h>

#include <map>
#include <string>

/** @class mcRootReaderAlg
 * @brief Reads Monte Carlo data from a persistent ROOT file and stores the
 * the data in the TDS.
 *
 * @author Heather Kelly
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/mcRootReaderAlg.cxx,v 1.74.42.1.2.1 2010/04/07 13:25:17 heather Exp $
 */


class mcRootReaderAlg : public Algorithm, virtual public IIncidentListener
{	
public:
    
    mcRootReaderAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    /// Handles setup by opening ROOT file in read mode and creating a new TTree
    StatusCode initialize();
    
    /// Orchastrates reading from ROOT file and storing the data on the TDS 
    /// for each event
    StatusCode execute();

    /// handle "incidents"
    void handle(const Incident &inc) {
        if( inc.type()=="BeginEvent")beginEvent();
        else if(inc.type()=="EndEvent")endEvent();
    }
    
    void beginEvent() { }
    void endEvent();

    /// Closes the ROOT file and cleans up
    StatusCode finalize();

private:
    
    /// Retrieves event and run Id from the ROOT file and stores them on the TDS
    StatusCode readMcEvent();
    
    /// Retrieves McParticles from the ROOT file and stores them on the TDS
    StatusCode readMcParticles();
    
    /// Retrieves McPositionHits from the ROOT file and stores them on the TDS
    StatusCode readMcPositionHits();
    
    /// Retrieves McIntegratingHits from the ROOT file and stores them on the
    /// TDS
    StatusCode readMcIntegratingHits();
    
    /// Retrieves McTrajectorys from the ROOT file and stores them on the TDS
    StatusCode readMcTrajectories();

    /// Top-level Monte Carlo ROOT object
    McEvent *m_mcEvt;
//    /// name of the input ROOT file
    std::string m_fileName;
    /// Array of file names for TChain
    StringArrayProperty m_fileList;
    /// name of the Monte Carlo TTree stored in the ROOT file
    std::string m_treeName;
    /// Branch name for events
    std::string m_branchName;
    /// Option string which will be passed to McEvent::Clear
    std::string m_clearOption;
    // Branches to exclude from reading
    StringArrayProperty m_excludeBranchList;

    commonData m_common;

    /// Keep track of MC particles as we retrieve them from the ROOT file
    /// The id is the unique id assigned by ROOT for every TObject.
    std::map<UInt_t, Event::McParticle*> m_particleMap;

    IRootIoSvc*   m_rootIoSvc;
};

static const AlgFactory<mcRootReaderAlg>  Factory;
const IAlgFactory& mcRootReaderAlgFactory = Factory;


mcRootReaderAlg::mcRootReaderAlg(const std::string& name, 
                                 ISvcLocator* pSvcLocator) 
                                 : Algorithm(name, pSvcLocator), m_mcEvt(0)
{
    // Input pararmeters that may be set via the jobOptions file
    // Input ROOT file name, this will be overridden in RootIoSvc is provided a meta ROOT file for reading
    // Retain for backward compatibility, mcRootFileList is preferred
    declareProperty("mcRootFile",m_fileName="");

    // mcRootFileList to fill TChain, this will be overriden if RootIoSvc is provided a meta ROOT file for reading
    StringArrayProperty initList;
    std::vector<std::string> initVec;
    initList.setValue(initVec);
    declareProperty("mcRootFileList", m_fileList=initList);

    // ROOT TTree name
    declareProperty("mcTreeName", m_treeName="Mc");
    declareProperty("mcBranchName", m_branchName="McEvent");
    declareProperty("clearOption", m_clearOption="");
    declareProperty("ExcludeBranches",m_excludeBranchList=initList);
    
    m_particleMap.clear();


}

StatusCode mcRootReaderAlg::initialize()
{
    // Purpose and Method:  Called once before the run begins.
    //	  This method opens a new ROOT file and prepares for reading.
    
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    // Use the Job options service to set the Algorithm's parameters
    // This will retrieve parameters set in the job options file
    setProperties();

    m_rootIoSvc = 0 ;
    if ( service("RootIoSvc", m_rootIoSvc, true).isFailure() ){
        log << MSG::INFO << "Couldn't find the RootIoSvc!" << endreq;
        log << MSG::INFO << "Reading cannot continue" << endreq;
        // Now need RootIoSvc for both reading and writing, we cannot continue without it
        m_rootIoSvc = 0;
        return StatusCode::FAILURE;
    }   

    if ( (m_fileList.value().size() > 0) && ( !m_fileName.empty() )) {
        log << MSG::WARNING << "Both mcRootFile and mcRootFileList have "
            << "been specified, mcRootFile is deprecated, please use "
            << "mcRootFileList" << endreq;
         return StatusCode::FAILURE;
    } else if ( (m_fileList.value().size() == 0) && ( !m_fileName.empty() ) )
        m_rootIoSvc->appendFileList(m_fileList, m_fileName);
    else if (m_fileList.value().size() == 0)
        m_rootIoSvc->appendFileList(m_fileList, "mc.root");

    m_mcEvt = 0;
    m_common.m_mcEvt = m_mcEvt;

    // Set up new school system...
    // Use the name of this TTree (default "Mc") as key type 
    m_rootIoSvc->prepareRootInput("mc", m_treeName, m_branchName, 0, m_fileList);
    if (m_excludeBranchList.value().size() > 0) {
        std::vector<std::string>::const_iterator excludeListItr;
        for (excludeListItr = m_excludeBranchList.value().begin();
             excludeListItr != m_excludeBranchList.value().end();
             excludeListItr++ ) {
             std::string branchName = *excludeListItr;
             bool foundFlag  = m_rootIoSvc->setBranchStatus("mc",branchName,0);
             if (!foundFlag)
                 log << MSG::WARNING << "Did  not find any matching branch"
                     << " names for " << branchName << endreq;
             else
                 log << MSG::INFO << "Set BranchStatus to 0 (off) for " 
                     << "branch " << branchName << endreq;
        }

    }

    // use the incident service to register begin, end events
    IIncidentSvc* incsvc = 0;
    sc = service ("IncidentSvc", incsvc, true);

    if( sc.isFailure() ) return sc;

    incsvc->addListener(this, "BeginEvent", 100);
    incsvc->addListener(this, "EndEvent", 0);

 
    return sc;
    
}

StatusCode mcRootReaderAlg::execute()
{
    // Purpose and Method:  Called once per event.  This method calls
    //	 the appropriate methods to read data from the ROOT file and store
    //	 data on the TDS.
    
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    // Try reading the event this way... 
    // use name of TTree as key type
    m_mcEvt = dynamic_cast<McEvent*>(m_rootIoSvc->getNextEvent("mc"));

    if (!m_mcEvt) {
        // Do not fail if there was no MC data to read - this may be an Event Display run - where the user 
        // did not provide an MC input file
        log << MSG::INFO << "No MC Data Available" << endreq;
        return StatusCode::SUCCESS;
    }

    sc = readMcEvent();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to read top level McEvent" << endreq;
        return sc;
    }

    // Clear the root to TDS maps
    m_common.m_rootMcPartMap.clear();
    m_common.m_rootMcPosHitMap.clear();
    m_common.m_rootMcIntHitMap.clear();
    m_common.m_rootMcTrajectoryMap.clear();
    m_common.m_rootMcTrajectoryPointMap.clear();

    // Clear the TDS to root maps
    m_common.m_mcPartMap.clear();
    m_common.m_mcPosHitMap.clear();
    m_common.m_mcIntHitMap.clear();
    m_common.m_mcTrajectoryMap.clear();
    m_common.m_mcTrajectoryPointMap.clear();
    
    m_particleMap.clear();
    
    Event::MCEvent* mcEventObj = 
        SmartDataPtr<Event::MCEvent>(eventSvc(), EventModel::MC::Event);
    if (!mcEventObj) {
        log << "Failed to retrieve McEvent" << endreq;
        return sc;
    }
    
    sc = readMcParticles();
    if (sc.isFailure()) return sc;
    
    sc = readMcPositionHits();
    if (sc.isFailure()) return sc;
    
    sc = readMcIntegratingHits();
    if (sc.isFailure()) return sc;
    
    sc = readMcTrajectories();
    if (sc.isFailure()) return sc;

    return sc;
}

StatusCode mcRootReaderAlg::readMcEvent() {
    // Purpose and Method:  Retrieves the event and run ids from the ROOT file
    //	  Stores them on the TDS in the EventHeader
    // Input:  ROOT McEvent object
    // TDS Output:  EventModel::EventHeader
    
    MsgStream log(msgSvc(), name());
    
    StatusCode sc = StatusCode::SUCCESS;
    
    // Retrieve the Event data for this event
    SmartDataPtr<Event::EventHeader> evt(eventSvc(), EventModel::EventHeader);
    if (!evt) return sc;
    
    unsigned int eventIdTds = evt->event();
    unsigned int runIdTds = evt->run();
    
    unsigned int eventIdRoot = m_mcEvt->getEventId();
    unsigned int runIdRoot = m_mcEvt->getRunId();
    int sourceIdRoot = m_mcEvt->getSourceId();
    unsigned int sequenceRoot = m_mcEvt->getSequence();
    TimeStamp timeTds(m_mcEvt->getTimeStamp());
    
    // Check to see if the event and run ids have already been set.
    if (eventIdTds != eventIdRoot) evt->setEvent(eventIdRoot);
    if (runIdTds != runIdRoot) evt->setRun(runIdRoot);
    
    log << MSG::DEBUG << "Reading Event (run, event): (" << runIdRoot
        << ", " << eventIdRoot << ")" << endreq;

    SmartDataPtr<Event::MCEvent> mcEvt(eventSvc(), EventModel::MC::Event);
    if (!mcEvt) return sc;
    mcEvt->initialize(runIdRoot, sourceIdRoot, sequenceRoot, timeTds);

    mcEvt->setSourceName(m_mcEvt->getSourceName());
    
    return sc;
}

StatusCode mcRootReaderAlg::readMcParticles() {
    // Purpose and Method:  Retrieve the McParticle collection from the ROOT
    //	  file and fill the TDS McParticle collection
    // Input:  ROOT McParticle collection
    // TDS Output:  EventModel::EVENT::McParticleCol
    
    MsgStream log(msgSvc(), name());
    
    StatusCode sc = StatusCode::SUCCESS;
    
    TObjArray *particles = m_mcEvt->getMcParticleCol();
    if (!particles) return sc;
    TIter partIter(particles);
    
    // create the TDS location for the McParticle Collection
    Event::McParticleCol* pTdsCol = new Event::McParticleCol;
    sc = eventSvc()->registerObject(EventModel::MC::McParticleCol, pTdsCol);
    if (sc.isFailure()) {
        log << "Failed to register McParticle Collection" << endreq;
        return sc;
    }
    
    McParticle *pRoot;
    // Create the map of ROOT unique ids and TDS McParticle objects
    while ((pRoot = (McParticle*)partIter.Next()) != 0) {
        Event::McParticle *pTds = new Event::McParticle();
        m_particleMap[pRoot->GetUniqueID()] = pTds;
        m_common.m_rootMcPartMap[pRoot] = pTds;
        m_common.m_mcPartMap[pTds]      = pRoot;
    }
    
    // reset the iterator to the beginning of the list of McParticles
    partIter.Reset();
    
    // Now that the map is available, we initialize all of the TDS McParticles
    while ((pRoot = (McParticle*)partIter.Next())!=0) {
        
        Event::McParticle *pTds = m_particleMap[pRoot->GetUniqueID()];
        
        int idTds = pRoot->getParticleId();
        
        unsigned int statusBitsTds = pRoot->getStatusFlags();
        
        TLorentzVector initialMomRoot = pRoot->getInitialFourMomentum();
        
        CLHEP::HepLorentzVector initialMomTds(initialMomRoot.X(), initialMomRoot.Y(),
            initialMomRoot.Z(), initialMomRoot.T());
        
        TLorentzVector finalMomRoot = pRoot->getFinalFourMomentum();
        CLHEP::HepLorentzVector finalMomTds(finalMomRoot.X(), finalMomRoot.Y(), 
            finalMomRoot.Z(), finalMomRoot.T());
        
        TVector3 initPosRoot = pRoot->getInitialPosition();
        HepPoint3D initPosTds(initPosRoot.X(), 
            initPosRoot.Y(), initPosRoot.Z());
        
        TVector3 finalPosRoot = pRoot->getFinalPosition();
        HepPoint3D finalPosTds(finalPosRoot.X(), 
            finalPosRoot.Y(), finalPosRoot.Z());
        
        const McParticle *momRoot = pRoot->getMother();
        
        Event::McParticle *momTds = 0;
        
        if (momRoot != 0) {
            // The case where this is the primary particle
            if (momRoot->GetUniqueID() == pRoot->GetUniqueID()) {
                momTds = pTds;
            } else if (m_particleMap.find(momRoot->GetUniqueID()) != m_particleMap.end()) {
                momTds = m_particleMap[momRoot->GetUniqueID()];
            } else {
                log << MSG::INFO << "Failed to find the McParticle in the"
                    << " map!!" << endreq;
            }
        } else {
            static bool warn = false;
            if (!warn) {
                log << MSG::WARNING << "Cannot read in mother McParticles using "
                    << "\nSetting all McParticle mother pointers to"
                    << " point at themselves" << endreq;
                warn = true;
            }
            momTds = pTds;
        }

        std::string processTdsStr(pRoot->getProcess().Data());
        // Setup the TDS version fo the McParticle
        pTds->init(momTds, idTds, statusBitsTds, initialMomTds, 
                    finalMomTds, initPosTds, finalPosTds, processTdsStr );
        
        // Removed here some code to check root versions and prevent reading of 
        // daughter particles for versions older than root 3.04.02
        
        // Add the TDS McParticle to the TDS collection of McParticles
        pTdsCol->push_back(pTds);
        
    }
    
    return sc;
}

StatusCode mcRootReaderAlg::readMcPositionHits() {
    // Purpose and Method:  Retrieve the McPositionHit collection from the ROOT
    //	file and fill the TDS McPositionHit collection.
    // Input:  ROOT McPositionHit collectoin
    // TDS Output:  EventModel::MC::McPositionHitCol
    
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    const TObjArray *posHits = m_mcEvt->getMcPositionHitCol();
    if (!posHits) return sc;
    
    // this test and alternate return needed to feed to G4Generator, which wants to register its own collection
    if(  posHits->GetLast() ==-1)  return sc;

    TIter hitIter(posHits);
    
    // create the TDS location for the McParticle Collection
    Event::McPositionHitVector* pTdsCol = new Event::McPositionHitVector;
    sc = eventSvc()->registerObject(EventModel::MC::McPositionHitCol, pTdsCol);
    if (sc.isFailure()) {
        log << "Failed to register McPositionHit Collection" << endreq;
        return sc;
    }
    
    const McPositionHit *posHitRoot;
    while ((posHitRoot = (McPositionHit*)hitIter.Next())!=0) {
        
        Event::McPositionHit *posHitTds = new Event::McPositionHit();

        m_common.m_rootMcPosHitMap[posHitRoot] = posHitTds;
        m_common.m_mcPosHitMap[posHitTds]      = const_cast<McPositionHit*>(posHitRoot);
        
        RootPersistence::convert(*posHitRoot,*posHitTds) ;
        
        const McParticle* mcPartRoot = posHitRoot->getMcParticle();
        Event::McParticle *mcPartTds = 0;
        if (mcPartRoot != 0) {
            mcPartTds = m_particleMap[mcPartRoot->GetUniqueID()];
            posHitTds->setMcParticle(mcPartTds);
        }
        
        const McParticle *originRoot = posHitRoot->getOriginMcParticle();
        Event::McParticle *originTds = 0;
        if (originRoot != 0) {
            originTds = m_particleMap[originRoot->GetUniqueID()];
            posHitTds->setOriginMcParticle(originTds);
        }
        
        // add the McPositionHit to the TDS collection of McPositionHits
        pTdsCol->push_back(posHitTds);
    }
    
    return sc;
}

StatusCode mcRootReaderAlg::readMcIntegratingHits() {
    // Purpose and Method:  Retrieve the McIntegratingHit collection from the
    //	   ROOT file and fill the TDS McIntegratingHit collection.
    // Input:  ROOT McIntegratingHit collection
    // TDS Output:  EventModel::EVENT::McIntegratingHitCol
    
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    const TObjArray *intHits = m_mcEvt->getMcIntegratingHitCol();
    if (!intHits) return sc;

    // this test and alternate return needed to feed to G4Generator, which wants to register its own collection
    if(  intHits->GetLast() ==-1)  return sc;

    TIter hitIter(intHits);
    
    // create the TDS location for the McParticle Collection
    Event::McIntegratingHitVector* pTdsCol = new Event::McIntegratingHitVector;
    sc = eventSvc()->registerObject(EventModel::MC::McIntegratingHitCol, pTdsCol);
    if (sc.isFailure()) {
        log << "Failed to register McIntegratingHit" << endreq;
        return sc;
    }
    
    McIntegratingHit *intHitRoot;
    while ((intHitRoot = (McIntegratingHit*)hitIter.Next())!=0) {
        
        Event::McIntegratingHit *intHitTds = new Event::McIntegratingHit();

        m_common.m_rootMcIntHitMap[intHitRoot] = intHitTds;
        m_common.m_mcIntHitMap[intHitTds]      = intHitRoot;
        
        const VolumeIdentifier idRoot = intHitRoot->getVolumeId();
        idents::VolumeIdentifier idTds;
        idTds = RootPersistence::convert(idRoot) ;
        
        intHitTds->setVolumeID(idTds);
        
        //	  const McIntegratingHit::energyDepositMap mcPartMapRoot = intHitRoot->getItemizedEnergy();
        //	  McIntegratingHit::energyDepositMap::const_iterator rootMapIt;
        //	  log << MSG::DEBUG << "EnergyMap size: " << mcPartMapRoot.size() << endreq;
        // Can't seem to read this back in due to TRef problem in Root 3.02.03
        /*
        for (rootMapIt = mcPartMapRoot.begin(); rootMapIt != mcPartMapRoot.end(); rootMapIt++){
        McParticle* mcPartRoot = rootMapIt->first;
        Event::McParticle *mcPartTds = m_particleMap[mcPartRoot->GetUniqueID()];
        double e = rootMapIt->second;
        TVector3 posRoot = mcPartRoot->getFinalPosition();
        HepPoint3D posTds(posRoot.X(), posRoot.Y(), posRoot.Z());
        intHitTds->addEnergyItem(e, mcPartTds, posTds);
        }
        */
        
        double totalEnergyRoot = intHitRoot->getTotalEnergy();
        const TVector3 moment1Root = intHitRoot->getMoment1();
        HepPoint3D moment1Tds(moment1Root.X(), moment1Root.Y(), moment1Root.Z());
        const TVector3 moment2Root = intHitRoot->getMoment2();
        HepPoint3D moment2Tds(moment2Root.X(), moment2Root.Y(), moment2Root.Z());
        
        double energyArr[3] = { intHitRoot->getMcParticleEnergy(McIntegratingHit::PRIMARY),
            intHitRoot->getMcParticleEnergy(McIntegratingHit::ELECTRON),
            intHitRoot->getMcParticleEnergy(McIntegratingHit::POSITRON) };
        
        intHitTds->setEnergyItems(totalEnergyRoot, energyArr, moment1Tds, moment2Tds);
        
        // Add the TDS McIntegratingHit to the TDS McIntegratingHit collection
        pTdsCol->push_back(intHitTds);
    }
    
    return sc;
}

StatusCode mcRootReaderAlg::readMcTrajectories() 
{
    // Purpose and Method:  Retrieve the McTrajectory collection from the ROOT
    //	file and fill the TDS McTrajectory collection.
    //  NOTE: This collection is usually NOT present
    // Input:  ROOT McTrajectory collection
    // TDS Output:  EventModel::MC::McTrajectoryCol
    
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    const TObjArray *trajectories = m_mcEvt->getMcTrajectoryCol();
    if (!trajectories) return sc;

    // this test and alternate return needed to feed to G4Generator, which wants to register its own collection
    if(  trajectories->GetLast() ==-1)  return sc;

    TIter hitIter(trajectories);
    
    // create the TDS location for the McTrajectory Collection
    Event::McTrajectoryCol* pTdsCol = new Event::McTrajectoryCol();  
    eventSvc()->registerObject(EventModel::MC::McTrajectoryCol,pTdsCol);
    if (sc.isFailure()) {
        log << "Failed to register McTrajectory Collection" << endreq;
        return sc;
    }
    
    // Loop through the root trajectories and do the conversion
    for(int rootIdx = 0; rootIdx < trajectories->GetEntries(); rootIdx++)
    {
        McTrajectory* trajectoryRoot = (McTrajectory*)trajectories->At(rootIdx);

        Event::McTrajectory* trajectoryTds = new Event::McTrajectory();

        m_common.m_rootMcTrajectoryMap[trajectoryRoot] = trajectoryTds;
        m_common.m_mcTrajectoryMap[trajectoryTds]      = trajectoryRoot;
        
        RootPersistence::convert(*trajectoryRoot,*trajectoryTds);

        // Now loop through and set up the McTrajectoryPoints
        const TObjArray* points = trajectoryRoot->getMcPointCol();

        for(int idx=0; idx < points->GetEntries(); idx++)
        {
            // Pointer to the McTrajectoryPoint
            McTrajectoryPoint* pointRoot = (McTrajectoryPoint*)(points->At(idx));

            // Pointer to the TDS McTrajectoryPoint
            Event::McTrajectoryPoint* pointTds = new Event::McTrajectoryPoint();

            // Store in map
            m_common.m_rootMcTrajectoryPointMap[pointRoot] = pointTds;
            m_common.m_mcTrajectoryPointMap[pointTds]      = pointRoot;

            RootPersistence::convert(*pointRoot, *pointTds);

            trajectoryTds->addPoint(pointTds);
        }
        
        // add the McTrajectory to the TDS collection of McTrajectories
        pTdsCol->push_back(trajectoryTds);
    }

    return sc;
}

void mcRootReaderAlg::endEvent() {
    if (m_mcEvt) m_mcEvt->Clear(m_clearOption.c_str());
    m_mcEvt = 0;
}


StatusCode mcRootReaderAlg::finalize()
{
    StatusCode sc = StatusCode::SUCCESS;
    setFinalized();
    return sc;
}

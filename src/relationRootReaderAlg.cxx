#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/DigiEvent.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McTrajectory.h"
#include "Event/MonteCarlo/McPositionHit.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/MonteCarlo/McRelTableDefs.h"
#include "Event/Digi/AcdDigi.h"
#include "Event/Digi/CalDigi.h"
#include "Event/Digi/TkrDigi.h"

#include "Event/RelTable/RelTable.h"

#include "idents/CalXtalId.h"
#include "idents/TowerId.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TChainIndex.h"
#include "TDirectory.h"
#include "TObjArray.h"
#include "TCollection.h"  // Declares TIter
#ifdef WIN32
#include "TSystem.h" // To get TreePlayer loaded
#endif

#include "mcRootData/McEvent.h"
#include "digiRootData/DigiEvent.h"
#include "reconRootData/ReconEvent.h"
#include "commonRootData/RelTable.h"

#include "facilities/Util.h"

#include "commonData.h"
#include "RootIo/IRootIoSvc.h"
#include "RootConvert/Utilities/RootReaderUtil.h"

/** @class relationRootReaderAlg
 * @brief Reads Digitization data from a persistent ROOT file and stores the
 * the data in the TDS.
 *
 * @author Heather Kelly
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/relationRootReaderAlg.cxx,v 1.24 2007/02/15 19:30:06 usher Exp $
 */

class relationRootReaderAlg : public Algorithm
{	
public:
    
    relationRootReaderAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    /// Handles setup by opening ROOT file in read mode and creating a new TTree
    StatusCode initialize();
   
    /// Orchastrates reading from ROOT file and storing the data on the TDS for each event
    StatusCode execute();
    
    /// Closes the ROOT file and cleans up
    StatusCode finalize();
            
private:

    /// Creates TDS tables
    StatusCode createTDSTables();

    /// Reads top-level DigiEvent
    StatusCode readRelations();
    

    /// Handles converting TkrDigi relations back to TDS
    StatusCode readTkrDigiRelations(Relation*);
    /// Handles converting CalDigi relations back to TDS
    StatusCode readCalDigiRelations(Relation*);
    /// Handles converting track relations back to TDS
    StatusCode readTkrTrackRelations(Relation*);
    /// Handles converting vertex relations back to TDS
    StatusCode readTkrVertexRelations(Relation*);
    /// Handles converting McParticle <--> McTrajectory relations back to TDS
    StatusCode readMcParticleRelations(Relation*);
    /// Handles converting McTrajectoryPoint to hit relations back to TDS
    StatusCode readMcTrajectoryPointRelations(Relation*);

    //StatusCode readTkrDigiRelations();  
    //StatusCode readCalDigiRelations();

    /// Closes the ROOT file
    void close();
   
    /// ROOT file pointer
    TFile *m_relFile;
    /// ROOT tree pointer
    TChain *m_relTree;
    /// Top-level Monte Carlo ROOT object
    RelTable *m_relTab;
    /// name of the output ROOT file
    std::string m_fileName;
    /// List of files
    StringArrayProperty m_fileList;
    /// name of the Monte Carlo TTree stored in the ROOT file
    std::string m_treeName;
    /// Stores number of events available in the input ROOT TTree
    Long64_t m_numEvents;

    commonData m_common;
    IRootIoSvc*   m_rootIoSvc;

/// typedefs for tables
    typedef Event::RelTable<Event::TkrDigi,Event::McPositionHit>    TkrDigiRelTab;
    typedef Event::Relation<Event::TkrDigi,Event::McPositionHit>    TkrDigiRelType;
    typedef Event::RelTable<Event::CalDigi,Event::McIntegratingHit> CalDigiRelTab;
    typedef Event::Relation<Event::CalDigi,Event::McIntegratingHit> CalDigiRelType;

    /// Internal pointer to TDS TkrDigi table
    TkrDigiRelTab*            m_tkrDigiRelTab;
    /// Internal pointer to tDS CalDigi table
    CalDigiRelTab*            m_calDigiRelTab;

    /// Internal pointer to McParticle object list
    Event::McPartToTrajectoryTabList* m_mcPartToTrajList;

    /// Internal pointer to McTrajectoryPoint to Pos hit list
    Event::McPointToPosHitTabList*    m_mcPointToPosHitList;

    /// Internal pointer to McTrajectoryPoint to Pos hit list
    Event::McPointToIntHitTabList*    m_mcPointToIntHitList;
};

static const AlgFactory<relationRootReaderAlg>  Factory;
const IAlgFactory& relationRootReaderAlgFactory = Factory;


relationRootReaderAlg::relationRootReaderAlg(const std::string& name, ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator)
{
    // Input pararmeters that may be set via the jobOptions file
    // Input ROOT file name
    declareProperty("rootFile",m_fileName="");
    StringArrayProperty initList;
    std::vector<std::string> initVec;
    initVec.push_back("relation.root");
    initList.setValue(initVec);
    declareProperty("relationRootFileList", m_fileList=initList);
    initVec.clear();
    // Input TTree name
    declareProperty("treeName", m_treeName="Relations");

    m_mcPartToTrajList    = 0;
    m_mcPointToPosHitList = 0;
    m_mcPointToIntHitList = 0;
}

StatusCode relationRootReaderAlg::initialize()
{
    // Purpose and Method:  Called once before the run begins.  This method
    //    opens a new ROOT file and prepares for reading.

    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    // Use the Job options service to set the Algorithm's parameters
    // This will retrieve parameters set in the job options file
    setProperties();

     m_rootIoSvc = 0 ;
   if ( service("RootIoSvc", m_rootIoSvc, true).isFailure() ){
        log << MSG::INFO << "Couldn't find the RootIoSvc!" << endreq;
        log << MSG::INFO << "Event loop will not terminate gracefully" << endreq;
        m_rootIoSvc = 0;
        //return StatusCode::FAILURE;
    }

    facilities::Util::expandEnvVar(&m_fileName);

#ifdef WIN32
	gSystem->Load("libTreePlayer.dll");
#endif  
   
    // Save the current directory for the ntuple writer service
    TDirectory *saveDir = gDirectory;   
    m_relTree = new TChain(m_treeName.c_str());

    
    // add root files to TChain and check if files exist
    StatusCode openSc= RootPersistence::addFilesToChain(m_relTree, m_fileName, m_fileList, log);
    if (openSc.isFailure()) {
      return  openSc;
    }
    

    m_relTab = 0;
    m_relTree->SetBranchAddress("RelTable", &m_relTab);

    m_numEvents = m_relTree->GetEntries();
    if (m_rootIoSvc) {
        m_rootIoSvc->setRootEvtMax(m_numEvents);
        if(!m_relTree->GetTreeIndex()) {
          //  log << MSG::INFO << "Input file does not contain new style index, rebuilding" << endreq;
          //  m_relTree->BuildIndex("m_runId", "m_eventId");
          TChainIndex *chInd = new TChainIndex(m_relTree, "m_runId", "m_eventId");
          
          // could be zombie if files in chain are not in runId order
          if (!chInd->IsZombie())
              m_relTree->SetTreeIndex(chInd);
        }
        m_rootIoSvc->registerRootTree(m_relTree);
    }

    saveDir->cd();
    return sc;
    
}

StatusCode relationRootReaderAlg::execute()
{
    // Purpose and Method:  Called once per event.  This method calls
    //   the appropriate methods to read data from the ROOT file and store
    //   data on the TDS.

    MsgStream log(msgSvc(), name());

    StatusCode sc = StatusCode::SUCCESS;
    
    static Long64_t evtId = 0;
    Long64_t readInd;
    int numBytes;
	std::pair<int,int> runEventPair = (m_rootIoSvc) ? m_rootIoSvc->runEventPair() : std::pair<int,int>(-1,-1);
	
	if ((m_rootIoSvc) && (m_rootIoSvc->useIndex())) {
		readInd = m_rootIoSvc->index();
	} else if ((m_rootIoSvc) && (m_rootIoSvc->useRunEventPair())) {
		int run = runEventPair.first;
		int evt = runEventPair.second;
                // returns -1 if index does not exist
		readInd = m_relTree->GetEntryNumberWithIndex(run, evt);
	} else {
		readInd = evtId;
	}

    if ((readInd >= m_numEvents) || (readInd < 0)) {
        log << MSG::WARNING << "Requested index is out of bounds - no relation data retrieved" << endreq;
        return StatusCode::SUCCESS;
    }

    if (m_rootIoSvc) m_rootIoSvc->setActualIndex(readInd);

    numBytes = m_relTree->GetEvent(readInd);
	
    if ((numBytes <= 0) || (!m_relTab)) {
        log << MSG::WARNING << "Failed to Relational Table" << endreq;
        return StatusCode::SUCCESS;
    }


    sc = createTDSTables();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to create TDS relational tables" << endreq;
        return sc;
    }

    sc = readRelations();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to read Relational Table" << endreq;
        return sc;
    }


    /// Register our relational tables in the TDS (and turn over ownership)
    sc = eventSvc()->registerObject(EventModel::Digi::TkrDigiHitTab, m_tkrDigiRelTab->getAllRelations() );
    sc = eventSvc()->registerObject(EventModel::Digi::CalDigiHitTab, m_calDigiRelTab->getAllRelations() );
    sc = eventSvc()->registerObject("/Event/MC/McPartToTrajectory",  m_mcPartToTrajList );
    sc = eventSvc()->registerObject("/Event/MC/McPointToPosHit",     m_mcPointToPosHitList );
    sc = eventSvc()->registerObject("/Event/MC/McPointToIntHit",     m_mcPointToIntHitList );

    m_relTab->Clear();

	evtId = readInd+1;

    m_common.clear();
    
    return sc;
}

StatusCode relationRootReaderAlg::createTDSTables()
{
    // Purpose and Method:  Creates the TDS relational tables for output
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    m_tkrDigiRelTab = new TkrDigiRelTab();
    m_tkrDigiRelTab->init();

    m_calDigiRelTab = new CalDigiRelTab();
    m_calDigiRelTab->init();

    m_mcPartToTrajList    = new Event::RelationList<Event::McParticle,Event::McTrajectory>;
    m_mcPointToPosHitList = new Event::RelationList<Event::McTrajectoryPoint,Event::McPositionHit>;
    m_mcPointToIntHitList = new Event::RelationList<Event::McTrajectoryPoint,Event::McIntegratingHit>;

    return sc;
}

StatusCode relationRootReaderAlg::readRelations() {
    // Purpose and Method:  "Reads" the individual relations and, depending upon type,
    //                      calls the appropriate method for converting back to TDS relation
    MsgStream log(msgSvc(), name());

    StatusCode sc = StatusCode::SUCCESS;

    // Get pointer to the root relations
    const TObjArray* relations = m_relTab->getRelationTable();

    // Loop over these relations
    for(int idx = 0; idx < relations->GetEntries(); idx++)
    {
        TObject*  relObj   = relations->At(idx);
        Relation* relation = dynamic_cast<Relation*>(relObj);

        const TObject* key = relation->getKey();

        // Dynamic cast the key back to its original object and decide how to read the relation
        if      (const TkrDigi*           digi    = dynamic_cast<const TkrDigi*>(key))           readTkrDigiRelations(relation);
        else if (const CalDigi*           digi    = dynamic_cast<const CalDigi*>(key))           readCalDigiRelations(relation);
        else if (const TkrVertex*         vert    = dynamic_cast<const TkrVertex*>(key))         readTkrVertexRelations(relation);
        else if (const McParticle*        mcPart  = dynamic_cast<const McParticle*>(key))        readMcParticleRelations(relation);
        else if (const McTrajectoryPoint* mcPoint = dynamic_cast<const McTrajectoryPoint*>(key)) readMcTrajectoryPointRelations(relation);
        else
        {
            log << MSG::WARNING << "Unrecognized Relation key" << endreq;
        }
    }
    
    return sc;
}

StatusCode relationRootReaderAlg::readTkrDigiRelations(Relation* relation) 
{
    // Purpose and Method:  Converts a root TkrDigi/McPositionHit relation to TDS version
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    const TkrDigi*        digiRoot  = dynamic_cast<const TkrDigi*>(relation->getKey());
    TObject*              mcHitRoot = relation->getValueCol().At(0);

    Event::TkrDigi*       digiTds   = 0;
    Event::McPositionHit* mcHitTds  = 0;

    // Look up TDS to root relation for digis
    if (m_common.m_rootTkrDigiMap.find(digiRoot) != m_common.m_rootTkrDigiMap.end()) {
        digiTds = const_cast<Event::TkrDigi*>(m_common.m_rootTkrDigiMap[digiRoot]);
    } else {
        log << MSG::WARNING << "Could not located TkrDigi TDS/ROOT pair" << endreq;
    }

    // Look up TDS to root relation for position hits
    if (m_common.m_rootMcPosHitMap.find(mcHitRoot) != m_common.m_rootMcPosHitMap.end()) {
        mcHitTds = const_cast<Event::McPositionHit*>(m_common.m_rootMcPosHitMap[mcHitRoot]);
    } else {
        log << MSG::WARNING << "Could not located McPositionHit TDS/ROOT pair" << endreq;
    }

    // Make the event relation here
    typedef Event::Relation<Event::TkrDigi,Event::McPositionHit> relType;
    relType* rel = new relType(digiTds, mcHitTds);
    m_tkrDigiRelTab->addRelation(rel);

    return sc;
}


StatusCode relationRootReaderAlg::readCalDigiRelations(Relation* relation) 
{
    // Purpose and Method:  Converts a root CalDigi/McIntegratingHit relation to TDS version

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    const CalDigi*           digiRoot  = dynamic_cast<const CalDigi*>(relation->getKey());
    TObject*                 mcHitRoot = relation->getValueCol().At(0);

    Event::CalDigi*          digiTds   = 0;
    Event::McIntegratingHit* mcHitTds  = 0;

    // Look up TDS to root relation for digis
    if (m_common.m_rootCalDigiMap.find(digiRoot) != m_common.m_rootCalDigiMap.end()) {
        digiTds = const_cast<Event::CalDigi*>(m_common.m_rootCalDigiMap[digiRoot]);
    } else {
        log << MSG::WARNING << "Could not located CalDigi TDS/ROOT pair" << endreq;
    }

    // Look up TDS to root relation for position hits
    if (m_common.m_rootMcIntHitMap.find(mcHitRoot) != m_common.m_rootMcIntHitMap.end()) {
        mcHitTds = const_cast<Event::McIntegratingHit*>(m_common.m_rootMcIntHitMap[mcHitRoot]);
    } else {
        log << MSG::WARNING << "Could not located McIntegratingHit TDS/ROOT pair" << endreq;
    }

    // Make the event relation here
    typedef Event::Relation<Event::CalDigi,Event::McIntegratingHit> relType;
    relType* rel = new relType(digiTds, mcHitTds);
    m_calDigiRelTab->addRelation(rel);

    return sc;
}

/*
StatusCode relationRootReaderAlg::readTkrTrackRelations(Relation* relation) 
{
    // Purpose and Method:  Converts a root TkrPatCand/TkrFitTrack relation to TDS version

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    const TkrCandTrack*     tkrCandRoot  = dynamic_cast<const TkrCandTrack*>(relation->getKey());

    const TRefArray&        refArray     = relation->getValueCol();
    TObject*                tkrTrackRoot = refArray.At(0);
    Event::TkrPatCand*      tkrCandTds   = 0;
    Event::TkrFitTrackBase* tkrTrackTds  = 0;

    // Look up TDS to root relation for candidate tracks
    if (m_common.m_rootTkrCandMap.find(tkrCandRoot) != m_common.m_rootTkrCandMap.end()) {
        tkrCandTds = const_cast<Event::TkrPatCand*>(m_common.m_rootTkrCandMap[tkrCandRoot]);
    } else {
        log << MSG::WARNING << "Could not located TkrPatCand TDS/ROOT pair" << endreq;
    }

    // Look up TDS to root relation for candidate tracks
    if (m_common.m_rootTkrTrackMap.find(tkrTrackRoot) != m_common.m_rootTkrTrackMap.end()) {
        tkrTrackTds = const_cast<Event::TkrFitTrackBase*>(m_common.m_rootTkrTrackMap[tkrTrackRoot]);
    } else {
        log << MSG::WARNING << "Could not located TkrFitTrackBase TDS/ROOT pair" << endreq;
    }

    // Make the event relation here
    Event::TkrFitTrackRel* rel = new Event::TkrFitTrackRel(tkrCandTds, tkrTrackTds);
    m_trackRelTab->addRelation(rel);
    
    return sc;
}
*/
StatusCode relationRootReaderAlg::readTkrVertexRelations(Relation* relation) 
{
    // Purpose and Method:  Converts a root TkrVertex/TkrFitTrack relation to TDS version

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    const TkrVertex*  vertexRoot    = dynamic_cast<const TkrVertex*>(relation->getKey());

    const TRefArray&  refArray     = relation->getValueCol();
    TObject*          tkrTrackRoot = refArray.At(0);
    Event::TkrVertex* vertexTds    = 0;
    Event::TkrTrack*  tkrTrackTds  = 0;

    // Look up TDS to root relation for vertices
    if (m_common.m_rootTkrVertexMap.find(vertexRoot) != m_common.m_rootTkrVertexMap.end()) {
        vertexTds = const_cast<Event::TkrVertex*>(m_common.m_rootTkrVertexMap[vertexRoot]);
    } else {
        log << MSG::WARNING << "Could not located TkrVertex TDS/ROOT pair" << endreq;
    }

    // Look up TDS to root relation for candidate tracks
    if (m_common.m_rootTkrTrackMap.find(tkrTrackRoot) != m_common.m_rootTkrTrackMap.end()) {
        tkrTrackTds = const_cast<Event::TkrTrack*>(m_common.m_rootTkrTrackMap[tkrTrackRoot]);
    } else {
        log << MSG::WARNING << "Could not located TkrTrack TDS/ROOT pair" << endreq;
    }

    return sc;
}


StatusCode relationRootReaderAlg::readMcParticleRelations(Relation* relation) 
{
    // Purpose and Method:  Converts a root McParticle <--> McTrajectory relation to TDS version
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    const McParticle*    mcPartRoot = dynamic_cast<const McParticle*>(relation->getKey());
    TObject*             mcTrajRoot = relation->getValueCol().At(0);

    Event::McParticle*   mcPartTds  = 0;
    Event::McTrajectory* mcTrajTds  = 0;

    // Look up TDS to root relation for McParticles
    if (m_common.m_rootMcPartMap.find(mcPartRoot) != m_common.m_rootMcPartMap.end()) {
        mcPartTds = const_cast<Event::McParticle*>(m_common.m_rootMcPartMap[mcPartRoot]);
    } else {
        log << MSG::WARNING << "Could not located McParticle TDS/ROOT pair" << endreq;
    }

    // Look up TDS to root relation for position hits
    if (m_common.m_rootMcTrajectoryMap.find(mcTrajRoot) != m_common.m_rootMcTrajectoryMap.end()) {
        mcTrajTds = const_cast<Event::McTrajectory*>(m_common.m_rootMcTrajectoryMap[mcTrajRoot]);
    } else {
        log << MSG::WARNING << "Could not located McTrajectory TDS/ROOT pair" << endreq;
    }

    // Make the event relation here
    Event::McPartToTrajectoryRel* rel = new Event::McPartToTrajectoryRel(mcPartTds, mcTrajTds);
    m_mcPartToTrajList->push_back(rel);

    return sc;
}

StatusCode relationRootReaderAlg::readMcTrajectoryPointRelations(Relation* relation) 
{
    // Purpose and Method:  Converts a root TkrDigi/McPositionHit relation to TDS version
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    const McTrajectoryPoint*  mcPointRoot  = dynamic_cast<const McTrajectoryPoint*>(relation->getKey());
    TObject*                  mcHitRoot = relation->getValueCol().At(0);

    Event::McTrajectoryPoint* mcPointTds   = 0;
    
    Event::McPositionHit* mcHitTds  = 0;

    // Look up TDS to root relation for digis
    if (m_common.m_rootMcTrajectoryPointMap.find(mcPointRoot) != m_common.m_rootMcTrajectoryPointMap.end()) {
        mcPointTds = const_cast<Event::McTrajectoryPoint*>(m_common.m_rootMcTrajectoryPointMap[mcPointRoot]);
    } else {
        log << MSG::WARNING << "Could not located McTrajectoryPoint TDS/ROOT pair" << endreq;
    }

    // Look up TDS to root relation for position hits
    if (m_common.m_rootMcPosHitMap.find(mcHitRoot) != m_common.m_rootMcPosHitMap.end()) 
    {
        Event::McPositionHit* mcHitTds = const_cast<Event::McPositionHit*>(m_common.m_rootMcPosHitMap[mcHitRoot]);
        Event::McPointToPosHitRel* rel = new Event::McPointToPosHitRel(mcPointTds, mcHitTds);
        m_mcPointToPosHitList->push_back(rel);
    }
    else if (m_common.m_rootMcIntHitMap.find(mcHitRoot) != m_common.m_rootMcIntHitMap.end())
    {
        Event::McIntegratingHit* mcHitTds = const_cast<Event::McIntegratingHit*>(m_common.m_rootMcIntHitMap[mcHitRoot]);
        Event::McPointToIntHitRel* rel = new Event::McPointToIntHitRel(mcPointTds, mcHitTds);
        m_mcPointToIntHitList->push_back(rel);
    } 
    else 
    {
        log << MSG::WARNING << "Could not located McPositionHit/McIntegrating TDS/ROOT pair" << endreq;
    }

    return sc;
}


void relationRootReaderAlg::close() 
{
    // Purpose and Method:  Closes the ROOT file at the end of the run.

   // HMK skip deleting chain due to seg fault with TChainIndex
   // if (m_relTree) {
   //     delete m_relTree;
   //     m_relTree = 0;
   // }
}

StatusCode relationRootReaderAlg::finalize()
{
    close();
    
    StatusCode sc = StatusCode::SUCCESS;
    setFinalized();
    return sc;
}


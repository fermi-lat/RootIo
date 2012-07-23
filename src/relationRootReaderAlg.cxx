#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/IIncidentListener.h"


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
 * Note this algorithm should be run after all other root reading algorithms are run so that the objects pointed to by
 * the relational table exist when the relations are read in.
 *
 * @author Heather Kelly
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/relationRootReaderAlg.cxx,v 1.47 2011/12/12 20:55:41 heather Exp $
 */

class relationRootReaderAlg : public Algorithm, virtual public IIncidentListener
{	
public:
    
    relationRootReaderAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    /// Handles setup by opening ROOT file in read mode and creating a new TTree
    StatusCode initialize();
   
    /// Orchastrates reading from ROOT file and storing the data on the TDS for each event
    StatusCode execute();
    
    /// handle "incidents"
    void handle(const Incident &inc) {
        if( inc.type()=="BeginEvent")beginEvent();
        else if(inc.type()=="EndEvent")endEvent();
    }

    void beginEvent() {};
    void endEvent();

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
    /// Handles converting CalXtalRecData to CalCluster relations back to TDS
    StatusCode readCalXtalToClusterRelations(Relation*);
    /// Handles converting track relations back to TDS
    StatusCode readTkrTrackRelations(Relation*);
    /// Handles converting vertex relations back to TDS
    StatusCode readTkrVertexRelations(Relation*);
    /// Handles converting McParticle <--> McTrajectory relations back to TDS
    StatusCode readMcParticleRelations(Relation*);
    /// Handles converting McTrajectoryPoint to hit relations back to TDS
    StatusCode readMcTrajectoryPointRelations(Relation*);

    /// Closes the ROOT file
    void close();
   
    /// Top-level Monte Carlo ROOT object
    RelTable *m_relTab;
//    /// name of the input ROOT file
    std::string m_fileName;
    /// List of files
    StringArrayProperty m_fileList;
    /// name of the Monte Carlo TTree stored in the ROOT file
    std::string m_treeName;
    /// Stores number of events available in the input ROOT TTree
    /// Branch name for events
    std::string m_branchName;
    /// Option string which will be passed to McEvent::Clear
    std::string m_clearOption;
    /// Branch Exclusion List
    StringArrayProperty m_excludeBranchList;

    Long64_t m_numEvents;

    commonData m_common;
    IRootIoSvc*   m_rootIoSvc;

    bool m_terminateOnReadError;

/// typedefs for tables
    typedef Event::RelTable<Event::TkrDigi,Event::McPositionHit>         TkrDigiToPosHitTab;
    typedef Event::Relation<Event::TkrDigi,Event::McPositionHit>         TkrDigiToPosHitRel;
    typedef Event::RelationList<Event::TkrDigi,Event::McPositionHit>     TkrDigiToPosHitTabList;
    typedef Event::RelTable<Event::CalXtalRecData,Event::CalCluster>     CalClusterToXtalRecTab;
    typedef Event::Relation<Event::CalXtalRecData,Event::CalCluster>     CalClusterToXtalRecRel;
    typedef Event::RelationList<Event::CalXtalRecData,Event::CalCluster> CalClusterToXtalRecTabList;
    typedef Event::RelTable<Event::CalDigi,Event::McIntegratingHit>      CalDigiToIntHitTab;
    typedef Event::Relation<Event::CalDigi,Event::McIntegratingHit>      CalDigiToIntHitRel;
    typedef Event::RelationList<Event::CalDigi,Event::McIntegratingHit>  CalDigiToIntHitTabList;

    /// Internal pointer to TDS TkrDigi table
    TkrDigiToPosHitTabList*           m_tkrDigiToPosHitList;
    /// Internal pointer to TDS CalDigi table
    CalDigiToIntHitTabList*           m_calDigiToIntHitList;
    /// Internal pointer to the TDS CalCluster table
    CalClusterToXtalRecTabList*       m_calClusterToXtalRecList;
    /// Internal pointer to McParticle object list
    Event::McPartToTrajectoryTabList* m_mcPartToTrajList;

    /// Internal pointer to McTrajectoryPoint to Pos hit list
    Event::McPointToPosHitTabList*    m_mcPointToPosHitList;

    /// Internal pointer to McTrajectoryPoint to Pos hit list
    Event::McPointToIntHitTabList*    m_mcPointToIntHitList;

    /// Define a map of function pointers to handle conversions
    typedef StatusCode(relationRootReaderAlg::*ConvFuncPtr)(Relation*);
    std::map<std::string, ConvFuncPtr> m_convFuncMap;
    ConvFuncPtr                        m_currentConverter;
    std::string                        m_currentRelType;
};

//static const AlgFactory<relationRootReaderAlg>  Factory;
//const IAlgFactory& relationRootReaderAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(relationRootReaderAlg);

relationRootReaderAlg::relationRootReaderAlg(const std::string& name, ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator), m_relTab(0), m_branchName("RelTable")
{
    // Input pararmeters that may be set via the jobOptions file
    // Input ROOT file name, for backward compatibility with older JO files, relationRootFileList is preferred
    // This will be overridden if RootIoSvc is provided a meta ROOT file for reading
    declareProperty("rootFile",m_fileName="");

    // Provide a list of input files to init TChain, this will be overridden if RootIoSvc is provided a
    // meta ROOT file for reading
    StringArrayProperty initList;
    std::vector<std::string> initVec;
    initList.setValue(initVec);
    declareProperty("relationRootFileList", m_fileList=initList);
    initVec.clear();
    // Input TTree name
    declareProperty("relTreeName", m_treeName="Relations");
    declareProperty("relBranchName", m_branchName="RelTable");
    declareProperty("clearOption", m_clearOption="");
    declareProperty("ExcludeBranches",m_excludeBranchList=initList);

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
        log << MSG::INFO << "Reading cannot continue" << endreq;
        m_rootIoSvc = 0;
        return StatusCode::FAILURE;
    }

    if ( (m_fileList.value().size() > 0) && ( !m_fileName.empty() )) {
        log << MSG::WARNING << "Both relationRootFile and relationRootFileList "
            << " have been specified, relationRootFile is deprecated, "
            << "please use "
            << "relationRootFileList" << endreq;
         return StatusCode::FAILURE;
    } else if ( (m_fileList.value().size() == 0) && ( !m_fileName.empty() ) )
        m_rootIoSvc->appendFileList(m_fileList, m_fileName);
    else if (m_fileList.value().size() == 0)
        m_rootIoSvc->appendFileList(m_fileList, "relations.root");

    m_terminateOnReadError = m_rootIoSvc->terminateOnReadError();


    // Set up new school system...
    // Use treeName as key type
    m_rootIoSvc->prepareRootInput("rel", m_treeName, m_branchName, 0, m_fileList);

    if (m_excludeBranchList.value().size() > 0) {
        std::vector<std::string>::const_iterator excludeListItr;
        for (excludeListItr = m_excludeBranchList.value().begin();
             excludeListItr != m_excludeBranchList.value().end();
             excludeListItr++ ) {
             std::string branchName = *excludeListItr;
             bool foundFlag  = m_rootIoSvc->setBranchStatus("rel",branchName,0);
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

    // Set up the map between relation key types and the actual conversion routine
    m_convFuncMap["TkrDigi"]           = &relationRootReaderAlg::readTkrDigiRelations;
    m_convFuncMap["CalDigi"]           = &relationRootReaderAlg::readCalDigiRelations;
    m_convFuncMap["CalXtalRecData"]    = &relationRootReaderAlg::readCalXtalToClusterRelations;
    m_convFuncMap["TkrVertex"]         = &relationRootReaderAlg::readTkrVertexRelations;
    m_convFuncMap["McParticle"]        = &relationRootReaderAlg::readMcParticleRelations;
    m_convFuncMap["McTrajectoryPoint"] = &relationRootReaderAlg::readMcTrajectoryPointRelations;

    m_currentConverter = 0;
    m_currentRelType   = "";

    return sc;
    
}

StatusCode relationRootReaderAlg::execute()
{
    // Purpose and Method:  Called once per event.  This method calls
    //   the appropriate methods to read data from the ROOT file and store
    //   data on the TDS.

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS; 
    
    // Try reading the event this way... 
    // use treeName as key type
    m_relTab = dynamic_cast<RelTable*>(m_rootIoSvc->getNextEvent("rel"));

    if (!m_relTab) {
        if (m_terminateOnReadError) {
            log << MSG::ERROR << "Failed to read in Relation Data" << endreq;
            return StatusCode::FAILURE;
         }
         // Do not fail if there was no Relation data to read - this may be an Event Display run - where the user 
        // did not provide a relation input file
        log << MSG::WARNING << "No Relation Data Available" << endreq;
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
    sc = eventSvc()->registerObject(EventModel::Digi::TkrDigiHitTab,        m_tkrDigiToPosHitList );
    sc = eventSvc()->registerObject(EventModel::Digi::CalDigiHitTab,        m_calDigiToIntHitList );
    sc = eventSvc()->registerObject(EventModel::CalRecon::CalClusterHitTab, m_calClusterToXtalRecList );
    sc = eventSvc()->registerObject("/Event/MC/McPartToTrajectory",         m_mcPartToTrajList );
    sc = eventSvc()->registerObject("/Event/MC/McPointToPosHit",            m_mcPointToPosHitList );
    sc = eventSvc()->registerObject("/Event/MC/McPointToIntHit",            m_mcPointToIntHitList );
    
    return sc;
}

StatusCode relationRootReaderAlg::createTDSTables()
{
    // Purpose and Method:  Creates the TDS relational tables for output
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    m_tkrDigiToPosHitList     = new TkrDigiToPosHitTabList;
    m_calDigiToIntHitList     = new CalDigiToIntHitTabList;
    m_calClusterToXtalRecList = new CalClusterToXtalRecTabList;
    m_mcPartToTrajList        = new Event::RelationList<Event::McParticle,Event::McTrajectory>;
    m_mcPointToPosHitList     = new Event::RelationList<Event::McTrajectoryPoint,Event::McPositionHit>;
    m_mcPointToIntHitList     = new Event::RelationList<Event::McTrajectoryPoint,Event::McIntegratingHit>;

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

        const TObject* first = relation->getFirst();

        if (!first)
        {
            continue;
        }

        std::string className = first->ClassName();

        if (className != m_currentRelType)
        {
            m_currentRelType = className;

            std::map<std::string, ConvFuncPtr>::iterator converterIter = m_convFuncMap.find(className);

            if (converterIter != m_convFuncMap.end()) m_currentConverter = (*converterIter).second;
            else                                      m_currentConverter = 0;
        }

        if (m_currentConverter) sc = (this->*m_currentConverter)(relation);
        else
        {
            const TObject* tmp = relation->getSecond();
            log << MSG::WARNING << "Unrecognized Relation key: " << className << endreq;
        }
    }
    
    return sc;
}

StatusCode relationRootReaderAlg::readTkrDigiRelations(Relation* relation) 
{
    // Purpose and Method:  Converts a root TkrDigi/McPositionHit relation to TDS version
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    const TkrDigi*  digiRoot = dynamic_cast<const TkrDigi*>(relation->getFirst());
    Event::TkrDigi* digiTds  = 0;

    // Look up TDS to root relation for digis
    if (m_common.m_rootTkrDigiMap.find(digiRoot) != m_common.m_rootTkrDigiMap.end()) {
        digiTds = const_cast<Event::TkrDigi*>(m_common.m_rootTkrDigiMap[digiRoot]);
    } else {
        log << MSG::WARNING << "Could not located TkrDigi TDS/ROOT pair" << endreq;
        return sc;
    }

    const TObject* mcHitRoot = relation->getSecond();

    // Look up TDS to root relation for position hits
    if (m_common.m_rootMcPosHitMap.find(mcHitRoot) != m_common.m_rootMcPosHitMap.end()) 
    {
        Event::McPositionHit* mcHitTds = const_cast<Event::McPositionHit*>(m_common.m_rootMcPosHitMap[mcHitRoot]);

        // Make the event relation here
        typedef Event::Relation<Event::TkrDigi,Event::McPositionHit> relType;
        relType* rel = new relType(digiTds, mcHitTds);
        m_tkrDigiToPosHitList->push_back(rel);

        // Transfer the "infos" string(s)
        for(int idx = 0; idx < relation->getInfos().GetEntries(); idx++)
        {
            TObjString* objString  = (TObjString*)relation->getInfos().At(idx);
            std::string infoString = std::string(objString->GetString().Data());
            rel->addInfo(infoString);
        }
    } 
    else 
    {
        log << MSG::WARNING << "Could not located McPositionHit TDS/ROOT pair" << endreq;
    }

    return sc;
}


StatusCode relationRootReaderAlg::readCalDigiRelations(Relation* relation) 
{
    // Purpose and Method:  Converts a root CalDigi/McIntegratingHit relation to TDS version

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    const CalDigi*  digiRoot = dynamic_cast<const CalDigi*>(relation->getFirst());
    Event::CalDigi* digiTds  = 0;

    // Look up TDS to root relation for digis
    if (m_common.m_rootCalDigiMap.find(digiRoot) != m_common.m_rootCalDigiMap.end()) {
        digiTds = const_cast<Event::CalDigi*>(m_common.m_rootCalDigiMap[digiRoot]);
    } else {
        log << MSG::WARNING << "Could not located CalDigi TDS/ROOT pair" << endreq;
        return sc;
    }

    const TObject* mcHitRoot = relation->getSecond();

    // Look up TDS to root relation for position hits
    if (m_common.m_rootMcIntHitMap.find(mcHitRoot) != m_common.m_rootMcIntHitMap.end()) 
    {
        Event::McIntegratingHit* mcHitTds = const_cast<Event::McIntegratingHit*>(m_common.m_rootMcIntHitMap[mcHitRoot]);

        // Make the event relation here
        typedef Event::Relation<Event::CalDigi,Event::McIntegratingHit> relType;
        relType* rel = new relType(digiTds, mcHitTds);
        m_calDigiToIntHitList->push_back(rel);

        // Transfer the "infos" string(s)
        for(int idx = 0; idx < relation->getInfos().GetEntries(); idx++)
        {
            TObjString* objString  = (TObjString*)relation->getInfos().At(idx);
            std::string infoString = std::string(objString->GetString().Data());
            rel->addInfo(infoString);
        }
    } 
    else 
    {
        log << MSG::WARNING << "Could not located McIntegratingHit TDS/ROOT pair" << endreq;
    }

    return sc;
}

StatusCode relationRootReaderAlg::readCalXtalToClusterRelations(Relation* relation) 
{
    // Purpose and Method:  Converts a root CalDigi/McIntegratingHit relation to TDS version

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    const CalXtalRecData*  xtalRoot = dynamic_cast<const CalXtalRecData*>(relation->getFirst());
    Event::CalXtalRecData* xtalTds  = 0;

    // Look up TDS to root relation for digis
    if (m_common.m_rootCalXtalRecDataMap.find(xtalRoot) != m_common.m_rootCalXtalRecDataMap.end()) 
    {
        xtalTds = const_cast<Event::CalXtalRecData*>(m_common.m_rootCalXtalRecDataMap[xtalRoot]);
    } else {
        log << MSG::WARNING << "Could not located CalXtalRecData TDS/ROOT pair" << endreq;
        return sc;
    }

    const CalCluster* clusterRoot = dynamic_cast<const CalCluster*>(relation->getSecond());

    // Look up TDS to root relation for position hits
    if (m_common.m_rootCalClusterMap.find(clusterRoot) != m_common.m_rootCalClusterMap.end()) {
        Event::CalCluster* clusterTds = const_cast<Event::CalCluster*>(m_common.m_rootCalClusterMap[clusterRoot]);

        // Make the event relation here
        typedef Event::Relation<Event::CalXtalRecData,Event::CalCluster> relType;
        relType* rel = new relType(xtalTds,clusterTds);
        m_calClusterToXtalRecList->push_back(rel);

        // Transfer the "infos" string(s)
        for(int idx = 0; idx < relation->getInfos().GetEntries(); idx++)
        {
            TObjString* objString  = (TObjString*)relation->getInfos().At(idx);
            std::string infoString = std::string(objString->GetString().Data());
            rel->addInfo(infoString);
        }
    } 
    else 
    {
        log << MSG::WARNING << "Could not located CalCluster<-->CalXtalRecData TDS/ROOT pair" << endreq;
    }

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

    const TkrVertex*  vertexRoot   = dynamic_cast<const TkrVertex*>(relation->getFirst());

    const TkrTrack*   tkrTrackRoot = dynamic_cast<const TkrTrack*>(relation->getSecond());
    Event::TkrVertex* vertexTds    = 0;
    Event::TkrTrack*  tkrTrackTds  = 0;

    // Look up TDS to root relation for vertices
    if (m_common.m_rootTkrVertexMap.find(vertexRoot) != m_common.m_rootTkrVertexMap.end()) {
        vertexTds = const_cast<Event::TkrVertex*>(m_common.m_rootTkrVertexMap[vertexRoot]);
    } else {
        log << MSG::WARNING << "Could not located TkrVertex TDS/ROOT pair" << endreq;
    }

    // Look up TDS to root relation for candidate tracks
    if (m_common.m_rootTkrTrackMap.find(tkrTrackRoot) != m_common.m_rootTkrTrackMap.end()) 
    {
        tkrTrackTds = const_cast<Event::TkrTrack*>(m_common.m_rootTkrTrackMap[tkrTrackRoot]);

        // Make the event relation here
        //Event::TkrFitTrackRel* rel = new Event::TkrFitTrackRel(tkrCandTds, tkrTrackTds);
        //m_trackRelTab->addRelation(rel);
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

    const McParticle*  mcPartRoot = dynamic_cast<const McParticle*>(relation->getFirst());
    Event::McParticle* mcPartTds  = 0;

    // Look up TDS to root relation for McParticles
    if (m_common.m_rootMcPartMap.find(mcPartRoot) != m_common.m_rootMcPartMap.end()) {
        mcPartTds = const_cast<Event::McParticle*>(m_common.m_rootMcPartMap[mcPartRoot]);
    } else {
        log << MSG::WARNING << "Could not located McParticle TDS/ROOT pair" << endreq;
    }

    const McTrajectory* mcTrajRoot = dynamic_cast<const McTrajectory*>(relation->getSecond());

    // Look up TDS to root relation for position hits
    if (m_common.m_rootMcTrajectoryMap.find(mcTrajRoot) != m_common.m_rootMcTrajectoryMap.end()) 
    {
        Event::McTrajectory* mcTrajTds = const_cast<Event::McTrajectory*>(m_common.m_rootMcTrajectoryMap[mcTrajRoot]);

        // Make the event relation here
        Event::McPartToTrajectoryRel* rel = new Event::McPartToTrajectoryRel(mcPartTds, mcTrajTds);
        m_mcPartToTrajList->push_back(rel);

        // Transfer the "infos" string(s)
        for(int idx = 0; idx < relation->getInfos().GetEntries(); idx++)
        {
            TObjString* objString  = (TObjString*)relation->getInfos().At(idx);
            std::string infoString = std::string(objString->GetString().Data());
            rel->addInfo(infoString);
        }
    } 
    else 
    {
        log << MSG::WARNING << "Could not located McTrajectory TDS/ROOT pair" << endreq;
    }

    return sc;
}

StatusCode relationRootReaderAlg::readMcTrajectoryPointRelations(Relation* relation) 
{
    // Purpose and Method:  Converts a root TkrDigi/McPositionHit relation to TDS version
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    const McTrajectoryPoint*  mcPointRoot = dynamic_cast<const McTrajectoryPoint*>(relation->getFirst());
    Event::McTrajectoryPoint* mcPointTds  = 0;

    // Look up TDS to root relation for digis
    if (m_common.m_rootMcTrajectoryPointMap.find(mcPointRoot) != m_common.m_rootMcTrajectoryPointMap.end()) {
        mcPointTds = const_cast<Event::McTrajectoryPoint*>(m_common.m_rootMcTrajectoryPointMap[mcPointRoot]);
    } else {
        log << MSG::WARNING << "Could not located McTrajectoryPoint TDS/ROOT pair" << endreq;
    }

    const TObject* mcHitRoot = relation->getSecond();

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

}
 
void relationRootReaderAlg::endEvent() {
    if (m_relTab)  m_relTab->Clear(m_clearOption.c_str());
    m_relTab = 0;
}

StatusCode relationRootReaderAlg::finalize()
{
    close();
    
    StatusCode sc = StatusCode::SUCCESS;
    //setFinalized();
    return sc;
}


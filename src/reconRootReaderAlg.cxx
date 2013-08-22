#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/IIncidentListener.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "idents/CalXtalId.h"
#include "Event/Recon/AcdRecon/AcdRecon.h"
#include "Event/Recon/AcdRecon/AcdReconV2.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrTruncationInfo.h"  
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrTree.h"
#include "Event/Recon/TkrRecon/TkrFilterParams.h"
#include "Event/Recon/TkrRecon/TkrEventParams.h"
#include "Event/Recon/TkrRecon/TkrVecPointInfo.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/CalRecon/CalClusterMap.h"   
#include "Event/Recon/CalRecon/CalXtalRecData.h"  
#include "Event/Recon/CalRecon/CalEventEnergy.h"  
#include "Event/Recon/TreeClusterRelation.h"

#include "CLHEP/Matrix/Matrix.h"

#include "LdfEvent/EventSummaryData.h"

#include "AncillaryDataEvent/Recon.h"

#include "reconRootData/ReconEvent.h"

#include "facilities/Util.h"

#include "commonData.h"

#include "RootIo/IRootIoSvc.h"
#include "RootConvert/Utilities/RootReaderUtil.h"

// ADDED FOR THE FILE HEADERS DEMO
#include "RootIo/FhTool.h"

// low level converters
#include "RootConvert/Recon/TkrTruncationInfoConvert.h"
#include "RootConvert/Recon/CalClusterConvert.h"
#include "RootConvert/Recon/CalEventEnergyConvert.h"
#include "RootConvert/Recon/CalMipTrackConvert.h"
#include "RootConvert/Recon/CalXtalRecDataConvert.h"
#include "RootConvert/Recon/AcdReconConvert.h"
#include "RootConvert/Recon/AdfReconConvert.h"

#include <vector>
#include <map>
#include <stack>

/** @class reconRootReaderAlg
* @brief Reads Reconstruction data from a persistent ROOT file and stores the
* the data in the TDS.
*
* @author Heather Kelly
* $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/RootIo/src/reconRootReaderAlg.cxx,v 1.121 2013/06/19 05:10:19 lsrea Exp $
*/

class reconRootReaderAlg : public Algorithm, virtual public IIncidentListener
{   
public:

    reconRootReaderAlg(const std::string& name, ISvcLocator* pSvcLocator);

    /// Handles setup by opening ROOT file in read mode and creating a new TTree
    StatusCode initialize();

    /// Orchastrates reading from ROOT file and storing the data on the TDS for each event
    StatusCode execute();

    void handle(const Incident &inc) {
        if( inc.type()=="BeginEvent")beginEvent();
        else if(inc.type()=="EndEvent")endEvent();
    }

    void beginEvent() { };
    void endEvent();

    /// Closes the ROOT file and cleans up
    StatusCode finalize();

private:

    /// Reads top-level DigiEvent
    StatusCode readReconEvent();

    /// Reads TKR recon data from ROOT and puts data on TDS
    StatusCode readTkrRecon();

    StatusCode storeTkrTruncationInfo(TkrRecon *tkrRecRoot);
    StatusCode storeTkrClusterCol(TkrRecon *tkrRecRoot);

    StatusCode storeTrackAndVertexCol(TkrRecon *tkrRecRoot, bool vertexOnTdsFlag);
    StatusCode storeTkrVecPointCol(TkrRecon* tkrRecRoot);
    StatusCode rebuildTkrVecPointCol(Event::TkrVecPointCol* vecPointCol);
    StatusCode storeTkrVecPointsLinkCol(TkrRecon* tkrRecRoot);
    StatusCode storeTkrTreeCol(TkrRecon* tkrRecRoot);
    StatusCode storeTkrFilterParamsCol(TkrRecon* tkrRecoRoot);
    StatusCode storeTkrEventParams(TkrRecon* tkrRecRoot);
    StatusCode storeTkrVecPointInfo(TkrRecon* tkrRecRoot);

    /// convert a ROOT TkrFilterParams object to a TDS Event::TkrFilterParams object
    Event::TkrFilterParams* convertTkrFilterParams(const TkrFilterParams* filterParamsRoot);
    /// convert a ROOT TkrFitHit to a TDS Event::TkrFitHit
    Event::TkrTrackParams convertTkrTrackParams(const TkrTrackParams& paramsRoot);
    /// convert a ROOT TkrTrackHit to a TDS Event::TkrTrackHit
    Event::TkrTrackHit* convertTkrTrackHit(const TkrTrackHit* trackHitRoot);
    /// Stores TkrTrack objects into Event::TkrFitTrack
    Event::TkrTrack* convertTkrTrack(const TkrTrack* trackRoot);
    /// convert a ROOT TkrTree to a TDS Event::TkrTree
    Event::TkrTree* convertTkrTree(const TkrTree* treeRoot, Event::TkrVecNodeQueue* vecNodeQueue);
    /// convert at ROOT TkrTreeCompressed to a TDS Event::TkrTree
    Event::TkrTree* convertTkrTree(const TkrTreeCompressed* treeRoot, Event::TkrVecNodeQueue* vecNodeQueue);
    /// convert a ROOT TkrVecNode to a TDS Event::TkrVecNode
    Event::TkrVecNode* convertTkrVecNode(const TkrVecNode* node);
    /// convert a ROOT TkrVecNodeCompressed to a TDS Event::TkrVecNode
    Event::TkrVecNode* convertTkrVecNode(const TkrVecNodeCompressed* node);
    /// This to rebuild the sibling map after reading back in the TkrVecNode Tree
    int makeSiblingMap(Event::TkrVecNode* curNode, Event::TkrNodeSiblingMap* siblingMap);

    /// Reads CAL recon data from ROOT and puts data on TDS
    StatusCode readCalRecon();

    /// read in CAL xtal recon data from ROOT and store on the TDS
    StatusCode storeCalXtalRecDataCol(CalRecon *calRecRoot);

    /// read CAL cluster data from ROOT and store on TDS in CalClusterMap
    StatusCode storeCalClusterMap(CalRecon *calRecRoot);

    /// read CAL cluster data from ROOT and store on TDS
    StatusCode storeCalClusterCol(CalRecon *calRecRoot);

    /// read CAL cluster data from ROOT and store on TDS
    StatusCode storeCalMipTrackCol(CalRecon *calRecRoot);

    /// read CAL eventEnergy  data from ROOT and store on TDS
    StatusCode storeCalEventEnergyCol(CalRecon *calRecRoot);

    /// read CalEventEnergyMap data from ROOT and store on TDS
    StatusCode storeCalEventEnergyMap(CalRecon *calRecRoot);

    /// Reads ACD recon data from ROOT and puts data on the TDS
    StatusCode readAcdRecon();

    StatusCode readAdfRecon();

    /// The final stage is to handle any relations cross subsystem
    StatusCode readCalTkrAcdRelations();
    /// In particular the tree to cluster relations
    StatusCode storeTreeClusterRelations(TkrRecon* tkrRecRoot);

    /// Closes the ROOT file
    void close();

    /// Top-level Monte Carlo ROOT object
    ReconEvent *m_reconEvt;
    /// name of the input ROOT file
    std::string m_fileName;
    /// List of input files
    StringArrayProperty m_fileList;
    /// name of the Recon TTree stored in the ROOT file
    std::string m_treeName;
    std::string m_branchName;
    /// Option string which will be passed to McEvent::Clear
    std::string m_clearOption;
    /// Branch Exclusion List
    StringArrayProperty m_excludeBranchList;

    commonData m_common;

    IRootIoSvc*   m_rootIoSvc;

    bool FixAcdStreamerDone;

    bool m_terminateOnReadError;

    // ADDED FOR THE FILE HEADERS DEMO
    IFhTool * m_headersTool ;

    // The below is meant to provide mapping between clusters and TkrVecPoints in the event
    // we are reading in compressed TkrTree information
    typedef std::map<const Event::TkrCluster* , std::map<const Event::TkrCluster* , Event::TkrVecPoint*> > TkrVecPointLookUpMap;

    TkrVecPointLookUpMap m_tkrVecPointLookUpMap;
};

//static const AlgFactory<reconRootReaderAlg>  Factory;
//const IAlgFactory& reconRootReaderAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(reconRootReaderAlg);


reconRootReaderAlg::reconRootReaderAlg(const std::string& name, ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator), m_reconEvt(0)
{
    // Input pararmeters that may be set via the jobOptions file
    // Input ROOT file name, this will be overridden if RootIoSvc is provided a meta ROOT file for reading
    // Provided for backward-compatibility with older JO files, reconRootFileList is preferred JO parameter
    declareProperty("reconRootFile",m_fileName="");

    // Allows init of a TChain of files for reading, this will be overridden if RootIoSvc is provided
    // a meta ROOT file for reading
    StringArrayProperty initList;
    std::vector<std::string> initVec;
    initList.setValue(initVec);
    declareProperty("reconRootFileList", m_fileList=initList);
    initVec.clear();
    // Input TTree name
    declareProperty("reconTreeName", m_treeName="Recon");
    declareProperty("reconBranchName", m_branchName="ReconEvent");
    declareProperty("clearOption", m_clearOption="");
    declareProperty("ExcludeBranches", m_excludeBranchList=initList);

    m_tkrVecPointLookUpMap.clear();
}

StatusCode reconRootReaderAlg::initialize()
{
    // Purpose and Method:  Called once before the run begins.  This method
    //    opens a new ROOT file and prepares for reading.

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
        log << MSG::INFO << "Couldn't find the RootIoSvc!" << endreq;
        log << MSG::DEBUG << "Event loop will not terminate gracefully" << endreq;
        m_rootIoSvc = 0;
        // RootIoSvc is required for reading and writing, cannot continue without it
        return StatusCode::FAILURE;
    }

    if ( (m_fileList.value().size() > 0) && ( !m_fileName.empty() )) {
        log << MSG::WARNING << "Both reconRootFile and reconRootFileList have "
            << "been specified, reconRootFile is deprecated, please use "
            << "reconRootFileList" << endreq;
        return StatusCode::FAILURE;
    } else if ( (m_fileList.value().size() == 0) && ( !m_fileName.empty() ) )
        m_rootIoSvc->appendFileList(m_fileList, m_fileName);
    else if (m_fileList.value().size() == 0)
        m_rootIoSvc->appendFileList(m_fileList, "recon.root");

    m_terminateOnReadError = m_rootIoSvc->terminateOnReadError();


    // Set up new school system...
    // Use treeName as key type
    m_rootIoSvc->prepareRootInput("recon", m_treeName, m_branchName, 0, m_fileList);

    if (m_excludeBranchList.value().size() > 0) {
        std::vector<std::string>::const_iterator excludeListItr;
        for (excludeListItr = m_excludeBranchList.value().begin();
            excludeListItr != m_excludeBranchList.value().end();
            excludeListItr++ ) {
                std::string branchName = *excludeListItr;
                bool foundFlag  = m_rootIoSvc->setBranchStatus("recon",branchName,0);
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

StatusCode reconRootReaderAlg::execute()
{
    // Purpose and Method:  Called once per event.  This method calls
    //   the appropriate methods to read data from the ROOT file and store
    //   data on the TDS.

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Try reading the event this way... 
    // Use treeName as key type
    m_reconEvt = dynamic_cast<ReconEvent*>(m_rootIoSvc->getNextEvent("recon")) ;

    if (!m_reconEvt) {
        if (m_terminateOnReadError) {
            log << MSG::ERROR << "Failed to read in Recon data" << endreq;
            return StatusCode::FAILURE;
        }
        // Do not fail if there was no RECON data to read - this may be an Event Display run - where the user 
        // did not provide an RECON input file
        log << MSG::WARNING << "No Recon data Available" << endreq;
        return StatusCode::SUCCESS;
    }

    sc = readReconEvent();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to read top level ReconEvent" << endreq;
        return sc;
    }

    sc = readTkrRecon();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to load Tkr Recon" << endreq;
        return sc;
    }

    sc = readCalRecon();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to load Cal Recon" << endreq;
        return sc;
    }

    sc = readAcdRecon();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to load Acd Recon" << endreq;
        return sc;
    }

    sc = readAdfRecon();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to load Adf Recon" << endreq;
        return sc;
    }
    sc = readCalTkrAcdRelations();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to load Acd, Cal, Tkr relations" << endreq;
        return sc;
    }

    //    evtId = readInd+1;
    //    saveDir->cd();

    return sc;
}


StatusCode reconRootReaderAlg::readReconEvent() {

    MsgStream log(msgSvc(), name());

    StatusCode sc = StatusCode::SUCCESS;

    // Retrieve the Event data for this event
    SmartDataPtr<Event::EventHeader> evt(eventSvc(), EventModel::EventHeader);
    if (!evt) {
        log << MSG::ERROR << "Failed to retrieve Event" << endreq;
        return StatusCode::FAILURE;
    }

    unsigned int eventIdTds = evt->event();
    unsigned int runIdTds = evt->run();

    unsigned int eventIdRoot = m_reconEvt->getEventId();
    unsigned int runIdRoot = m_reconEvt->getRunId();

    // Check to see if the event and run ids have already been set.
    if (eventIdTds != eventIdRoot) evt->setEvent(eventIdRoot);
    if (runIdTds != runIdRoot) evt->setRun(runIdRoot);

    log << MSG::DEBUG << "Reading Event (run, event): (" << runIdRoot
        << ", " << eventIdRoot << ")" << endreq;

    evt->setGleamEventFlags(m_reconEvt->getGleamEventFlags());

    // Only update the eventflags on the TDS if the /Event/EventSummary
    // does not yet exist (digiRootReader may fill this for us)
    SmartDataPtr<LdfEvent::EventSummaryData> summaryTds(eventSvc(), "/Event/EventSummary");
    if (!summaryTds) {
        LdfEvent::EventSummaryData *evtSumTds = new LdfEvent::EventSummaryData();
        evtSumTds->initEventFlags(m_reconEvt->getEventFlags());
        sc = eventSvc()->registerObject("/Event/EventSummary", evtSumTds);
        if( sc.isFailure() ) {
            log << MSG::ERROR << "could not register /Event/EventSummary " << endreq
                ;
            return sc;
        }
    }

    return sc;
}

StatusCode reconRootReaderAlg::readTkrRecon() {
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Let's just always make sure this gets cleared
    m_tkrVecPointLookUpMap.clear();

    TkrRecon *tkrRecRoot = m_reconEvt->getTkrRecon();
    if(!tkrRecRoot) return sc;

    // Make sure the TkrRecon branch exists on the TDS
    DataObject* pnode =0;
    sc = eventSvc()->retrieveObject(EventModel::TkrRecon::Event, pnode);
    if( sc.isFailure() ) 
    {
        sc = eventSvc()->registerObject(EventModel::TkrRecon::Event,new DataObject);
        if( sc.isFailure() ) 
        {
            log << MSG::ERROR << "Could not create TkrRecon directory" << endreq;
            return sc;
        }
    }

    // check to see if TKR cluster collection already exists on TDS
    // If not, store cluster collection on the TDS.
    SmartDataPtr<Event::TkrClusterCol> clusterColTds(eventSvc(), EventModel::TkrRecon::TkrClusterCol);
    if (clusterColTds) {
        log << MSG::INFO << "Tkr Cluster Collection is already on the TDS" 
            << endreq;
    } else {
        sc = storeTkrClusterCol(tkrRecRoot);
        if (sc.isFailure()) {
            log << MSG::ERROR << "failed to store TKR cluster collection on TDS" << endreq;
            return sc;
        }
    }

    // check to see if TkrVecPoint collection already exists on TDS
    // If not, store cluster collection on the TDS.
    SmartDataPtr<Event::TkrVecPointCol> tkrVecPointColTds(eventSvc(), EventModel::TkrRecon::TkrVecPointCol);
    if (tkrVecPointColTds) {
        log << MSG::INFO << "TkrVecPoint Collection is already on the TDS" 
            << endreq;
    } else {
        sc = storeTkrVecPointCol(tkrRecRoot);
        if (sc.isFailure()) {
            log << MSG::ERROR << "failed to store TkrVecPoint collection on TDS" << endreq;
            return sc;
        }
    }

    // check to see if TkrVecPointsLink collection already exists on TDS
    // If not, store cluster collection on the TDS.
    SmartDataPtr<Event::TkrVecPointsLinkCol> tkrVecPointsLinkColTds(eventSvc(), EventModel::TkrRecon::TkrVecPointsLinkCol);
    if (tkrVecPointsLinkColTds) {
        log << MSG::INFO << "TkrVecPointsLink Collection is already on the TDS" 
            << endreq;
    } else {
        sc = storeTkrVecPointsLinkCol(tkrRecRoot);
        if (sc.isFailure()) {
            log << MSG::ERROR << "failed to store TkrVecPointsLink collection on TDS" << endreq;
            return sc;
        }
    }

    // check to see if TkrEventParams object already exists on TDS
    // If not, store a new one on the TDS.
    SmartDataPtr<Event::TkrEventParams> tkrEventParamsTds(eventSvc(), EventModel::TkrRecon::TkrEventParams);
    if (tkrEventParamsTds) {
        log << MSG::INFO << "TkrEventParamsis already on the TDS" 
            << endreq;
    } else {
        sc = storeTkrEventParams(tkrRecRoot);
        if (sc.isFailure()) {
            log << MSG::ERROR << "failed to store TkrEventParams on TDS" << endreq;
            return sc;
        }
    }

    // check to see if TkrFilterParams collection already exists on TDS
    // If not, store filter params collection on the TDS.
    SmartDataPtr<Event::TkrFilterParamsCol> tkrFilterParamsColTds(eventSvc(), EventModel::TkrRecon::TkrFilterParamsCol);
    if (tkrFilterParamsColTds) {
        log << MSG::INFO << "TkrFilterParams Collection is already on the TDS" 
            << endreq;
    } else {
        sc = storeTkrFilterParamsCol(tkrRecRoot);
        if (sc.isFailure()) {
            log << MSG::ERROR << "failed to store TkrFilterParams collection on TDS" << endreq;
            return sc;
        }
    }

    // check to see if TKR vertex collection exists on TDS already
    // Set a boolean flag if the vertex collection already exists
    bool vertexOnTdsFlag = false;
    SmartDataPtr<Event::TkrVertexCol> vertexColTds(eventSvc(), EventModel::TkrRecon::TkrVertexCol);
    if (vertexColTds) {
        log << MSG::INFO << "Tkr VertexCol is already on TDS" << endreq;
        vertexOnTdsFlag = true;
    }

    // check to see if TKR track collection exists on the TDS already
    SmartDataPtr<Event::TkrTrackCol> trackColTds(eventSvc(), EventModel::TkrRecon::TkrTrackCol);
    if (trackColTds) {
        log << MSG::INFO << "TkrTrackCol is already on TDS" << endreq;
    } else {
        sc = storeTrackAndVertexCol(tkrRecRoot, vertexOnTdsFlag);
        if (sc.isFailure()) {
            log << MSG::ERROR << "failed to store FitTrackCol on the TDS" << endreq;
            return sc;
        }
    }

    // check to see if TkrTree collection already exists on TDS
    // If not, store tree collection on the TDS.
    SmartDataPtr<Event::TkrTreeCol> tkrTreeColTds(eventSvc(), EventModel::TkrRecon::TkrTreeCol);
    if (tkrTreeColTds) {
        log << MSG::INFO << "TkrTree Collection is already on the TDS" 
            << endreq;
    } else {
        sc = storeTkrTreeCol(tkrRecRoot);
        if (sc.isFailure()) {
            log << MSG::ERROR << "failed to store TkrTree collection on TDS" << endreq;
            return sc;
        }
    }

    //check to see if TKR TruncationInfo exists on the TDS already
    SmartDataPtr<Event::TkrTruncationInfo> truncationInfoTds(eventSvc(),EventModel::TkrRecon::TkrTruncationInfo);
    if(truncationInfoTds) {
        log << MSG::INFO << "TkrTruncationInfo is already on TDS" << endreq;
    } else {
        sc = storeTkrTruncationInfo(tkrRecRoot);
        if (sc.isFailure()) {
            log << MSG::ERROR << "failed to store TkrtruncationInfo on the TDS" << endreq;
            return sc;
        }
    }

    //check to see if TkrVecPointInfo exists on the TDS already
    SmartDataPtr<Event::TkrVecPointInfo> tkrVecPointInfoTds(eventSvc(),EventModel::TkrRecon::TkrVecPointInfo);
    if(tkrVecPointInfoTds) {
        log << MSG::INFO << "TkrVecPointInfo is already on TDS" << endreq;
    } else {
        sc = storeTkrVecPointInfo(tkrRecRoot);
        if (sc.isFailure()) {
            log << MSG::ERROR << "failed to store TkrVecPointInfo on the TDS" << endreq;
            return sc;
        }
    }

    //    // check to see if TKR CRtrack collection exists on the TDS already
    //SmartDataPtr<Event::TkrTrackCol> crTrackColTds(eventSvc(), EventModel::TkrRecon::TkrCRTrackCol);
    //if (trackColTds) {
    //    log << MSG::INFO << "TkrCRTrackCol is already on TDS" << endreq;
    //} else {
    //    sc = storeTrackAndVertexCol(tkrRecRoot, vertexOnTdsFlag);
    //    if (sc.isFailure()) {
    //        log << MSG::ERROR << "failed to store FitTrackCol on the TDS" << endreq;
    //        return sc;
    //    }
    //}

    return sc;
}

StatusCode reconRootReaderAlg::storeTkrClusterCol(TkrRecon *tkrRecRoot) {
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    Event::TkrClusterCol *clusterTdsCol = new Event::TkrClusterCol;

    const TObjArray * clusterRootCol = tkrRecRoot->getClusterCol();
    TIter clusterIter(clusterRootCol);
    TkrCluster *clusterRoot = 0;

    Event::TkrIdClusterMap* clusMap = new Event::TkrIdClusterMap;
    sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrIdClusterMap,
        clusMap);
    if (sc.isFailure()) {
        log << MSG::ERROR << "failed to register TkrIdClusterMap on the TDS" << endreq;
        return sc;
    }

    while ((clusterRoot = (TkrCluster*)clusterIter.Next())!=0) 
    {
        commonRootData::TkrId tkrIdRoot = clusterRoot->getTkrId();

        idents::TkrId tkrId(tkrIdRoot.getTowerX(),
            tkrIdRoot.getTowerY(),
            tkrIdRoot.getTray(),
            tkrIdRoot.getBotTop() == commonRootData::TkrId::eTKRSiTop,
            tkrIdRoot.getView() );

        TVector3 posRoot = clusterRoot->getPosition();
        Point posTds(posRoot.X(), posRoot.Y(), posRoot.Z());

        UInt_t status = clusterRoot->getStatusWord(); // need this later

        Event::TkrCluster* clusterTds 
            = new Event::TkrCluster(tkrId, 
            clusterRoot->getFirstStrip(),
            clusterRoot->getLastStrip(),
            posTds,
            (int)clusterRoot->getRawToT(),
            clusterRoot->getMips(),
            status,
            clusterRoot->getNBad());

        clusterTdsCol->push_back(clusterTds);
        // 
        if((status&Event::TkrCluster::maskMERGERESULT)==0)(*clusMap)[tkrId].push_back(clusterTds);

        m_common.m_rootTkrClusterMap[clusterRoot] = clusterTds;
    }

    sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrClusterCol, clusterTdsCol);
    if (sc.isFailure()) {

        log << MSG::DEBUG;
        if( log.isActive()) log.stream() << "Failed to register TkrClusterCol";
        log << endreq;
        return StatusCode::FAILURE;
    }


    return sc;
}

StatusCode reconRootReaderAlg::storeTkrVecPointCol(TkrRecon *tkrRecRoot) 
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Create a TDS collection
    Event::TkrVecPointCol* tkrVecPointTdsCol = new Event::TkrVecPointCol;

    if (tkrRecRoot->nTkrVecPoints() > 0)
    {
        // Does the input data contain a vec point collection
        const TObjArray* tkrVecPointRootCol = tkrRecRoot->getTkrVecPointCol();

        TIter tkrVecPointIter(tkrVecPointRootCol);
        TkrVecPoint* tkrVecPointRoot = 0;

        while ((tkrVecPointRoot = (TkrVecPoint*)tkrVecPointIter.Next())!=0) 
        {
            const Event::TkrCluster* xCluster = m_common.m_rootTkrClusterMap[tkrVecPointRoot->getXCluster()];
            const Event::TkrCluster* yCluster = m_common.m_rootTkrClusterMap[tkrVecPointRoot->getYCluster()];

            Event::TkrVecPoint* tkrVecPointTds = new Event::TkrVecPoint();

            tkrVecPointTds->initialize(tkrVecPointRoot->getLayer(),
                tkrVecPointRoot->getStatusWord(),
                xCluster,
                yCluster);

            tkrVecPointTdsCol->push_back(tkrVecPointTds);

            m_common.m_rootTkrVecPointMap[tkrVecPointRoot] = tkrVecPointTds;
        }
    }
    // One other thing to look for, if we were storing in "compressed Tree" mode then 
    // we want to build the TkrVecPoints from the clusters
    else if (tkrRecRoot->nTkrTrees() > 0 && tkrRecRoot->getTkrTreeCompressed(0) != 0)
    {
        rebuildTkrVecPointCol(tkrVecPointTdsCol);
    }

    // Regardless of what happened, store the object on the TDS
    sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrVecPointCol, tkrVecPointTdsCol);
    if (sc.isFailure()) {

        log << MSG::DEBUG;
        if( log.isActive()) log.stream() << "Failed to register TkrVecPointCol";
        log << endreq;
        return StatusCode::FAILURE;
    }

    return sc;
}

StatusCode reconRootReaderAlg::rebuildTkrVecPointCol(Event::TkrVecPointCol* vecPointCol)
{
    StatusCode sc = StatusCode::SUCCESS;

    // First we need to retrieve the list of clusters in order to be sure that we have something to do here
    SmartDataPtr<Event::TkrClusterCol> clusterColTds(eventSvc(), EventModel::TkrRecon::TkrClusterCol);
    if (clusterColTds)
    {
        // To reconstitute we are going to double loop through cluster list and not care, for now, about order
        for(Event::TkrClusterCol::iterator outsideItr  = clusterColTds->begin();
            outsideItr != clusterColTds->end();
            outsideItr++)
        {
            Event::TkrCluster* outsideCluster = *outsideItr;

            // if this cluster was merged into a super cluster then skip it here
            if (outsideCluster->isSet(Event::TkrCluster::maskMERGED)) continue;

            for(Event::TkrClusterCol::iterator insideItr  = outsideItr + 1;
                insideItr != clusterColTds->end();
                insideItr++)
            {
                Event::TkrCluster* insideCluster = *insideItr;

                // Similarly, if this cluster was merged into a super cluster then skip it here
                if (insideCluster->isSet(Event::TkrCluster::maskMERGED)) continue;

                // Clusters must be same bilayer and same tower to make TkrVecPoint
                if (outsideCluster->getLayer() == insideCluster->getLayer() &&
                    outsideCluster->tower()    == insideCluster->tower())
                {
                    // Swap clusters to get an X and Y version
                    Event::TkrCluster* clusterX = outsideCluster;
                    Event::TkrCluster* clusterY = insideCluster;

                    if (insideCluster->getTkrId().getView() == idents::TkrId::eMeasureX)
                    {
                        clusterX = insideCluster;
                        clusterY = outsideCluster;
                    }

                    Event::TkrVecPoint* vecPoint =  new Event::TkrVecPoint(clusterX->getLayer(), clusterX, clusterY);

                    vecPointCol->push_back(vecPoint);

                    m_tkrVecPointLookUpMap[clusterX][clusterY] = vecPoint;
                }
            }
        }
    }

    return sc;
}

StatusCode reconRootReaderAlg::storeTkrVecPointInfo(TkrRecon* tkrRecRoot)
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Create a new TkrVecPointInfo object
    Event::TkrVecPointInfo* vecPointInfoTds = new Event::TkrVecPointInfo();

    // Fill what we can
    vecPointInfoTds->setMaxNumSkippedLayers(tkrRecRoot->getTkrVecPointInfo().getMaxNumLinkCombinations());
    vecPointInfoTds->setNumTkrVecPoints(tkrRecRoot->getTkrVecPointInfo().getNumTkrVecPoints());
    vecPointInfoTds->setNumBiLayersWVecPoints(tkrRecRoot->getTkrVecPointInfo().getNumBiLayersWVecPoints());
    vecPointInfoTds->setMaxNumLinkCombinations(tkrRecRoot->getTkrVecPointInfo().getMaxNumLinkCombinations());

    log << MSG::DEBUG << "Vec vars " << vecPointInfoTds->getMaxNumSkippedLayers() << " " 
        << vecPointInfoTds->getNumTkrVecPoints() << endreq;

    // Store on the TDS
    sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrVecPointInfo, vecPointInfoTds);
    if (sc.isFailure()) {

        log << MSG::DEBUG;
        if( log.isActive()) log.stream() << "Failed to register TkrVecPointInfo in TDS";
        log << endreq;
        return StatusCode::FAILURE;
    }

    return sc;
}

StatusCode reconRootReaderAlg::storeTkrVecPointsLinkCol(TkrRecon *tkrRecRoot) 
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    Event::TkrVecPointsLinkCol* tkrVecPointsLinkTdsCol = new Event::TkrVecPointsLinkCol;

    const TObjArray * tkrVecPointsLinkRootCol = tkrRecRoot->getTkrVecPointsLinkCol();
    TIter tkrVecPointsLinkIter(tkrVecPointsLinkRootCol);
    TkrVecPointsLink* tkrVecPointsLinkRoot = 0;

    while ((tkrVecPointsLinkRoot = (TkrVecPointsLink*)tkrVecPointsLinkIter.Next())!=0) 
    {
        const Event::TkrVecPoint* topPoint = m_common.m_rootTkrVecPointMap[tkrVecPointsLinkRoot->getFirstVecPoint()];
        const Event::TkrVecPoint* botPoint = m_common.m_rootTkrVecPointMap[tkrVecPointsLinkRoot->getSecondVecPoint()];

        Event::TkrVecPointsLink* tkrVecPointsLinkTds = new Event::TkrVecPointsLink(topPoint, 
            botPoint, 
            tkrVecPointsLinkRoot->getMaxScatAngle());

        // Fill in the rest of the parameters
        tkrVecPointsLinkTds->updateStatusBits(tkrVecPointsLinkRoot->getStatusBits());
        tkrVecPointsLinkTds->setMaxScatAngle(tkrVecPointsLinkRoot->getMaxScatAngle());

        tkrVecPointsLinkTdsCol->push_back(tkrVecPointsLinkTds);

        m_common.m_rootTkrVecPointsLinkMap[tkrVecPointsLinkRoot] = tkrVecPointsLinkTds;
    }

    sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrVecPointsLinkCol, tkrVecPointsLinkTdsCol);
    if (sc.isFailure()) {

        log << MSG::DEBUG;
        if( log.isActive()) log.stream() << "Failed to register TkrVecPointCol";
        log << endreq;
        return StatusCode::FAILURE;
    }

    return sc;
}

StatusCode reconRootReaderAlg::storeTkrTreeCol(TkrRecon* tkrRecRoot)
{
    // Purpose and Method:  Retrieve track and vertex collections from ROOT file and
    //  Store the transient versions on the TDS
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // If no Trees then no work
    if (tkrRecRoot->nTkrTrees())
    {
        // Create a Tree collection for the TDS
        Event::TkrTreeCol* treeColTds = new Event::TkrTreeCol();

        // We also want a container for the TkrVecNodes which make up the Tree...
        // Note that putting them in an object vector is good for making sure they
        // are cleared before the next event but is NOT the same container they were 
        // stored in during reconstruction. However, at that time we needed a container
        // which would automagically sort them, now we don't need that.
        Event::TkrVecNodeQueue* vecNodeQueue = new Event::TkrVecNodeQueue();

        // Loop through the Trees in the input root object
        for(int idx = 0; idx < tkrRecRoot->nTkrTrees(); idx++)
        {
            Event::TkrTree* treeTds = 0;

            // Recover pointer to root version of Tree
            TkrTree* treeRoot = tkrRecRoot->getTkrTree(idx);

            if (treeRoot != 0) 
            {
                treeTds = convertTkrTree(treeRoot, vecNodeQueue);

                m_common.m_rootTkrTreeMap[treeRoot] = treeTds;
            }
            else
            {
                TkrTreeCompressed* treeRoot = tkrRecRoot->getTkrTreeCompressed(idx);

                treeTds = convertTkrTree(treeRoot, vecNodeQueue);

                m_common.m_rootTkrTreeMap[treeRoot] = treeTds;
            }

            // Finally! Store in the TDS object
            treeColTds->push_back(treeTds);
        }

        // Store in the Tree collection TDS
        sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrTreeCol, treeColTds);
        if (sc.isFailure()) 
        {
            log << MSG::DEBUG;
            if( log.isActive()) log.stream() << "Failed to register TkrTreeCol";
            log << endreq;
            return StatusCode::FAILURE;
        }

        // Store in the TkrVecNode collection TDS
        sc = eventSvc()->registerObject("/Event/TkrRecon/TkrVecNodeCol", vecNodeQueue);
        if (sc.isFailure()) 
        {
            log << MSG::DEBUG;
            if( log.isActive()) log.stream() << "Failed to register TkrVecQueue";
            log << endreq;
            return StatusCode::FAILURE;
        }
    }

    return sc;
}

Event::TkrTree* reconRootReaderAlg::convertTkrTree(const TkrTree* treeRoot, Event::TkrVecNodeQueue* vecNodeQueue)
{
    // Recover pointer to the root version of the head node in the Tree
    const TkrVecNode* headNodeRoot = treeRoot->getHeadNode();

    // The head node is the full Tree, go through and convert to TDS equivalent
    Event::TkrVecNode* headNodeTds = convertTkrVecNode(headNodeRoot);

    // Keep track of the head node in our Queue structure
    vecNodeQueue->push(headNodeTds);

    // Now look up the pointers to the best and second best branches
    const Event::TkrVecNode* bestLeafTds = m_common.m_rootTkrVecNodeMap[treeRoot->getBestLeaf()];
    const Event::TkrVecNode* nextLeafTds = m_common.m_rootTkrVecNodeMap[treeRoot->getSecondLeaf()];

    // Do we really need this?
    Event::TkrNodeSiblingMap* nodeSiblingMapTds = new Event::TkrNodeSiblingMap();

    makeSiblingMap(headNodeTds, nodeSiblingMapTds);

    // Recover the tracks associated to this Tree
    const Event::TkrTrack* bestTrackTds = m_common.m_rootTkrTrackMap[treeRoot->getBestTrack()];
    const Event::TkrTrack* nextTrackTds = 0;

    if (treeRoot->GetEntries() > 1) nextTrackTds = m_common.m_rootTkrTrackMap[treeRoot->At(1)];

    // Recover the Tree axis params
    Event::TkrFilterParams* treeParamsTds = convertTkrFilterParams(treeRoot->getAxisParams());

    // Create a shiny new TkrTree to fill with information and store in the TDS
    Event::TkrTree* treeTds = new Event::TkrTree(headNodeTds, 
        const_cast<Event::TkrVecNode*>(bestLeafTds),
        const_cast<Event::TkrVecNode*>(nextLeafTds),
        nodeSiblingMapTds, 
        treeParamsTds, 
        const_cast<Event::TkrTrack*>(bestTrackTds) );

    // Set the two angles:
    treeTds->setBestBranchAngleToAxis(treeRoot->getBestBranchAngleToAxis());
    treeTds->setAxisSeededAngleToAxis(treeRoot->getAxisSeededAngleToAxis());

    // Add the second track if it exists
    if (nextTrackTds) treeTds->push_back(const_cast<Event::TkrTrack*>(nextTrackTds));

    return treeTds;
}

Event::TkrTree* reconRootReaderAlg::convertTkrTree(const TkrTreeCompressed* treeRoot, Event::TkrVecNodeQueue* vecNodeQueue)
{
    // Recover pointer to the root version of the head node in the Tree
    const TkrVecNodeCompressed* headNodeRoot = treeRoot->getHeadNode();

    // The head node is the full Tree, go through and convert to TDS equivalent
    Event::TkrVecNode* headNodeTds = convertTkrVecNode(headNodeRoot);

    // Keep track of the head node in our Queue structure
    vecNodeQueue->push(headNodeTds);

    // Now look up the pointers to the best and second best branches
    const Event::TkrVecNode* bestLeafTds = m_common.m_rootTkrVecNodeMap[treeRoot->getBestLeaf()];
    const Event::TkrVecNode* nextLeafTds = m_common.m_rootTkrVecNodeMap[treeRoot->getSecondLeaf()];

    // Do we really need this?
    Event::TkrNodeSiblingMap* nodeSiblingMapTds = new Event::TkrNodeSiblingMap();

    makeSiblingMap(headNodeTds, nodeSiblingMapTds);

    // Recover the tracks associated to this Tree
    const Event::TkrTrack* bestTrackTds = m_common.m_rootTkrTrackMap[treeRoot->getBestTrack()];
    const Event::TkrTrack* nextTrackTds = 0;

    if (treeRoot->GetEntries() > 1) nextTrackTds = m_common.m_rootTkrTrackMap[treeRoot->At(1)];

    // Recover the Tree axis params
    Event::TkrFilterParams* treeParamsTds = convertTkrFilterParams(treeRoot->getAxisParams());

    // Create a shiny new TkrTree to fill with information and store in the TDS
    Event::TkrTree* treeTds = new Event::TkrTree(headNodeTds, 
        const_cast<Event::TkrVecNode*>(bestLeafTds),
        const_cast<Event::TkrVecNode*>(nextLeafTds),
        nodeSiblingMapTds, 
        treeParamsTds, 
        const_cast<Event::TkrTrack*>(bestTrackTds) );

    // Set the two angles:
    treeTds->setBestBranchAngleToAxis(treeRoot->getBestBranchAngleToAxis());
    treeTds->setAxisSeededAngleToAxis(treeRoot->getAxisSeededAngleToAxis());

    // Add the second track if it exists
    if (nextTrackTds) treeTds->push_back(const_cast<Event::TkrTrack*>(nextTrackTds));

    return treeTds;
}

Event::TkrVecNode* reconRootReaderAlg::convertTkrVecNode(const TkrVecNode* tkrVecNodeRoot)
{
    Event::TkrVecNode* tkrVecNodeTds = 0;

    if (tkrVecNodeRoot)
    {
        // Start by getting a new TDS version of the TkrVecNode object
        tkrVecNodeTds = new Event::TkrVecNode();

        // We'll need pointers to the associated link and parent node
        const Event::TkrVecPointsLink* associatedLinkTds = 0;
        const Event::TkrVecNode*       parentNodeTds     = 0;

        if (tkrVecNodeRoot->getAssociatedLink()) associatedLinkTds = m_common.m_rootTkrVecPointsLinkMap[tkrVecNodeRoot->getAssociatedLink()];
        if (tkrVecNodeRoot->getParentNode())     parentNodeTds     = m_common.m_rootTkrVecNodeMap[tkrVecNodeRoot->getParentNode()];

        // Initialize the new TDS TkrVecNode object
        tkrVecNodeTds->initializeInfo(const_cast<Event::TkrVecPointsLink*>(associatedLinkTds),
            const_cast<Event::TkrVecNode*>(parentNodeTds),
            tkrVecNodeRoot->getStatusBits(),
            tkrVecNodeRoot->getRmsAngleSum(),
            tkrVecNodeRoot->getNumAnglesInSum(),
            tkrVecNodeRoot->getNumLeaves(),
            tkrVecNodeRoot->getNumBranches(),
            tkrVecNodeRoot->getDepth(),
            tkrVecNodeRoot->getBestNumBiLayers(),
            tkrVecNodeRoot->getBestRmsAngle() );

        // Store this in the TDS collection

        // Keep track of root/TDS relationship
        m_common.m_rootTkrVecNodeMap[tkrVecNodeRoot] = tkrVecNodeTds;

        // Now loop through the daughters of this node and continue to fill information
        for (int idx = 0; idx < tkrVecNodeRoot->GetEntries(); idx++)
        {
            // Recover the pointer to the root version
            TkrVecNode* daughterNodeRoot = (TkrVecNode*)tkrVecNodeRoot->At(idx);

            // Convert to TDS version
            Event::TkrVecNode* daughterNodeTds = convertTkrVecNode(daughterNodeRoot);

            // Add daughter to our current node
            if (daughterNodeTds) tkrVecNodeTds->push_back(daughterNodeTds);
        }
    }

    return tkrVecNodeTds;
}

Event::TkrVecNode* reconRootReaderAlg::convertTkrVecNode(const TkrVecNodeCompressed* tkrVecNodeRoot)
{
    Event::TkrVecNode* tkrVecNodeTds = 0;

    if (tkrVecNodeRoot)
    {
        // Start by getting a new TDS version of the TkrVecNode object
        tkrVecNodeTds = new Event::TkrVecNode();

        // We'll need pointers to the associated link and parent node
        Event::TkrVecPointsLink* associatedLinkTds = 0;
        const Event::TkrVecNode* parentNodeTds     = 0;

        // Recover pointer to parent (if one)
        if (tkrVecNodeRoot->getParentNode()) parentNodeTds = m_common.m_rootTkrVecNodeMap[tkrVecNodeRoot->getParentNode()];

        // Here we deal with recreating the associated link if there is one.(signaled by a non-zero pointer to a cluster
        if (tkrVecNodeRoot->getTopPointClusterX() != 0)
        {
            // Start by recovering the top point
            const Event::TkrCluster*  topClusterXTds = m_common.m_rootTkrClusterMap[tkrVecNodeRoot->getTopPointClusterX()];
            const Event::TkrCluster*  topClusterYTds = m_common.m_rootTkrClusterMap[tkrVecNodeRoot->getTopPointClusterY()];
            Event::TkrVecPoint*       topPoint       = m_tkrVecPointLookUpMap[topClusterXTds][topClusterYTds];

            // reset the top points status
            //topPoint->set

            // Now recover the bottom point
            const Event::TkrCluster*  botClusterXTds = m_common.m_rootTkrClusterMap[tkrVecNodeRoot->getBotPointClusterX()];
            const Event::TkrCluster*  botClusterYTds = m_common.m_rootTkrClusterMap[tkrVecNodeRoot->getBotPointClusterY()];
            Event::TkrVecPoint*       botPoint       = m_tkrVecPointLookUpMap[botClusterXTds][botClusterYTds];

            // Create a new TkrVecPointsLink and initialize it
            associatedLinkTds = new Event::TkrVecPointsLink(topPoint, botPoint, 0.1);

            // Reset the status bits
            associatedLinkTds->updateStatusBits(tkrVecNodeRoot->getLinkStatus());

            // Store link in the TDS collection
            SmartDataPtr<Event::TkrVecPointsLinkCol> vecPointsLinkCol(eventSvc(), EventModel::TkrRecon::TkrVecPointsLinkCol);
            if (vecPointsLinkCol) vecPointsLinkCol->push_back(associatedLinkTds);
        }

        // Initialize the new TDS TkrVecNode object
        tkrVecNodeTds->initializeInfo(const_cast<Event::TkrVecPointsLink*>(associatedLinkTds),
            const_cast<Event::TkrVecNode*>(parentNodeTds),
            tkrVecNodeRoot->getStatusBits(),
            tkrVecNodeRoot->getRmsAngleSum(),
            tkrVecNodeRoot->getNumAnglesInSum(),
            tkrVecNodeRoot->getNumLeaves(),
            tkrVecNodeRoot->getNumBranches(),
            tkrVecNodeRoot->getDepth(),
            tkrVecNodeRoot->getBestNumBiLayers(),
            tkrVecNodeRoot->getBestRmsAngle() );

        // Keep track of root/TDS relationship
        m_common.m_rootTkrVecNodeMap[tkrVecNodeRoot] = tkrVecNodeTds;

        // Now loop through the daughters of this node and continue to fill information
        for (int idx = 0; idx < tkrVecNodeRoot->GetEntries(); idx++)
        {
            // Recover the pointer to the root version
            TkrVecNodeCompressed* daughterNodeRoot = (TkrVecNodeCompressed*)tkrVecNodeRoot->At(idx);

            // Convert to TDS version
            Event::TkrVecNode* daughterNodeTds = convertTkrVecNode(daughterNodeRoot);

            // Add daughter to our current node
            if (daughterNodeTds) tkrVecNodeTds->push_back(daughterNodeTds);
        }
    }

    return tkrVecNodeTds;
}

int reconRootReaderAlg::makeSiblingMap(Event::TkrVecNode*        curNode, 
                                       Event::TkrNodeSiblingMap* siblingMap)
{
    //**************************************************************************************
    // NOTE: This is a copy of the same method that exists in TkrTreeBuilder! 
    //**************************************************************************************
    //
    // This method aims to set the "distance to the main branch" for each node in the tree
    // while it also builds out the "sibling map" which provides a list of all nodes at a given layer.
    // The "distance to the main branch" is the number of nodes from the nearest main branch.
    // A main branch is defined as the "first" branch, meaning that if you start with the head
    // node, then a "main" branch will be the first daughter in the list of nodes below the
    // current node. 

    // This is the non-recursive version

    // Set up a stack to handle the nodes that we will encounter
    // The idea is that we will visit each node first, the order
    // will be along the main branch first, then various branches 
    // off that as we proceed to process nodes in the stack.
    std::stack<Event::TkrVecNode*> nodeStack;

    // Seed it with the input node
    nodeStack.push(curNode);

    // Keep track of the number of links to the main branch
    int toMainBranch = 0;

    // Loop over stack until empty
    while(!nodeStack.empty())
    {
        // Recover pointer to node at the top of the stack
        Event::TkrVecNode* node = nodeStack.top();

        // Pop the top of the stack since we're done with this element
        nodeStack.pop();

        // Add this nodes daughters to the stack, provided we are not a leaf... 
        if (!node->empty())
        {
            // Add the daughters of this node to the stack
            // Note: need to add in reverse order so main branch will be accessed first
            for(Event::TkrVecNodeSet::reverse_iterator nodeItr = node->rbegin(); nodeItr != node->rend(); nodeItr++)
            {
                Event::TkrVecNode* nextNode = *nodeItr;

                nodeStack.push(nextNode);
            }
        }

        // If we are the head node then skip the rest
        if (!node->getAssociatedLink()) continue;

        // Find parent's distance from the main branch
        int parentDistance = node->getParentNode()->getBiLyrs2MainBrch();

        // Set the distance to the main branch
        node->setBiLyrs2MainBrch(parentDistance + toMainBranch);

        // While we are here, set the link to "associated" 
        const_cast<Event::TkrVecPointsLink*>(node->getAssociatedLink())->setAssociated();

        // Store this node in our sibling map 
        (*siblingMap)[node->getCurrentBiLayer()].push_back(node);

        // Otherwise we are at a leaf so we need to set the offset to get nonzero distances from the main branch
        if (node->empty()) toMainBranch = 1; // Note we set to one
    }

    return siblingMap->size();
}

StatusCode reconRootReaderAlg::storeTkrFilterParamsCol(TkrRecon* tkrRecRoot)
{
    // Purpose and Method:  Retrieve track and vertex collections from ROOT file and
    //  Store the transient versions on the TDS
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Is there a collection to deal with?
    if (tkrRecRoot->nTkrFilterParams())
    {
        // Create a collection for storage in the TDS
        Event::TkrFilterParamsCol* filterParamsColTds = new Event::TkrFilterParamsCol();
        // Loop through the entries
        for(int idx = 0; idx < tkrRecRoot->nTkrFilterParams(); idx++)
        {
            // Get a pointer to the root TkrFilterParams object
            TkrFilterParams* filterParamsRoot = tkrRecRoot->getTkrFilterParams(idx);

            // Convert to the TDS version
            Event::TkrFilterParams* filterParamsTds = convertTkrFilterParams(filterParamsRoot);

            // Store in TDS collection
            filterParamsColTds->push_back(filterParamsTds);

            // And keep track of relation between root and TDS
            m_common.m_rootTkrFilterParamsMap[filterParamsRoot] = filterParamsTds;
        }

        // Store in the TDS
        sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrFilterParamsCol, filterParamsColTds);
        if (sc.isFailure()) 
        {
            log << MSG::DEBUG;
            if( log.isActive()) log.stream() << "Failed to register TkrFilterParam";
            log << endreq;
            return StatusCode::FAILURE;
        }
    }

    return sc;
}

Event::TkrFilterParams* reconRootReaderAlg::convertTkrFilterParams(const TkrFilterParams* filterParamsRoot)
{
    // Must convert the position and direction from rootspeak to Gaudispeak
    Point  eventPos(filterParamsRoot->getEventPosition().x(),
        filterParamsRoot->getEventPosition().y(),
        filterParamsRoot->getEventPosition().z() );
    Vector eventAxis(filterParamsRoot->getEventAxis().x(),
        filterParamsRoot->getEventAxis().y(),
        filterParamsRoot->getEventAxis().z() );

    // Should now be able to populate the object
    Event::TkrFilterParams* filterParamsTds = new Event::TkrFilterParams(filterParamsRoot->getEventEnergy(), 
        eventPos,
        eventAxis);
    filterParamsTds->setStatusBit(filterParamsRoot->getStatusBits());
    filterParamsTds->setNumBiLayers(filterParamsRoot->getNumBiLayers());
    filterParamsTds->setNumIterations(filterParamsRoot->getNumIterations());
    filterParamsTds->setNumHitsTotal(filterParamsRoot->getNumHitsTotal());
    filterParamsTds->setNumDropped(filterParamsRoot->getNumDropped());
    filterParamsTds->setChiSquare(filterParamsRoot->getChiSquare());
    filterParamsTds->setAverageDistance(filterParamsRoot->getAverageDistance());
    filterParamsTds->setTransRms(filterParamsRoot->getTransRms());
    filterParamsTds->setLongRms(filterParamsRoot->getLongRms());
    filterParamsTds->setLongRmsAsym(filterParamsRoot->getLongRmsAsym());

    return filterParamsTds;
}

StatusCode reconRootReaderAlg::storeTkrEventParams(TkrRecon* tkrRecRoot)
{
    // Purpose and Method:  Retreive the TkrEventParams object from the 
    //                      TkrRecon object and store on the TDS
    // 
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Is there an object to store?
    if (tkrRecRoot->getTkrEventParams())
    {
        // Get the pointer in root
        TkrEventParams* eventParamsRoot = tkrRecRoot->getTkrEventParams();

        // Must convert the position and direction from rootspeak to Gaudispeak
        Point  eventPos(eventParamsRoot->getEventPosition().x(),
            eventParamsRoot->getEventPosition().y(),
            eventParamsRoot->getEventPosition().z() );
        Vector eventAxis(eventParamsRoot->getEventAxis().x(),
            eventParamsRoot->getEventAxis().y(),
            eventParamsRoot->getEventAxis().z() );

        // Should now be able to populate the object
        Event::TkrEventParams* eventParamsTds = new Event::TkrEventParams(eventParamsRoot->getEventEnergy(), 
            eventPos,
            eventAxis);
        eventParamsTds->setStatusBit(eventParamsRoot->getStatusBits());
        eventParamsTds->setNumBiLayers(eventParamsRoot->getNumBiLayers());
        eventParamsTds->setNumIterations(eventParamsRoot->getNumIterations());
        eventParamsTds->setNumHitsTotal(eventParamsRoot->getNumHitsTotal());
        eventParamsTds->setNumDropped(eventParamsRoot->getNumDropped());
        eventParamsTds->setChiSquare(eventParamsRoot->getChiSquare());
        eventParamsTds->setTransRms(eventParamsRoot->getTransRms());
        eventParamsTds->setLongRmsAve(eventParamsRoot->getLongRmsAve());

        // Store in the TDS
        sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrEventParams, eventParamsTds);
        if (sc.isFailure()) 
        {
            log << MSG::DEBUG;
            if( log.isActive()) log.stream() << "Failed to register TkrEventsParam";
            log << endreq;
            return StatusCode::FAILURE;
        }
    }

    return sc;
}

StatusCode reconRootReaderAlg::storeTrackAndVertexCol(
    TkrRecon *tkrRecRoot, bool vertexOnTdsFlag) 
{
    // Purpose and Method:  Retrieve track and vertex collections from ROOT file and
    //  Store the transient versions on the TDS
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Keep track of TkrFitTracks and the ROOT id so we can fill TkrVertex
    std::map<int, Event::TkrTrack*> trackMap;
    trackMap.clear();

    // Create TDS version of track collections
    Event::TkrTrackCol *trackTdsCol   = new Event::TkrTrackCol;
    Event::TkrTrackCol *crTrackTdsCol = new Event::TkrTrackCol;
    Event::TkrTrack    *trackTds    = 0;

    // Retrieve ROOT version of track collection
    const TObjArray *trackRootCol = tkrRecRoot->getTrackCol();
    TIter trackIter(trackRootCol);
    TObject *trackObj = 0;

    while ((trackObj = trackIter.Next())!=0) 
    {
        // int trkIdx = -1;
        TkrTrack* trackRoot = dynamic_cast<TkrTrack*>(trackObj);

        trackTds = convertTkrTrack(trackRoot);

        // Keep relation between Event and Root fit tracks
        m_common.m_rootTkrTrackMap[trackObj] = trackTds;

        trackTdsCol->push_back(trackTds);       
    }

    // ADW: Retrieve ROOT version of cosmic-ray track collection
    // Pretty ugly duplicate code.  Eventually move to map of collections.
    const TObjArray *crTrackRootCol = tkrRecRoot->getCRTrackCol();
    TIter crTrackIter(crTrackRootCol);
    trackObj = 0;

    while ((trackObj = crTrackIter.Next())!=0) 
    {
        // int trkIdx = -1;
        TkrTrack* trackRoot = dynamic_cast<TkrTrack*>(trackObj);

        trackTds = convertTkrTrack(trackRoot);

        // Keep relation between Event and Root fit tracks
        m_common.m_rootTkrTrackMap[trackObj] = trackTds;

        crTrackTdsCol->push_back(trackTds);
    }

    sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrTrackCol, trackTdsCol);
    if (sc.isFailure()) {

        log << MSG::DEBUG;
        if( log.isActive()) log.stream() << "Failed to register TkrTrackCol";
        log << endreq;
        return StatusCode::FAILURE;
    }

    sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrCRTrackCol, crTrackTdsCol);
    if (sc.isFailure()) {

        log << MSG::DEBUG;
        if( log.isActive()) log.stream() << "Failed to register TkrCRTrackCol";
        log << endreq;
        return StatusCode::FAILURE;
    }

    // If the vertex collection is not already on the TDS, proceed
    if (vertexOnTdsFlag == true) return sc;

    // Create the TDS version of the vertex collection
    Event::TkrVertexCol *vertexColTds = new Event::TkrVertexCol;

    // Retrieve ROOT version of vertex collection
    const TObjArray *vertexColRoot = tkrRecRoot->getVertexCol();
    TIter vertexIter(vertexColRoot);
    TkrVertex *vertexRoot = 0;

    while ((vertexRoot = (TkrVertex*)vertexIter.Next())!=0) 
    {
        commonRootData::TkrId hitId  =  vertexRoot->getTkrId();
        idents::TkrId         hitTdsId( hitId.getTowerX(),hitId.getTowerY(),hitId.getTray(),
            hitId.getBotTop(), hitId.getView());
        Point                 position( vertexRoot->getPosition().x(),
            vertexRoot->getPosition().y(),
            vertexRoot->getPosition().z() );

        Event::TkrVertex *vertexTds = new Event::TkrVertex(hitTdsId, 
            vertexRoot->getEnergy(),
            vertexRoot->getQuality(), 
            vertexRoot->getChiSquare(), 
            vertexRoot->getAddedRadLen(),
            vertexRoot->getDOCA(), 
            vertexRoot->getTkr1ArcLen(), 
            vertexRoot->getTkr2ArcLen(), 
            position.z(), 
            convertTkrTrackParams(vertexRoot->getVertexParams()));

        vertexTds->setStatusBit(vertexRoot->getStatusBits());

        // Keep relation between Event and Root vertices
        m_common.m_rootTkrVertexMap[vertexRoot] = vertexTds;

        unsigned int numTracks = vertexRoot->getNumTracks();
        unsigned int iTrack;
        for (iTrack = 0; iTrack < numTracks; iTrack++) 
        {
            const Event::TkrTrack* trackTds = m_common.m_rootTkrTrackMap[vertexRoot->getTrack(iTrack)];
            vertexTds->addTrack(trackTds);
        }

        vertexColTds->push_back(vertexTds);
    }

    sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrVertexCol, vertexColTds);
    if (sc.isFailure()) {

        log << MSG::DEBUG;
        if( log.isActive()) log.stream() << "Failed to register Tkr VertexCol";
        log << endreq;
        return StatusCode::FAILURE;
    }

    return sc;
}

/// Stores TkrTrack objects into Event::TkrFitTrack
Event::TkrTrack* reconRootReaderAlg::convertTkrTrack(const TkrTrack* trackRoot)
{
    Event::TkrTrack *trackTds = new Event::TkrTrack();

    // Convert these to TVector3's 
    Point  position( trackRoot->getInitialPosition().x(),
        trackRoot->getInitialPosition().y(),
        trackRoot->getInitialPosition().z());
    Vector direction(trackRoot->getInitialDirection().x(),
        trackRoot->getInitialDirection().y(),
        trackRoot->getInitialDirection().z());

    // Begin with filling the track member variables
    trackTds->setStatusBit(trackRoot->getStatusBits());
    trackTds->setInitialEnergy(trackRoot->getInitialEnergy());
    trackTds->setInitialPosition(position);
    trackTds->setInitialDirection(direction);
    trackTds->setChiSquareFilter(trackRoot->getChiSquareFilter());
    trackTds->setChiSquareSmooth(trackRoot->getChiSquareSmooth());
    trackTds->setNDegreesOfFreedom(trackRoot->getNDegreesOfFreedom());
    trackTds->setScatter(trackRoot->getScatter());
    trackTds->setQuality(trackRoot->getQuality());
    trackTds->setKalEnergy(trackRoot->getKalEnergy());
    trackTds->setKalEnergyError(trackRoot->getKalEnergyError());
    trackTds->setKalThetaMS(trackRoot->getKalThetaMS());
    trackTds->setNumXGaps(trackRoot->getNumXGaps());
    trackTds->setNumYGaps(trackRoot->getNumYGaps());
    trackTds->setNumXFirstGaps(trackRoot->getNumXFirstGaps());
    trackTds->setNumYFirstGaps(trackRoot->getNumYFirstGaps());
    trackTds->setNumSegmentPoints(trackRoot->getNumSegmentPoints());
    trackTds->setChiSqSegment(trackRoot->chiSquareSegment());
    trackTds->setNumXHits(trackRoot->getNumXHits());
    trackTds->setNumYHits(trackRoot->getNumYHits());
    trackTds->setTkrCalRadLen(trackRoot->getTkrCalRadlen());
    trackTds->setRangeEnergy(trackRoot->getRangeEnergy());

    // Now loop over the hit planes and fill that information
    TIterator* hitIter = trackRoot->Iterator();
    TObject*   hitObj  = 0;

    while ((hitObj = hitIter->Next())!=0)
    {
        TkrTrackHit*        trackHitRoot = dynamic_cast<TkrTrackHit*>(hitObj);
        Event::TkrTrackHit* trackHit     = convertTkrTrackHit(trackHitRoot);

        trackTds->push_back(trackHit);
    }

    return trackTds;
}
/// convert a ROOT TkrFitHit to a TDS Event::TkrFitHit
Event::TkrTrackParams reconRootReaderAlg::convertTkrTrackParams(const TkrTrackParams& paramsRoot)
{
    Event::TkrTrackParams params(paramsRoot.getxPosition(), paramsRoot.getxSlope(), 
        paramsRoot.getyPosition(), paramsRoot.getySlope(),
        paramsRoot.getxPosxPos(),  paramsRoot.getxPosxSlp(), 
        paramsRoot.getxPosyPos(),  paramsRoot.getxPosySlp(),
        paramsRoot.getxSlpxSlp(),  paramsRoot.getxSlpyPos(), 
        paramsRoot.getxSlpySlp(),  paramsRoot.getyPosyPos(), 
        paramsRoot.getyPosySlp(),  paramsRoot.getySlpySlp() );

    return params;
}

/// convert a ROOT TkrTrackHit to a TDS Event::TkrTrackHit
Event::TkrTrackHit* reconRootReaderAlg::convertTkrTrackHit(const TkrTrackHit* trackHitRoot)
{
    Event::TkrTrackHit* trackHitTds = new Event::TkrTrackHit();
    commonRootData::TkrId hitId  = trackHitRoot->getTkrId();

    // Fill cluster id only if it exists
    if (const TkrCluster* clusRoot   = trackHitRoot->getClusterPtr())
    {
        const Event::TkrCluster* clusTds = m_common.m_rootTkrClusterMap[clusRoot];
        trackHitTds->setClusterPtr(clusTds);
    }

    if(hitId.hasTowerX()&&hitId.hasTowerY()&&hitId.hasTray()&&hitId.hasBotTop()) {
        if(hitId.hasView()) {
            idents::TkrId hitTdsId(hitId.getTowerX(), hitId.getTowerY(),
                hitId.getTray(), hitId.getBotTop(), hitId.getView());
            trackHitTds->setTkrID(hitTdsId);
        } else { // no view if no cluster
            idents::TkrId hitTdsId(hitId.getTowerX(), hitId.getTowerY(),
                hitId.getTray(), hitId.getBotTop());
            trackHitTds->setTkrID(hitTdsId);
        }
    }
    else {
        idents::TkrId hitTdsId = idents::TkrId();
        trackHitTds->setTkrID(hitTdsId);
    }

    // Fill in the "easy" stuff first 
    trackHitTds->setStatusBit(trackHitRoot->getStatusBits());
    trackHitTds->setZPlane(trackHitRoot->getZPlane());
    trackHitTds->setEnergy(trackHitRoot->getEnergy());
    trackHitTds->setRadLen(trackHitRoot->getRadLen());
    trackHitTds->setActiveDist(trackHitRoot->getActiveDist());
    trackHitTds->setChiSquareFilter(trackHitRoot->getChiSquareFilter());
    trackHitTds->setChiSquareSmooth(trackHitRoot->getChiSquareSmooth());

    // Now set the track params
    if (trackHitRoot->validMeasuredHit()) trackHitTds->getTrackParams(Event::TkrTrackHit::MEASURED) = 
        convertTkrTrackParams(trackHitRoot->getTrackParams(TkrTrackHit::MEASURED));
    if (trackHitTds->validPredictedHit()) trackHitTds->getTrackParams(Event::TkrTrackHit::PREDICTED) = 
        convertTkrTrackParams(trackHitRoot->getTrackParams(TkrTrackHit::PREDICTED));
    if (trackHitTds->validFilteredHit())  trackHitTds->getTrackParams(Event::TkrTrackHit::FILTERED) = 
        convertTkrTrackParams(trackHitRoot->getTrackParams(TkrTrackHit::FILTERED));
    if (trackHitTds->validFilteredHit())  trackHitTds->getTrackParams(Event::TkrTrackHit::REVFIT)   = 
        convertTkrTrackParams(trackHitRoot->getTrackParams(TkrTrackHit::REVFIT)); 
    if (trackHitTds->validSmoothedHit())  trackHitTds->getTrackParams(Event::TkrTrackHit::SMOOTHED) = 
        convertTkrTrackParams(trackHitRoot->getTrackParams(TkrTrackHit::SMOOTHED));
    if (trackHitTds->validRevFit())       trackHitTds->getTrackParams(Event::TkrTrackHit::REVFIT) = 
        convertTkrTrackParams(trackHitRoot->getTrackParams(TkrTrackHit::REVFIT));
    if (trackHitTds->validMaterial())     trackHitTds->getTrackParams(Event::TkrTrackHit::QMATERIAL) = 
        convertTkrTrackParams(trackHitRoot->getTrackParams(TkrTrackHit::QMATERIAL));

    return trackHitTds;
}

StatusCode reconRootReaderAlg::storeTkrTruncationInfo(TkrRecon *tkrRecRoot) {
    StatusCode sc = StatusCode::SUCCESS;
    const TObjArray *truncationDataColRoot = tkrRecRoot->getTruncationDataCol();
    Event::TkrTruncationInfo* truncationInfoTds = new Event::TkrTruncationInfo();
    RootPersistence::convert(*truncationDataColRoot,*truncationInfoTds) ;
    sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrTruncationInfo, truncationInfoTds);    
    return sc;    
}

StatusCode reconRootReaderAlg::readCalRecon() {
    // Purpose and Method:: Read in CAL recon data from ROOT and store on the TDS

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    CalRecon *calRecRoot = m_reconEvt->getCalRecon();
    if (!calRecRoot) return StatusCode::SUCCESS;

    // Make sure the CalRecon branch exists on the TDS
    DataObject* pnode =0;
    sc = eventSvc()->retrieveObject(EventModel::CalRecon::Event, pnode);
    if( sc.isFailure() ) 
    {
        sc = eventSvc()->registerObject(EventModel::CalRecon::Event,new DataObject);
        if( sc.isFailure() ) 
        {
            log << MSG::ERROR << "Could not create CalRecon directory" << endreq;
            return sc;
        }
    }

    SmartDataPtr<Event::CalXtalRecCol> xtalRecColTds(eventSvc(),EventModel::CalRecon::CalXtalRecCol);
    if (xtalRecColTds){
        log << MSG::INFO << "XtalRecCol data is already on the TDS" << endreq;
    } else {
        sc = storeCalXtalRecDataCol(calRecRoot);
    }

    if (sc.isFailure()) {
        log << MSG::INFO << "Failed to store CalXtalRecCol on the TDS" << endreq;
        return sc;
    }

    SmartDataPtr<Event::CalEventEnergyCol> checkCalEventEnergyColTds(eventSvc(), EventModel::CalRecon::CalEventEnergyCol);
    if (checkCalEventEnergyColTds) {
        log << MSG::INFO << "CalEventEnergy data is already on the TDS" << endreq;
    } else {
        sc = storeCalEventEnergyCol(calRecRoot);
    }
    if (sc.isFailure()) {
        log << MSG::INFO << "Failed to store CalEventEnergy on the TDS" << endreq;
        return sc;
    }

    // Look up the CalClusterCol, which is the "owner" of the CalCluster objects in the TDS
    SmartDataPtr<Event::CalClusterCol> checkCalClusterColTds(eventSvc(),EventModel::CalRecon::CalClusterCol);

    // If they are already there then this must be some sort of error...
    if (checkCalClusterColTds)
    {
        log << MSG::INFO << "CalClusterCol data is already on the TDS" << endreq;
        return sc;
    } else {
        sc = storeCalClusterCol(calRecRoot);
    }

    // New school: Cal Clusters are stored in the CalClusterMap
    SmartDataPtr<Event::CalClusterMap> checkCalClusterMapTds(eventSvc(),EventModel::CalRecon::CalClusterMap);
    if (checkCalClusterMapTds){
        log << MSG::INFO << "CalClusterMap data is already on the TDS" << endreq;
        return sc;
    } else {
        sc = storeCalClusterMap(calRecRoot);
    }

    SmartDataPtr<Event::CalEventEnergyMap> checkCalEventEnergyMapTds(eventSvc(), EventModel::CalRecon::CalEventEnergyMap);
    if (checkCalEventEnergyMapTds) {
        log << MSG::INFO << "CalEventEnergyMap data is already on the TDS" << endreq;
    } else {
        sc = storeCalEventEnergyMap(calRecRoot);
    }
    if (sc.isFailure()) {
        log << MSG::INFO << "Failed to store CalEventEnergyMap on the TDS" << endreq;
        return sc;
    }

    SmartDataPtr<Event::CalMipTrackCol> checkCalMipTrackColTds(eventSvc(),EventModel::CalRecon::CalMipTrackCol);
    if (checkCalMipTrackColTds){
        log << MSG::INFO << "CalMipTrackCol data is already on the TDS" << endreq;
        return sc;
    } else {
        sc = storeCalMipTrackCol(calRecRoot);
    }


    return sc;
}


StatusCode reconRootReaderAlg::storeCalXtalRecDataCol(CalRecon *calRecRoot) {
    // Purpose and Method:  Retrieve CAL xtal recon data from ROOT and store on the TDS
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Create new TDS xtalRec collection
    Event::CalXtalRecCol *calXtalRecColTds = new Event::CalXtalRecCol();
    // Retrieve ROOT CAL xtal Rec collection
    const TObjArray *calXtalRecColRoot = calRecRoot->getCalXtalRecCol();
    TIter calXtalIter(calXtalRecColRoot);
    const CalXtalRecData *calXtalRecRoot = 0;
    while ((calXtalRecRoot = (CalXtalRecData*)calXtalIter.Next())!=0) {

        if (calXtalRecRoot->getMode() == CalXtalId::ALLRANGE) {
            unsigned int range ;
            for ( range = idents::CalXtalId::LEX8 ;
                range < idents::CalXtalId::HEX1 ;
                ++range ) {    
                    if (!calXtalRecRoot->getRangeRecData(range)) {
                        log << MSG::DEBUG;
                        if( log.isActive()) log.stream() << "Readout for Range " << range << " does not exist";
                        log << endreq;
                    }
            }
        }

        Event::CalXtalRecData * calXtalRecDataTds
            = new Event::CalXtalRecData ;
        RootPersistence::convert(*calXtalRecRoot,*calXtalRecDataTds) ;

        calXtalRecColTds->push_back(calXtalRecDataTds);

        m_common.m_rootCalXtalRecDataMap[calXtalRecRoot] = calXtalRecDataTds;
    }

    //register output data collection as a TDS object
    sc = eventSvc()->registerObject(EventModel::CalRecon::CalXtalRecCol, calXtalRecColTds);
    if (sc.isFailure()) {

        log << MSG::DEBUG;
        if( log.isActive()) log.stream() << "Failed to register CalXtalRecCol";
        log << endreq;
        return StatusCode::FAILURE;
    }

    return sc;
}

StatusCode reconRootReaderAlg::storeCalClusterMap(CalRecon *calRecRoot) 
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // First recover the pointer to the TMap
    const TMap* calClusterMapRoot = calRecRoot->getCalClusterMap();

    // Get a shiny new TDS version of the CalClusterMap
    Event::CalClusterMap *calClusterMapTds = new Event::CalClusterMap();

    // Set up an iterator to go through the entries
    TMapIter calClusterMapIter(calClusterMapRoot);

    // As we iterate through the map, we'll be obtaining the "next" key
    const TObject* nextKey = 0;

    // Ok, now iterate through the input map
    while((nextKey = calClusterMapIter.Next()) != 0)
    {
        // Recover the key string
        const TObjString* keyRoot = (TObjString*)nextKey;

        // Convert this to a std::string
        std::string keyTds(keyRoot->GetName());

        // Recover the TObjArray which will hold the Cal Clusters
        TObjArray* clusterVecRoot = (TObjArray*)calClusterMapRoot->GetValue(nextKey);

        // Make sure something is there
        if (clusterVecRoot->GetEntries() > 0)
        {
            // Get an iterator over this array
            TObjArrayIter calClusterVecIter(clusterVecRoot);

            CalCluster *clusterRoot = 0;

            while ((clusterRoot = (CalCluster*)calClusterVecIter.Next())!=0) 
            {        
                // Check to be sure there is an entry in the root to TDS map
                if (m_common.m_rootCalClusterMap.find(clusterRoot) != m_common.m_rootCalClusterMap.end()) 
                {
                    // Recover the pointer to the cluster in question
                    Event::CalCluster* clusterTds = const_cast<Event::CalCluster*>(m_common.m_rootCalClusterMap[clusterRoot]);

                    // Add this to the CalClusterMap
                    (*calClusterMapTds)[keyTds].push_back(clusterTds) ;
                }
                else
                {
                    int whatshouldiputhere = 0;
                }
            }
        }
    }

    if (!calClusterMapTds->empty()) sc = eventSvc()->registerObject(EventModel::CalRecon::CalClusterMap, calClusterMapTds);
    else                            delete calClusterMapTds;

    return sc;
}

StatusCode reconRootReaderAlg::storeCalClusterCol(CalRecon *calRecRoot) {
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    const TObjArray *calClusterColRoot = calRecRoot->getCalClusterCol();
    TIter calClusterIter(calClusterColRoot);
    CalCluster *calClusterRoot = 0;

    Event::CalClusterCol *calClusterColTds = new Event::CalClusterCol();

    while ((calClusterRoot = (CalCluster*)calClusterIter.Next())!=0) {        
        Event::CalCluster * calClusterTds = new Event::CalCluster() ;
        RootPersistence::convert(*calClusterRoot,*calClusterTds) ;
        calClusterColTds->push_back(calClusterTds) ;

        m_common.m_rootCalClusterMap[calClusterRoot] = calClusterTds;
    }

    sc = eventSvc()->registerObject(EventModel::CalRecon::CalClusterCol, calClusterColTds);

    return sc;
}

StatusCode reconRootReaderAlg::storeCalMipTrackCol(CalRecon *calRecRoot) 
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    const TObjArray *calMipTrackColRoot = calRecRoot->getCalMipTrackCol();
    TIter calMipTrackIter(calMipTrackColRoot);
    CalMipTrack *calMipTrackRoot = 0;

    Event::CalMipTrackCol *calMipTrackColTds = new Event::CalMipTrackCol();

    while ((calMipTrackRoot = (CalMipTrack*)calMipTrackIter.Next())!=0) 
    {        
        Event::CalMipTrack* calMipTrackTds = new Event::CalMipTrack() ;
        RootPersistence::convert(*calMipTrackRoot,*calMipTrackTds) ;
        calMipTrackColTds->push_back(calMipTrackTds) ;
    }

    sc = eventSvc()->registerObject(EventModel::CalRecon::CalMipTrackCol, calMipTrackColTds);

    return sc;
}

StatusCode reconRootReaderAlg::storeCalEventEnergyCol(CalRecon *calRecRoot) {

    // David C. : currently, there is only one CalEventEnergy in the TDS,
    // yet, I prefered to consider CalEventEnergy as a usual objet on
    // the ROOT side, that is why in the root tree there is a collection
    // of CalEventEnergy. This collection will always have a single element,
    // until a change is made in the TDS.

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    //    const TObjArray * calEventEnergyColRoot = calRecRoot->getCalEventEnergyCol();
    //    if (calEventEnergyColRoot->GetEntries()>1) {
    //        // this should not happen !!
    //        log<<MSG::ERROR ;
    //        if (log.isActive()) log.stream()<<"Several CalEventEnergy in ROOT file" ;
    //        log<<endreq ;
    //        return StatusCode::FAILURE;
    //    }


    //    Event::CalEventEnergyCol * calEventEnergyColTds = new Event::CalEventEnergyCol();
    //    TIter calEventEnergyIter(calEventEnergyColRoot) ;
    //    CalEventEnergy * calEventEnergyRoot = 0 ;
    //    while ((calEventEnergyRoot = (CalEventEnergy*)calEventEnergyIter.Next())!=0) {        
    //        Event::CalEventEnergy * calEventEnergyTds = new Event::CalEventEnergy() ;
    //        RootPersistence::convert(*calEventEnergyRoot,*calEventEnergyTds) ;
    //        calEventEnergyColTds->push_back(calEventEnergyTds) ;

    //        m_common.m_rootCalEventEnergyMap[calEventEnergyRoot] = calEventEnergyTds;
    //    }

    //    sc = eventSvc()->registerObject(EventModel::CalRecon::CalEventEnergyCol, calEventEnergyColTds);

    return sc;

}

StatusCode reconRootReaderAlg::storeCalEventEnergyMap(CalRecon *calRecRoot) 
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // First recover the pointer to the TMap
    const TMap* calEventEnergyMapRoot = calRecRoot->getCalEventEnergyMap();

    // Get a shiny new TDS version of the CalClusterMap
    Event::CalEventEnergyMap *calEventEnergyMapTds = new Event::CalEventEnergyMap();

    // Set up an iterator to go through the entries
    TMapIter calEventEnergyMapIter(calEventEnergyMapRoot);

    // As we iterate through the map, we'll be obtaining the "next" key
    const TObject* nextKey = 0;

    // Ok, now iterate through the input map
    while((nextKey = calEventEnergyMapIter.Next()) != 0)
    {
        // Recover the key cluster
        const CalCluster* clusterRoot = (CalCluster*)nextKey;

        // Recover the TDS version of the cluster
        if (m_common.m_rootCalClusterMap.find(clusterRoot) != m_common.m_rootCalClusterMap.end())
        {
            Event::CalCluster* clusterTds = const_cast<Event::CalCluster*>(m_common.m_rootCalClusterMap[clusterRoot]);

            // Recover the TObjArray which will hold the CalEventEnergy objects
            TObjArray* eventEnergyVecRoot = (TObjArray*)calEventEnergyMapRoot->GetValue(nextKey);

            // Make sure something is there
            if (eventEnergyVecRoot->GetEntries() > 0)
            {
                // Get an iterator over this array
                TObjArrayIter calEventEnergyVecIter(eventEnergyVecRoot);

                CalEventEnergy* calEventEnergyRoot = 0;

                while ((calEventEnergyRoot = (CalEventEnergy*)calEventEnergyVecIter.Next())!=0) 
                {        
                    Event::CalEventEnergy * calEventEnergyTds = new Event::CalEventEnergy() ;
                    RootPersistence::convert(*calEventEnergyRoot,*calEventEnergyTds) ;

                    // Add this to the CalClusterMap
                    (*calEventEnergyMapTds)[clusterTds].push_back(calEventEnergyTds) ;
                }
            }
        }
    }

    sc = eventSvc()->registerObject(EventModel::CalRecon::CalEventEnergyMap, calEventEnergyMapTds);

    return sc;
}


StatusCode reconRootReaderAlg::readAcdRecon() {
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    const AcdReconV2 *acdRecRootV2 = m_reconEvt->getAcdReconV2();
    if (!acdRecRootV2) {

        log << MSG::DEBUG;
        if( log.isActive()) log.stream() << "No AcdReconV2 found in ROOT file";
        log << endreq;
        return StatusCode::SUCCESS;
    }

    Event::AcdReconV2 * acdRecTdsV2 = new Event::AcdReconV2();    
    RootPersistence::convert(*acdRecRootV2,*acdRecTdsV2) ;

    sc = eventSvc()->registerObject(EventModel::AcdReconV2::Event, acdRecTdsV2);
    if (sc.isFailure()) {

        log << MSG::DEBUG;
        if( log.isActive()) log.stream() << "Failed to register AcdReconV2";
        log << endreq;
        return StatusCode::FAILURE;
    }

    return sc;    
}




StatusCode reconRootReaderAlg::readAdfRecon()
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    const reconRootData::AdfRecon *adfRecRoot = m_reconEvt->getAdfRecon();
    if (!adfRecRoot) 
        return StatusCode::SUCCESS;

    AncillaryData::Recon *adfRecTds = new AncillaryData::Recon();
    RootPersistence::convert(*adfRecRoot, *adfRecTds);

    sc = eventSvc()->registerObject("/Event/AncillaryEvent/Recon", adfRecTds);

    if (sc.isFailure()) {

        log << MSG::ERROR << "Failed to register AdfRecon" << endreq;
        return StatusCode::FAILURE;
    }

    return sc;    
}

StatusCode reconRootReaderAlg::readCalTkrAcdRelations()
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // The Tree to Cluster relations are stored in the TkrRecon branch
    TkrRecon *tkrRecRoot = m_reconEvt->getTkrRecon();
    if(!tkrRecRoot) return sc;

    sc = storeTreeClusterRelations(tkrRecRoot);

    return sc;
}

StatusCode reconRootReaderAlg::storeTreeClusterRelations(TkrRecon* tkrRecRoot)
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    Event::TreeClusterRelationCol* relationTdsCol          = new Event::TreeClusterRelationCol();
    Event::TreeToRelationMap*      treeToRelationTdsMap    = new Event::TreeToRelationMap();
    Event::ClusterToRelationMap*   clusterToRelationTdsMap = new Event::ClusterToRelationMap();

    const TObjArray * relationRootCol = tkrRecRoot->getTreeClusterRelationCol();
    TIter relationIter(relationRootCol);
    TreeClusterRelation* relationRoot = 0;

    while ((relationRoot = (TreeClusterRelation*)relationIter.Next())!=0) 
    {
        // Recover the pointers to the root versions of the related objects
        Event::TkrTree*    tree    = const_cast<Event::TkrTree*>(m_common.m_rootTkrTreeMap[relationRoot->getTree()]);
        Event::CalCluster* cluster = const_cast<Event::CalCluster*>(m_common.m_rootCalClusterMap[relationRoot->getCluster()]);

        // Create a new TDS version
        Event::TreeClusterRelation* relationTds = new Event::TreeClusterRelation(tree, 
            cluster, 
            relationRoot->getTreeClusDoca(),
            relationRoot->getTreeClusCosAngle(),
            relationRoot->getTreeClusDistAtZ(),
            relationRoot->getClusEnergy() );

        relationTdsCol->push_back(relationTds);
        (*treeToRelationTdsMap)[tree].push_back(relationTds);
        (*clusterToRelationTdsMap)[cluster].push_back(relationTds);
    }

    // First we need to follow through on some craziness to create our subdirectory...
    DataObject* pnode =0;

    if( (eventSvc()->retrieveObject(EventModel::Recon::Event, pnode)).isFailure() ) 
    {
        sc = eventSvc()->registerObject(EventModel::Recon::Event, new DataObject);
        if( sc.isFailure() ) 
        {
            log << MSG::ERROR << "Could not create Recon directory in storeTreeClusterRelations" 
                << endreq;
            return sc;
        }
    }

    sc = eventSvc()->registerObject(EventModel::Recon::TreeClusterRelationCol, relationTdsCol);
    if (sc.isFailure()) 
    {
        log << MSG::DEBUG;
        if( log.isActive()) log.stream() << "Failed to register TreeClusterRelationCol";
        log << endreq;
        return StatusCode::FAILURE;
    }

    sc = eventSvc()->registerObject(EventModel::Recon::TreeToRelationMap, treeToRelationTdsMap);
    if (sc.isFailure()) 
    {
        log << MSG::DEBUG;
        if( log.isActive()) log.stream() << "Failed to register TreeClusterRelationCol";
        log << endreq;
        return StatusCode::FAILURE;
    }

    sc = eventSvc()->registerObject(EventModel::Recon::ClusterToRelationMap, clusterToRelationTdsMap);
    if (sc.isFailure()) 
    {
        log << MSG::DEBUG;
        if( log.isActive()) log.stream() << "Failed to register TreeClusterRelationCol";
        log << endreq;
        return StatusCode::FAILURE;
    }

    return sc;
}


void reconRootReaderAlg::close() 
{
    // Purpose and Method:  Writes the ROOT file at the end of the run.
    //    The TObject::kOverWrite parameter is used in the Write method
    //    since ROOT will periodically write to the ROOT file when the bufSize
    //    is filled.  Writing would create 2 copies of the same tree to be
    //    stored in the ROOT file, if we did not specify kOverwrite.

}

void reconRootReaderAlg::endEvent() {
    if (m_reconEvt)  m_reconEvt->Clear(m_clearOption.c_str()) ; 
    m_reconEvt = 0 ;
}

StatusCode reconRootReaderAlg::finalize()
{
    close();

    StatusCode sc = StatusCode::SUCCESS;
    //setFinalized(); No longer available in Gaudi v21r7
    return sc;
}

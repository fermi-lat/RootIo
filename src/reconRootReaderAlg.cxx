#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "idents/CalXtalId.h"
#include "Event/Recon/AcdRecon/AcdRecon.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/CalRecon/CalCluster.h"   
#include "Event/Recon/CalRecon/CalXtalRecData.h"  
#include "Event/Recon/CalRecon/CalEventEnergy.h"  

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
#include "RootConvert/Recon/CalClusterConvert.h"
#include "RootConvert/Recon/CalEventEnergyConvert.h"
#include "RootConvert/Recon/CalMipTrackConvert.h"
#include "RootConvert/Recon/CalXtalRecDataConvert.h"
#include "RootConvert/Recon/AcdReconConvert.h"
#include "RootConvert/Recon/AdfReconConvert.h"

#include <vector>
#include <map>

/** @class reconRootReaderAlg
* @brief Reads Reconstruction data from a persistent ROOT file and stores the
* the data in the TDS.
*
* @author Heather Kelly
* $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/reconRootReaderAlg.cxx,v 1.80 2007/05/11 23:00:02 usher Exp $
*/

class reconRootReaderAlg : public Algorithm
{   
public:
    
    reconRootReaderAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    /// Handles setup by opening ROOT file in read mode and creating a new TTree
    StatusCode initialize();
    
    /// Orchastrates reading from ROOT file and storing the data on the TDS for each event
    StatusCode execute();
    
    /// Closes the ROOT file and cleans up
    StatusCode finalize();
    
private:
    
    /// Reads top-level DigiEvent
    StatusCode readReconEvent();
    
    /// Reads TKR recon data from ROOT and puts data on TDS
    StatusCode readTkrRecon();
    
    StatusCode storeTkrClusterCol(TkrRecon *tkrRecRoot);
    
    StatusCode storeTrackAndVertexCol(TkrRecon *tkrRecRoot, bool vertexOnTdsFlag);
    
    /// convert a ROOT TkrFitHit to a TDS Event::TkrFitHit
    Event::TkrTrackParams convertTkrTrackParams(const TkrTrackParams& paramsRoot);
    /// convert a ROOT TkrTrackHit to a TDS Event::TkrTrackHit
    Event::TkrTrackHit* convertTkrTrackHit(const TkrTrackHit* trackHitRoot);
    /// Stores TkrTrack objects into Event::TkrFitTrack
    Event::TkrTrack* convertTkrTrack(const TkrTrack* trackRoot);
    
    /// Reads CAL recon data from ROOT and puts data on TDS
    StatusCode readCalRecon();
    
    /// read in CAL xtal recon data from ROOT and store on the TDS
    StatusCode storeCalXtalRecDataCol(CalRecon *calRecRoot);
    
    /// read CAL cluster data from ROOT and store on TDS
    StatusCode storeCalClusterCol(CalRecon *calRecRoot);
    
    /// read CAL cluster data from ROOT and store on TDS
    StatusCode storeCalMipTrackCol(CalRecon *calRecRoot);

    /// read CAL eventEnergy  data from ROOT and store on TDS
    StatusCode storeCalEventEnergyCol(CalRecon *calRecRoot);
    
    /// Reads ACD recon data from ROOT and puts data on the TDS
    StatusCode readAcdRecon();

    StatusCode readAdfRecon();
    
    /// Closes the ROOT file
    void close();
    
    /// Top-level Monte Carlo ROOT object
    ReconEvent *m_reconEvt;
//    /// name of the input ROOT file
//    std::string m_fileName;
    /// List of input files
    StringArrayProperty m_fileList;
    /// name of the Recon TTree stored in the ROOT file
    std::string m_treeName;
    std::string m_branchName;
    /// Option string which will be passed to McEvent::Clear
    std::string m_clearOption;

    commonData m_common;

    IRootIoSvc*   m_rootIoSvc;

    bool FixAcdStreamerDone;

    // ADDED FOR THE FILE HEADERS DEMO
    IFhTool * m_headersTool ;
};

static const AlgFactory<reconRootReaderAlg>  Factory;
const IAlgFactory& reconRootReaderAlgFactory = Factory;


reconRootReaderAlg::reconRootReaderAlg(const std::string& name, ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator), m_reconEvt(0)
{
    // Input pararmeters that may be set via the jobOptions file
    // Input ROOT file name
//    declareProperty("reconRootFile",m_fileName="");
    StringArrayProperty initList;
    std::vector<std::string> initVec;
    initVec.push_back("recon.root");
    initList.setValue(initVec);
    declareProperty("reconRootFileList", m_fileList=initList);
    initVec.clear();
    // Input TTree name
    declareProperty("reconTreeName", m_treeName="Recon");
    declareProperty("reconBranchName", m_branchName="ReconEvent");
    declareProperty("clearOption", m_clearOption="");
    
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
    }

    // Set up new school system...
    std::string type = "RECON";
    m_rootIoSvc->prepareRootInput(type, m_treeName, m_branchName, m_fileList);

    return sc;
    
}

StatusCode reconRootReaderAlg::execute()
{
    // Purpose and Method:  Called once per event.  This method calls
    //   the appropriate methods to read data from the ROOT file and store
    //   data on the TDS.
    
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    if (m_reconEvt)
     {
      m_reconEvt->Clear(m_clearOption.c_str()) ; 
      m_reconEvt = 0 ;
     }

    // Try reading the event this way... 
    std::string type = "RECON";
    m_reconEvt = dynamic_cast<ReconEvent*>(m_rootIoSvc->getNextEvent(type)) ;

    if (!m_reconEvt) return StatusCode::FAILURE;

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

        Event::TkrCluster* clusterTds 
            = new Event::TkrCluster(tkrId, 
            clusterRoot->getFirstStrip(),
            clusterRoot->getLastStrip(),
            posTds,
            (int)clusterRoot->getRawToT(),
            clusterRoot->getMips(),
            clusterRoot->getStatusWord(),
            clusterRoot->getNBad()
    );

        clusterTdsCol->push_back(clusterTds);
        (*clusMap)[tkrId].push_back(clusterTds);
 
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

StatusCode reconRootReaderAlg::storeTrackAndVertexCol(TkrRecon *tkrRecRoot, bool vertexOnTdsFlag) {
    // Purpose and Method:  Retrieve track and vertex collections from ROOT file and
    //  Store the transient versions on the TDS
    
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    // Keep track of TkrFitTracks and the ROOT id so we can fill TkrVertex
    std::map<int, Event::TkrTrack*> trackMap;
    trackMap.clear();
    
    // Create TDS version of track collection
    Event::TkrTrackCol *trackTdsCol = new Event::TkrTrackCol;
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
    
    sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrTrackCol, trackTdsCol);
    if (sc.isFailure()) {
        
        log << MSG::DEBUG;
        if( log.isActive()) log.stream() << "Failed to register TkrTrackCol";
        log << endreq;
        return StatusCode::FAILURE;
    }
    
    // If the vertex collection is not already on the TDS, continue
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
    if (trackHitTds->validMaterial())     trackHitTds->getTrackParams(Event::TkrTrackHit::QMATERIAL) = 
               convertTkrTrackParams(trackHitRoot->getTrackParams(TkrTrackHit::QMATERIAL));

    return trackHitTds;
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
    
    SmartDataPtr<Event::CalClusterCol> checkCalClusterColTds(eventSvc(),EventModel::CalRecon::CalClusterCol);
    if (checkCalClusterColTds){
        log << MSG::INFO << "CalClusterCol data is already on the TDS" << endreq;
        return sc;
    } else {
        sc = storeCalClusterCol(calRecRoot);
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
    
    const TObjArray * calEventEnergyColRoot = calRecRoot->getCalEventEnergyCol();
//    if (calEventEnergyColRoot->GetEntries()>1) {
//        // this should not happen !!
//        log<<MSG::ERROR ;
//        if (log.isActive()) log.stream()<<"Several CalEventEnergy in ROOT file" ;
//        log<<endreq ;
//        return StatusCode::FAILURE;
//    }
        
    
    Event::CalEventEnergyCol * calEventEnergyColTds = new Event::CalEventEnergyCol();
    TIter calEventEnergyIter(calEventEnergyColRoot) ;
    CalEventEnergy * calEventEnergyRoot = 0 ;
    while ((calEventEnergyRoot = (CalEventEnergy*)calEventEnergyIter.Next())!=0) {        
//    if ((calEventEnergyRoot = (CalEventEnergy*)calEventEnergyIter.Next())!=0) {        
        Event::CalEventEnergy * calEventEnergyTds = new Event::CalEventEnergy() ;
        RootPersistence::convert(*calEventEnergyRoot,*calEventEnergyTds) ;
        calEventEnergyColTds->push_back(calEventEnergyTds) ;
//        sc = eventSvc()->registerObject(EventModel::CalRecon::CalEventEnergy, calEventEnergyTds);
    }
    
    sc = eventSvc()->registerObject(EventModel::CalRecon::CalEventEnergyCol, calEventEnergyColTds);
    
    return sc;

}


StatusCode reconRootReaderAlg::readAcdRecon() {
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    const AcdRecon *acdRecRoot = m_reconEvt->getAcdRecon();
    if (!acdRecRoot) {
        
        log << MSG::DEBUG;
        if( log.isActive()) log.stream() << "No AcdRecon found in ROOT file";
        log << endreq;
        return StatusCode::SUCCESS;
    }

    Event::AcdRecon * acdRecTds = new Event::AcdRecon();    
    RootPersistence::convert(*acdRecRoot,*acdRecTds) ;

    sc = eventSvc()->registerObject(EventModel::AcdRecon::Event, acdRecTds);
    if (sc.isFailure()) {
        
        log << MSG::DEBUG;
        if( log.isActive()) log.stream() << "Failed to register AcdRecon";
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


void reconRootReaderAlg::close() 
{
    // Purpose and Method:  Writes the ROOT file at the end of the run.
    //    The TObject::kOverWrite parameter is used in the Write method
    //    since ROOT will periodically write to the ROOT file when the bufSize
    //    is filled.  Writing would create 2 copies of the same tree to be
    //    stored in the ROOT file, if we did not specify kOverwrite.
    
    // skip deletion of chain due to seg fault with TChainIndex
    //if (m_reconTree) {
    //    delete m_reconTree;
    //    m_reconTree = 0;
    //}
}

StatusCode reconRootReaderAlg::finalize()
{
    close();
    
    StatusCode sc = StatusCode::SUCCESS;
    setFinalized();
    return sc;
}

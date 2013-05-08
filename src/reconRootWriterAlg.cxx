#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"

#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrTree.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/TkrRecon/TkrFilterParams.h"
#include "Event/Recon/TkrRecon/TkrEventParams.h"

#include "Event/Recon/TreeClusterRelation.h"

#include "Event/Recon/CalRecon/CalClusterMap.h"

#include "LdfEvent/EventSummaryData.h"

#include "AncillaryDataEvent/Recon.h"

#include "idents/CalXtalId.h"

#include "facilities/Util.h"
#include "commonData.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TVector3.h"
#include "TMatrixD.h"

#include "reconRootData/ReconEvent.h"
#include "RootIo/IRootIoSvc.h"

// ADDED FOR THE FILE HEADERS DEMO
#include "RootIo/FhTool.h"

// low level converters
#include "RootConvert/Recon/TkrTruncationInfoConvert.h"

#include "RootConvert/Recon/CalClusterConvert.h"
#include "RootConvert/Recon/CalXtalRecDataConvert.h"   
#include "RootConvert/Recon/CalEventEnergyConvert.h"
#include "RootConvert/Recon/CalMipTrackConvert.h"
#include "RootConvert/Recon/AdfReconConvert.h"
#include "RootConvert/Recon/AcdReconConvert.h"

#include "RootConvert/Recon/GcrXtalConvert.h"
#include "RootConvert/Recon/GcrTrackConvert.h"

#include <cstdlib>

/** @class reconRootWriterAlg
* @brief Writes Recon TDS data to a persistent ROOT file.
*
* @author Heather Kelly and Tracy Usher
* $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/RootIo/src/reconRootWriterAlg.cxx,v 1.103 2013/02/19 04:24:51 usher Exp $
*/

class reconRootWriterAlg : public Algorithm
{   
public:
    
    reconRootWriterAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    /// Handles setup by opening ROOT file in write mode and creating a new TTree
    StatusCode initialize();
    
    /// Orchastrates reading from TDS and writing to ROOT for each event
    StatusCode execute();
    
    /// Closes the ROOT file and cleans up
    StatusCode finalize();
    
private:
    
    /// Retrieves event Id and run Id from TDS and fills the Recon ROOT object
    StatusCode writeReconEvent();
    
    /// Retrieves the TKR reconstruction data from the TDS and fills the TkrRecon
    /// ROOT object
    StatusCode writeTkrRecon();
    
    /// These are the methods specific to filling the pieces of the TkrRecon stuff
    void fillTkrClusterCol( TkrRecon* recon, Event::TkrClusterCol* clusterColTds);
    void fillTruncationInfo(TkrRecon* recon, Event::TkrTruncationInfo* truncationInfoTds);
    void fillFitTracks( TkrRecon* recon, Event::TkrTrackCol* tracksTds);
    void fillTkrVecPoints(TkrRecon* recon, const Event::TkrVecPointCol* vecPointsTds);
    void fillTkrVecPointsLinks(TkrRecon* recon, const Event::TkrVecPointsLinkCol* vecPointsLinksTds);
    void fillTkrEventParams(TkrRecon* recon, const Event::TkrEventParams* eventParamsTds);
    void fillTkrFilterParams(TkrRecon* recon, const Event::TkrFilterParamsCol* filterParamsTds);
    void fillTkrTrees(TkrRecon* recon, const Event::TkrTreeCol* treeTds);
    void fillTkrTreeCompressed(TkrRecon* recon, const Event::TkrTreeCol* treeTds);
    void fillVertices( TkrRecon* recon, Event::TkrVertexCol*   verticesTds, Event::TkrTrackCol* tracksTds);
    
    // Special recursive method for filling TkrVecNodes (since they are THE Tree!)
    TkrVecNode*           convertTkrVecNode(TkrRecon* recon, const Event::TkrVecNode* vecNode);
    TkrVecNodeCompressed* convertTkrVecNodeCompressed(TkrRecon* recon, const Event::TkrVecNode* vecNode);

    // Note that the TkrFilterParams not only hold Hough Filter results (in the TkrFilterParamsCol)
    // but are also the container of the Tree axis parameters. Hence we need a separate conversion
    TkrFilterParams* convertTkrFilterParams(const Event::TkrFilterParams* filterParams);

    TkrTrackHit*   convertTkrTrackHit(const Event::TkrTrackHit* trackHitTds);
    TkrTrackParams convertTkrTrackParams(const Event::TkrTrackParams& trackParamsTds);

    /// This will be for writing the objects that relate subsystems
    StatusCode writeTkrCalAcdRelations();
    /// Use this to deal with the TreeClusterRelations if they exist
    void fillTreeClusterRelations(const Event::TreeClusterRelationCol* treeClusterRelationCol);
    
    /// Retrieves the CAL reconstruction data from the TDS and fills the CalRecon
    /// ROOT object
    StatusCode writeCalRecon();
    
    /// These are the methods specific to filling the pieces of the CalRecon stuff
    void fillCalCluster(CalRecon *calRec, Event::CalClusterCol* clusterColTds);   
    void fillCalCluster(CalRecon* calRec, Event::CalClusterMap* clusterMapTds);
    void fillCalXtalRec(CalRecon *calRec, Event::CalXtalRecCol* xtalColTds); 
    void fillCalMipTrack(CalRecon *calRec, Event::CalMipTrackCol* calMipTrackColTds); 
    void fillCalEventEnergy(CalRecon *calRec, Event::CalEventEnergyCol* calEventEnergyCol);
    void fillCalEventEnergyMap(CalRecon* calRec, Event::CalEventEnergyMap* calEventEnergyMap);
    
    //CL: 08/22/06:
    void fillGcrXtal(CalRecon *calRec, Event::GcrXtalCol* gcrXtalColTds); 
    //void fillGcrSelectedXtal(CalRecon *calRec, Event::GcrSelectedXtalsCol* gcrSelectedXtalColTds); 
    void fillGcrTrack(CalRecon *calRec, Event::GcrTrack* gcrTrackTds); 

    StatusCode writeAcdRecon();

    StatusCode writeAdfRecon();
    
    /// Calls TTree::Fill for each event and clears m_mcEvt
    void writeEvent();
    
    /// Performs the final write to the ROOT file and closes
    void close();
    
    /// ROOT file pointer
    TFile *m_reconFile;
    /// ROOT tree pointer
    TTree *m_reconTree;
    /// Top-level Monte Carlo ROOT object
    ReconEvent *m_reconEvt;
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

    /// Keep track of relation between TDS objs and ROOT counterparts
    commonData m_common;
    IRootIoSvc* m_rootIoSvc;

    // ADDED FOR THE FILE HEADERS DEMO
    IFhTool * m_headersTool ;

    // writeTruncationInfo flag, just to be safe
    bool m_writeTruncationInfo;
    /// Should we write the Tree Based collection to PDS?
    bool m_writeTreeBasedInfo;
    /// If above is true, should we write out the full detail or summary?
    bool m_writeFullTreeOutput;
};


//static const AlgFactory<reconRootWriterAlg>  Factory;
//const IAlgFactory& reconRootWriterAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(reconRootWriterAlg);

reconRootWriterAlg::reconRootWriterAlg(const std::string& name, 
                                       ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator)
{
    // Input parameters available to be set via the jobOptions file
    declareProperty("reconRootFile",m_fileName="recon.root");
    declareProperty("splitMode", m_splitMode=1);
    declareProperty("bufferSize", m_bufSize=64000);
    // ROOT default compression
    declareProperty("compressionLevel", m_compressionLevel=1);
    declareProperty("treeName", m_treeName="Recon");

    declareProperty("writeTruncationInfo", m_writeTruncationInfo=true);

    declareProperty("writeTreeBasedInfo",  m_writeTreeBasedInfo=true);
    declareProperty("WriteFullTreeOutput", m_writeFullTreeOutput=true);
    
    // ADDED FOR THE FILE HEADERS DEMO
    m_headersTool = 0 ;
}

StatusCode reconRootWriterAlg::initialize()
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
        // RootIoSvc is now required for reading/writing, cannot continue without it
        return StatusCode::FAILURE;
    } 

    m_reconTree = m_rootIoSvc->prepareRootOutput("recon", m_fileName, m_treeName, 
        m_compressionLevel, "GLAST Reconstruction Data");

    m_reconEvt = new ReconEvent();
    m_rootIoSvc->setupBranch("recon", "ReconEvent", "ReconEvent", &m_reconEvt, m_bufSize, m_splitMode);
    m_common.m_reconEvt = m_reconEvt;
    
    return sc;
    
}

StatusCode reconRootWriterAlg::execute()
{
    // Purpose and Method:  Called once per event.  This method calls
    //   the appropriate methods to read data from the TDS and write data
    //   to the ROOT file.
    
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    if(!m_reconTree->GetCurrentFile()->IsOpen()) {
        log << MSG::ERROR << "ROOT file " << m_fileName 
            << " could not be opened for writing." << endreq;
        return StatusCode::FAILURE;
    }
    
    if (m_reconEvt) m_reconEvt->Clear();

    sc = writeReconEvent();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to write ReconEvent" << endreq;
        return sc;
    }
    
    sc = writeTkrRecon();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to write Tkr Recon data" << endreq;
        return sc;
    }
    
    sc = writeCalRecon();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to write Cal Recon Data" << endreq;
        return sc;
    }
    
    sc = writeAcdRecon();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to write Acd Recon Data" << endreq;
        return sc;
    }

    sc = writeTkrCalAcdRelations();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to write relation Recon Data" << endreq;
        return sc;
    }

    sc = writeAdfRecon();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to write Adf Recon Data" << endreq;
        return sc;
    }

    
    writeEvent();
    return sc;
}


StatusCode reconRootWriterAlg::writeReconEvent() {
    // Purpose and Method:  Retrieve the Event object from the TDS and set the
    //    event and run numbers in the ReconEvent ROOT object
    
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    // Retrieve the Event data for this event
    SmartDataPtr<Event::EventHeader> evtTds(eventSvc(), EventModel::EventHeader);
    
    if (!evtTds) return sc;
    
    UInt_t evtId = evtTds->event();
    UInt_t runId = evtTds->run();
    
    log << MSG::DEBUG;
    if( log.isActive())evtTds->fillStream(log.stream());
    log << endreq;
    
    m_reconEvt->initialize(evtId, runId, new TkrRecon, new CalRecon, new AcdRecon, new AcdReconV2);

    // For simulated data - this may not exist on the TDS and that is ok
    // no need to fail for that
    SmartDataPtr<LdfEvent::EventSummaryData> summaryTds(eventSvc(), "/Event/EventSummary");
    if (summaryTds) m_reconEvt->initEventFlags(summaryTds->eventFlags());
    
    m_reconEvt->initGleamEventFlags(evtTds->gleamEventFlags());

    return sc;
}

StatusCode reconRootWriterAlg::writeTkrRecon() {
    // Purpose and Method:  Retrieve the Tkr Recon data from the TDS
    
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    // Get pointer to TkrRecon part of ReconEvent
    TkrRecon* recon = m_reconEvt->getTkrRecon();
    if (!recon) return StatusCode::FAILURE;
    recon->initialize();

    SmartDataPtr<Event::TkrClusterCol> clusterColTds(eventSvc(), EventModel::TkrRecon::TkrClusterCol);
    if (clusterColTds) fillTkrClusterCol(recon, clusterColTds);

    if(m_writeTruncationInfo) {
        SmartDataPtr<Event::TkrTruncationInfo> truncationInfoTds(eventSvc(), EventModel::TkrRecon::TkrTruncationInfo);
        if (truncationInfoTds) {
            if(truncationInfoTds->isTruncated() ) fillTruncationInfo(recon, truncationInfoTds);
        }
    }

    // Retrieve the information on fit tracks
    SmartDataPtr<Event::TkrTrackCol> tracksTds(eventSvc(), EventModel::TkrRecon::TkrTrackCol);
    
    // Fill the fit tracks
    if (tracksTds) fillFitTracks(recon, tracksTds);
    
    // Retrieve the information on vertices
    SmartDataPtr<Event::TkrVertexCol> verticesTds(eventSvc(), EventModel::TkrRecon::TkrVertexCol);
    
    // Fill the vertices
    if (verticesTds && tracksTds) fillVertices(recon, verticesTds, tracksTds);

    // now add the cosmic-ray tracks
    SmartDataPtr<Event::TkrTrackCol> crTracksTds(eventSvc(), EventModel::TkrRecon::TkrCRTrackCol);
    if (crTracksTds) fillFitTracks(recon, crTracksTds);

    // If there are TkrEventParams be sure to add them
    SmartDataPtr<Event::TkrEventParams> eventParamsTds(eventSvc(), EventModel::TkrRecon::TkrEventParams);
    if (eventParamsTds) fillTkrEventParams(recon, eventParamsTds);

    // If there are TkrFilterParams be sure to add them
    SmartDataPtr<Event::TkrFilterParamsCol> filterParamsColTds(eventSvc(), EventModel::TkrRecon::TkrFilterParamsCol);
    if (filterParamsColTds) fillTkrFilterParams(recon, filterParamsColTds);

    // If requested (and its there), write out the full Tree Based info tree
    if (m_writeTreeBasedInfo)
    {
        // Write out the full output of all the Tree Based Tracking
        if (m_writeFullTreeOutput)
        {
            // Recover the collection of TkrVecPoints
            SmartDataPtr<Event::TkrVecPointCol> vecPointCol(eventSvc(), EventModel::TkrRecon::TkrVecPointCol);

            // Fill away
            if (vecPointCol) fillTkrVecPoints(recon, vecPointCol);

            // Recover the collection of TkrVecPointsLinks
            SmartDataPtr<Event::TkrVecPointsLinkCol> vecPointsLinkCol(eventSvc(), EventModel::TkrRecon::TkrVecPointsLinkCol);

            // Fill away
            if (vecPointsLinkCol) fillTkrVecPointsLinks(recon, vecPointsLinkCol);

            // Next up is to fill the TkrVecNodes... note that these are filled while we also handle the TkrTrees
            // so here we simply look up the TkrTrees and proceed that way
            SmartDataPtr<Event::TkrTreeCol> treeColTds(eventSvc(), EventModel::TkrRecon::TkrTreeCol);

            if (treeColTds) fillTkrTrees(recon, treeColTds);
        }
        // Otherwise we are trying to write the compressed version
        else
        {
            // Find the tree collection...
            SmartDataPtr<Event::TkrTreeCol> treeColTds(eventSvc(), EventModel::TkrRecon::TkrTreeCol);

            if (treeColTds) fillTkrTreeCompressed(recon, treeColTds);
        }
    }

    // Future Diagnostics will plug in here
    //Event::TkrDiagnostics* diagnosticsTds = SmartDataPtr<Event::TkrDiagnostics>(eventSvc(), EventModel::TkrRecon::TkrDiagnostics);
    //
    //TkrDiagnostics* diagnosticsRoot = new TkrDiagnostics();
    //
    //if (diagnosticsTds != 0)
    //{
    //    diagnosticsRoot->initializeInfo(diagnosticsTds->getNumClusters(),
    //                                    diagnosticsTds->getNumVecPoints(),
    //                                    diagnosticsTds->getNumVecLinks(),
    //                                    diagnosticsTds->getnLinksNonZeroLayers(),
    //                                    diagnosticsTds->getAveNumLinksLayer(),
    //                                    diagnosticsTds->getNumLinkCombinations(),
    //                                    diagnosticsTds->getNumTrackElements(),
    //                                    diagnosticsTds->getNumTkrTracks());
    //}
    //else
    //{
    //    diagnosticsRoot->Clear();
    //}
    //
    //recon->addDiagnostics(diagnosticsRoot);

    return sc;
}

void reconRootWriterAlg::fillTkrClusterCol(TkrRecon* recon, Event::TkrClusterCol* clusterColTds) {
    
    // Purpose and Method:  Reads collection of clusters from TDS and copies their
    //   contents to a ROOT version.
    
    //int numClusters = clusterColTds->size();

    Event::TkrClusterColConItr clusterIter;
    for (clusterIter = clusterColTds->begin(); clusterIter != clusterColTds->end(); clusterIter++)
    {
        Event::TkrCluster*    clusterTds = *clusterIter;
        idents::TkrId         clusTdsId  = clusterTds->getTkrId();
        commonRootData::TkrId clusId(clusTdsId.getTowerX(),clusTdsId.getTowerY(),clusTdsId.getTray(),
                                     clusTdsId.getBotTop(), clusTdsId.getView());
        Point posTds = clusterTds->position();
        TVector3 posRoot(posTds.x(), posTds.y(), posTds.z());
       
        TkrCluster *clusterRoot = new TkrCluster(clusId, 
                                                 clusterTds->firstStrip(), 
                                                 clusterTds->lastStrip(), 
                                                 posRoot, 
                                                 (int)clusterTds->getRawToT(),
                                                 clusterTds->getMips(), 
                                                 clusterTds->getStatusWord(),
                                                 clusterTds->getNBad()
                                                 );
        
        recon->addCluster(clusterRoot);

        TRef ref = clusterRoot;
        m_common.m_tkrClusterMap[clusterTds] = ref;
    }
}

void reconRootWriterAlg::fillTruncationInfo(TkrRecon* recon, Event::TkrTruncationInfo* truncationInfoTds)
{
    //Purpose and Method : fill ROOT persistent version of the truncation info with the TDS class.
    if(truncationInfoTds==0) return;
    Event::TkrTruncationInfo::TkrTruncationMap* truncMap = truncationInfoTds->getTruncationMap();
    if(truncMap==0) return;
    Event::TkrTruncationInfo::TkrTruncationMap::const_iterator iter = truncMap->begin();
    for(; iter!=truncMap->end(); ++iter) {
        TkrTruncationData* truncationDataRoot = new TkrTruncationData();
        RootPersistence::convert(*iter, *truncationDataRoot);
        recon->addTruncationData(truncationDataRoot);

    }
}
//*/

void reconRootWriterAlg::fillFitTracks(TkrRecon* recon, Event::TkrTrackCol* tracksTds)
{
    // Purpose and Method:  This creates root tracks from tds tracks 
    //                      and adds them to the list kept in TkrRecon

    // Iterate over the tracks in the TDS
    // Purpose and Method:  This converts Event::TkrKalFitTrack objects into root TkrKalFitTrack objects
    for(Event::TkrTrackColConPtr trkPtr = tracksTds->begin(); trkPtr != tracksTds->end(); trkPtr++)
    {
        Event::TkrTrack* trackTds  = *trkPtr;
        TkrTrack*        trackRoot = new TkrTrack();

        // Convert these to TVector3's 
        TVector3 position(trackTds->getInitialPosition().x(),
                          trackTds->getInitialPosition().y(),
                          trackTds->getInitialPosition().z());
        TVector3 direction(trackTds->getInitialDirection().x(),
                           trackTds->getInitialDirection().y(),
                           trackTds->getInitialDirection().z());

        // Begin with filling the track member variables
        trackRoot->setStatusBit(trackTds->getStatusBits());
        trackRoot->setInitialEnergy(trackTds->getInitialEnergy());
        trackRoot->setInitialPosition(position);
        trackRoot->setInitialDirection(direction);
        trackRoot->setChiSquareFilter(trackTds->getChiSquareFilter());
        trackRoot->setChiSquareSmooth(trackTds->getChiSquareSmooth());
        trackRoot->setNDegreesOfFreedom(trackTds->getNDegreesOfFreedom());
        trackRoot->setScatter(trackTds->getScatter());
        trackRoot->setQuality(trackTds->getQuality());
        trackRoot->setKalEnergy(trackTds->getKalEnergy());
        trackRoot->setKalEnergyError(trackTds->getKalEnergyError());
        trackRoot->setKalThetaMS(trackTds->getKalThetaMS());
        trackRoot->setNumXGaps(trackTds->getNumXGaps());
        trackRoot->setNumYGaps(trackTds->getNumYGaps());
        trackRoot->setNumXFirstGaps(trackTds->getNumXFirstGaps());
        trackRoot->setNumYFirstGaps(trackTds->getNumYFirstGaps());
        trackRoot->setNumSegmentPoints(trackTds->getNumSegmentPoints());
        trackRoot->setChiSqSegment(trackTds->chiSquareSegment());
        trackRoot->setNumXHits(trackTds->getNumXHits());
        trackRoot->setNumYHits(trackTds->getNumYHits());
        trackRoot->setTkrCalRadLen(trackTds->getTkrCalRadlen());
        trackRoot->setRangeEnergy(trackTds->getRangeEnergy());

        // Now loop over the hit planes and fill that information
        Event::TkrTrackHitVecConItr trkHitTdsItr = trackTds->begin();
        while(trkHitTdsItr != trackTds->end())
        {
            const Event::TkrTrackHit*  trackHitTds = *trkHitTdsItr++;
            TkrTrackHit*               trackHit    = convertTkrTrackHit(trackHitTds);
        
            trackRoot->Add(trackHit);
        }

        // Ok, I think this makes TkrTrack the owner of the TkrTrackHits... I think...
        trackRoot->SetOwner();

        // Keep relation between Event and Root candidate tracks
        TRef ref = trackRoot;
        m_common.m_tkrTrackMap[trackTds] = ref;

        // Ok, now add the track to the right collection!
        if((trackRoot->getStatusBits()&Event::TkrTrack::COSMICRAY)==0) {
            recon->addTrack(trackRoot);
        } else {
            recon->addCRTrack(trackRoot);
        }
    }
    
    return;
}

void reconRootWriterAlg::fillTkrVecPoints(TkrRecon* recon, const Event::TkrVecPointCol* vecPointsTds)
{
    // Purpose and Method:  This takes the TkrVecPoint collection from the TDS,
    //                      creates root equivalents and stores them in the world of root

    // Iterate over the TkrVecPoint collection in the TDS
    // Purpose and Method:  This converts Event::TkrKalFitTrack objects into root TkrKalFitTrack objects
    for(Event::TkrVecPointCol::const_iterator vecPtItr = vecPointsTds->begin(); vecPtItr != vecPointsTds->end(); vecPtItr++)
    {
        // Recover pointer to the TkrVecPoint object in TDS
        Event::TkrVecPoint* vecPointTds  = *vecPtItr;

        // Create a shiny new, unused, in the orignal packaging, root version of TkrVecPoint
        TkrVecPoint*        vecPointRoot = new TkrVecPoint();

        // Recover the root equivalent pointers to the clusters which make up this point
        TRef        tkrXClusterRef  = m_common.m_tkrClusterMap[vecPointTds->getXCluster()];
        TkrCluster* tkrXClusterRoot = (TkrCluster*) tkrXClusterRef.GetObject();
        TRef        tkrYClusterRef  = m_common.m_tkrClusterMap[vecPointTds->getYCluster()];
        TkrCluster* tkrYClusterRoot = (TkrCluster*) tkrYClusterRef.GetObject();

        // Initialize the root TkrVecPoint
        vecPointRoot->initializeInfo(vecPointTds->getLayer(),
                                     vecPointTds->getStatusWord(),
                                     tkrXClusterRoot,
                                     tkrYClusterRoot);

        // Keep relation between Event and Root TkrVecPoints
        TRef ref = vecPointRoot;
        m_common.m_tkrVecPointMap[vecPointTds] = ref;

        // Ok, now add the TkrVecPoint to the right collection!
        recon->addTkrVecPoint(vecPointRoot);
    }
    
    return;
}

void reconRootWriterAlg::fillTkrVecPointsLinks(TkrRecon* recon, const Event::TkrVecPointsLinkCol* vecPointsLinkColTds)
{
    // Purpose and Method:  This takes the TkrVecPoint collection from the TDS,
    //                      creates root equivalents and stores them in the world of root

    // Iterate over the TkrVecPoint collection in the TDS
    // Purpose and Method:  This converts Event::TkrKalFitTrack objects into root TkrKalFitTrack objects
    for(Event::TkrVecPointsLinkCol::const_iterator vecPtItr  = vecPointsLinkColTds->begin(); 
                                                   vecPtItr != vecPointsLinkColTds->end(); 
                                                   vecPtItr++)
    {
        // Recover pointer to the TkrVecPoint object in TDS
        Event::TkrVecPointsLink* vecPointsLinkTds  = *vecPtItr;

        // Create a shiny new, unused, in the orignal packaging, root version of TkrVecPoint
        TkrVecPointsLink*        vecPointsLinkRoot = new TkrVecPointsLink();

        // Recover the root equivalent pointers to the clusters which make up this point
        TRef         topVecPointRef  = m_common.m_tkrVecPointMap[vecPointsLinkTds->getFirstVecPoint()];
        TkrVecPoint* topVecPointRoot = (TkrVecPoint*) topVecPointRef.GetObject();
        TRef         botVecPointRef  = m_common.m_tkrVecPointMap[vecPointsLinkTds->getSecondVecPoint()];
        TkrVecPoint* botVecPointRoot = (TkrVecPoint*) botVecPointRef.GetObject();

        // Initialize the root TkrVecPoint
        vecPointsLinkRoot->initializeInfo(topVecPointRoot, 
                                          botVecPointRoot, 
                                          vecPointsLinkTds->getStatusBits(),
                                          vecPointsLinkTds->getMaxScatAngle());

        // Keep relation between Event and Root TkrVecPoints
        TRef ref = vecPointsLinkRoot;
        m_common.m_tkrVecPointsLinkMap[vecPointsLinkTds] = ref;

        // Ok, now add the TkrVecPoint to the right collection!
        recon->addTkrVecPointsLink(vecPointsLinkRoot);
    }
    
    return;
}

TkrFilterParams* reconRootWriterAlg::convertTkrFilterParams(const Event::TkrFilterParams* filterParamsTds)
{
    TkrFilterParams* filterParamsRoot = 0;

    if (filterParamsTds)
    {
        // Create the root equivalent object
        filterParamsRoot = new TkrFilterParams();

        // Convert the Point and Vector to TVector3 objects
        TVector3 eventPosition(filterParamsTds->getEventPosition().x(),
                               filterParamsTds->getEventPosition().y(),
                               filterParamsTds->getEventPosition().z());
        TVector3 eventAxis(filterParamsTds->getEventAxis().x(),
                           filterParamsTds->getEventAxis().y(),
                           filterParamsTds->getEventAxis().z());

        // Initialize the root TkrVecPoint
        filterParamsRoot->initializeInfo(eventPosition, 
                                         eventAxis,
                                         filterParamsTds->getStatusBits(),
                                         filterParamsTds->getEventEnergy(),
                                         filterParamsTds->getNumBiLayers(),
                                         filterParamsTds->getNumIterations(),
                                         filterParamsTds->getNumHitsTotal(),
                                         filterParamsTds->getNumDropped(),
                                         filterParamsTds->getChiSquare(),
                                         filterParamsTds->getAverageDistance(),
                                         filterParamsTds->getTransRms(), 
                                         filterParamsTds->getLongRms(),
                                         filterParamsTds->getLongRmsAsym() );
    }

    return filterParamsRoot;
}

void reconRootWriterAlg::fillTkrFilterParams(TkrRecon* recon, const Event::TkrFilterParamsCol* filterParamsColTds)
{
    // Purpose and Method:  This takes the TkrFilterParams collection from the TDS,
    //                      creates root equivalents and stores them in the world of root

    // Iterate over the input TkrFilterParams collection in the TDS
    for(Event::TkrFilterParamsCol::const_iterator filterItr  = filterParamsColTds->begin(); 
                                                  filterItr != filterParamsColTds->end(); 
                                                  filterItr++)
    {
        // Recover pointer to the TkrVecPoint object in TDS
        Event::TkrFilterParams* filterParamsTds  = *filterItr;

        // Use the "converter" to get the object
        TkrFilterParams*        filterParamsRoot = convertTkrFilterParams(filterParamsTds);

        // Keep relation between Event and Root TkrVecPoints
        TRef ref = filterParamsRoot;
        m_common.m_tkrFilterParamsMap[filterParamsTds] = ref;

        // Ok, now add the TkrVecPoint to the right collection!
        recon->addTkrFilterParams(filterParamsRoot);
    }

     return;
}

void reconRootWriterAlg::fillTkrEventParams(TkrRecon* recon, const Event::TkrEventParams* eventParamsTds)
{
    // Purpose and Method:  This takes the TkrFilterParams collection from the TDS,
    //                      creates root equivalents and stores them in the world of root

    // This test should not fail...
    if (eventParamsTds)
    {
        // Obtain a shiny new TkrEventParams root object from the universe
        TkrEventParams* eventParamsRoot = new TkrEventParams;

        // Convert the Point and Vector to TVector3 objects
        TVector3 eventPosition(eventParamsTds->getEventPosition().x(),
                               eventParamsTds->getEventPosition().y(),
                               eventParamsTds->getEventPosition().z());
        TVector3 eventAxis(eventParamsTds->getEventAxis().x(),
                           eventParamsTds->getEventAxis().y(),
                           eventParamsTds->getEventAxis().z());

        // Initialize the root TkrVecPoint
        eventParamsRoot->initializeInfo(eventPosition, 
                                        eventAxis,
                                        eventParamsTds->getStatusBits(),
                                        eventParamsTds->getEventEnergy(),
                                        eventParamsTds->getNumBiLayers(),
                                        eventParamsTds->getNumIterations(),
                                        eventParamsTds->getNumHitsTotal(),
                                        eventParamsTds->getNumDropped(),
                                        eventParamsTds->getChiSquare(),
                                        eventParamsTds->getTransRms(), 
                                        eventParamsTds->getLongRmsAve() );

        // Ok, now add the TkrVecPoint to the right collection!
        recon->addTkrEventParams(eventParamsRoot);
    }

    return;
}

TkrVecNode* reconRootWriterAlg::convertTkrVecNode(TkrRecon* recon, const Event::TkrVecNode* vecNodeTds)
{
    TkrVecNode* vecNodeRoot = 0;

    // If an input node then something to do! 
    if (vecNodeTds)
    {
        // Create the root equivalent node to match the TDS version
        vecNodeRoot = new TkrVecNode();

        // Find the reference to the parent node, if it exists!
        TkrVecNode* parentNodeRoot = 0;

        if (vecNodeTds->getParentNode())
        {
            TRef parentNodeRef  = m_common.m_tkrVecNodeMap[vecNodeTds->getParentNode()];

            parentNodeRoot = (TkrVecNode*) parentNodeRef.GetObject();
        }

        // Same song and dance for the associated link
        TkrVecPointsLink* associatedLinkRoot = 0;

        if (vecNodeTds->getAssociatedLink())
        {
            TRef associatedLinkRef = m_common.m_tkrVecPointsLinkMap[vecNodeTds->getAssociatedLink()];

            associatedLinkRoot = (TkrVecPointsLink*) associatedLinkRef.GetObject();
        }

        // Now initialize the root TkrVecNode with the values in the TDS version
        vecNodeRoot->initializeInfo(parentNodeRoot,
                                    associatedLinkRoot,
                                    vecNodeTds->getStatusBits(),
                                    vecNodeTds->getRmsAngleSum(),
                                    vecNodeTds->getNumAnglesInSum(),
                                    vecNodeTds->getNumLeaves(),
                                    vecNodeTds->getNumBranches(),
                                    vecNodeTds->getDepth(),
                                    vecNodeTds->getBestNumBiLayers(),
                                    vecNodeTds->getBestRmsAngle());

        // Add this node to the map of known nodes
        TRef ref = vecNodeRoot;
        m_common.m_tkrVecNodeMap[vecNodeTds] = ref;

        // And add to the recon object
        recon->addTkrVecNode(vecNodeRoot);

        // Loop over daughter nodes to set them
        for(Event::TkrVecNodeSet::const_iterator nodeItr = vecNodeTds->begin(); nodeItr != vecNodeTds->end(); nodeItr++)
        {
            // Recover the daughter node
            Event::TkrVecNode* daughterNodeTds = *nodeItr;

            // Use this to get a root daughter node
            TkrVecNode* daughterNodeRoot = convertTkrVecNode(recon, daughterNodeTds);

            // If we return with something then add to our daughter collection
            if (daughterNodeRoot) vecNodeRoot->Add(daughterNodeRoot);
        }
    }

    return vecNodeRoot;
}

TkrVecNodeCompressed* reconRootWriterAlg::convertTkrVecNodeCompressed(TkrRecon* recon, const Event::TkrVecNode* vecNodeTds)
{
    TkrVecNodeCompressed* vecNodeRoot = 0;

    // If an input node then something to do! 
    if (vecNodeTds)
    {
        // Create the root equivalent node to match the TDS version
        vecNodeRoot = new TkrVecNodeCompressed();

        // Find the reference to the parent node, if it exists!
        TkrVecNodeCompressed* parentNodeRoot = 0;

        if (vecNodeTds->getParentNode())
        {
            TRef parentNodeRef  = m_common.m_tkrVecNodeMap[vecNodeTds->getParentNode()];

            parentNodeRoot = (TkrVecNodeCompressed*) parentNodeRef.GetObject();
        }

        // Things are a bit different for the associated link. In compressed mode we are going to
        // store the pointers to the clusters associated to the 3D points that make the endpoints of the link.
        // And various status words...
        TkrCluster* topPointXClusterRoot = 0;
        TkrCluster* topPointYClusterRoot = 0;
        UInt_t      topPointStatus       = 0;
        TkrCluster* botPointXClusterRoot = 0;
        TkrCluster* botPointYClusterRoot = 0;
        UInt_t      botPointStatus       = 0;
        UInt_t      linkStatus           = 0;

        if (vecNodeTds->getAssociatedLink())
        {
            const Event::TkrVecPointsLink* linkTds             = vecNodeTds->getAssociatedLink();
            const Event::TkrVecPoint*      topPointTds         = linkTds->getFirstVecPoint();
            const Event::TkrCluster*       topPointXClusterTds = topPointTds->getXCluster();
            const Event::TkrCluster*       topPointYClusterTds = topPointTds->getYCluster();
            const Event::TkrVecPoint*      botPointTds         = linkTds->getSecondVecPoint();
            const Event::TkrCluster*       botPointXClusterTds = botPointTds->getXCluster();
            const Event::TkrCluster*       botPointYClusterTds = botPointTds->getYCluster();

            // recover the status words
            topPointStatus = const_cast<Event::TkrVecPoint*>(topPointTds)->getStatusWord();
            botPointStatus = const_cast<Event::TkrVecPoint*>(botPointTds)->getStatusWord();
            linkStatus     = linkTds->getStatusBits();

            // Now recover the pointers to the root versions of the clusters
            TRef topXClusterRef = m_common.m_tkrClusterMap[topPointXClusterTds];
            TRef topYClusterRef = m_common.m_tkrClusterMap[topPointYClusterTds];
            TRef botXClusterRef = m_common.m_tkrClusterMap[botPointXClusterTds];
            TRef botYClusterRef = m_common.m_tkrClusterMap[botPointYClusterTds];

            topPointXClusterRoot = (TkrCluster*) topXClusterRef.GetObject();
            topPointYClusterRoot = (TkrCluster*) topYClusterRef.GetObject();
            botPointXClusterRoot = (TkrCluster*) botXClusterRef.GetObject();
            botPointYClusterRoot = (TkrCluster*) botYClusterRef.GetObject();
        }

        // Now initialize the root TkrVecNode with the values in the TDS version
        vecNodeRoot->initializeInfo(parentNodeRoot,
                                    vecNodeTds->getStatusBits(),
                                    vecNodeTds->getRmsAngleSum(),
                                    vecNodeTds->getNumAnglesInSum(),
                                    vecNodeTds->getNumLeaves(),
                                    vecNodeTds->getNumBranches(),
                                    vecNodeTds->getDepth(),
                                    vecNodeTds->getBestNumBiLayers(),
                                    vecNodeTds->getBestRmsAngle());

        // Initialize the cluster/status word info
        vecNodeRoot->initializeLinkInfo(topPointXClusterRoot,
                                        topPointYClusterRoot,
                                        topPointStatus,
                                        botPointXClusterRoot,
                                        botPointYClusterRoot,
                                        botPointStatus,
                                        linkStatus);

        // Add this node to the map of known nodes
        TRef ref = vecNodeRoot;
        m_common.m_tkrVecNodeMap[vecNodeTds] = ref;

        // And add to the recon object
        recon->addTkrVecNodeCompressed(vecNodeRoot);

        // Loop over daughter nodes to set them
        for(Event::TkrVecNodeSet::const_iterator nodeItr = vecNodeTds->begin(); nodeItr != vecNodeTds->end(); nodeItr++)
        {
            // Recover the daughter node
            Event::TkrVecNode* daughterNodeTds = *nodeItr;

            // Use this to get a root daughter node
            TkrVecNodeCompressed* daughterNodeRoot = convertTkrVecNodeCompressed(recon, daughterNodeTds);

            // If we return with something then add to our daughter collection
            if (daughterNodeRoot) vecNodeRoot->Add(daughterNodeRoot);
        }
    }

    return vecNodeRoot;
}

void reconRootWriterAlg::fillTkrTrees(TkrRecon* recon, const Event::TkrTreeCol* treeColTds)
{
    // Purpose and Method:  This takes the TkrTree collection from the TDS,
    //                      creates root equivalents and stores them in the world of root

    // Iterate over the input TkrTree collection in the TDS
    for(Event::TkrTreeCol::const_iterator treeItr  = treeColTds->begin(); 
                                          treeItr != treeColTds->end(); 
                                          treeItr++)
    {
        // Recover pointer to the TkrVecPoint object in TDS
        Event::TkrTree* treeTds  = *treeItr;

        // Create the root equivalent object
        TkrTree*        treeRoot = new TkrTree();

        // Before anything else, convert all the TkrVecNodes comprising this Tree
        TkrVecNode* headNodeRoot = convertTkrVecNode(recon, treeTds->getHeadNode());

        // Once the above is done, then we can look up references to the best/second branch nodes
        TkrVecNode* bestLeafNodeRoot   = 0;
        TkrVecNode* secondLeafNodeRoot = 0;

        if (treeTds->getBestLeaf())
        {
            TRef bestLeafNodeRef = m_common.m_tkrVecNodeMap[treeTds->getBestLeaf()];

            bestLeafNodeRoot = (TkrVecNode*) bestLeafNodeRef.GetObject();
        }

        if (treeTds->getSecondLeaf())
        {
            TRef secondLeafNodeRef = m_common.m_tkrVecNodeMap[treeTds->getSecondLeaf()];

            secondLeafNodeRoot = (TkrVecNode*) secondLeafNodeRef.GetObject();
        }

        // Don't forget to do the same for the tracks!
        TRef trackRootRef = m_common.m_tkrTrackMap[treeTds->getBestTrack()];

        TkrTrack* bestTrack   = (TkrTrack*) trackRootRef.GetObject();
        TkrTrack* secondTrack = 0;

        if (treeTds->size() > 1)
        {
            trackRootRef = m_common.m_tkrTrackMap[treeTds->back()];
            secondTrack  = (TkrTrack*) trackRootRef.GetObject();
        }

        // Get the "filter params" (Tree Axis information)
        TkrFilterParams* filterParamsRoot = convertTkrFilterParams(treeTds->getAxisParams());

        // Initialize the root TkrVecPoint
        treeRoot->initializeInfo(headNodeRoot, 
                                 bestLeafNodeRoot, 
                                 secondLeafNodeRoot, 
                                 filterParamsRoot,
                                 treeTds->getBestBranchAngleToAxis(),
                                 treeTds->getAxisSeededAngleToAxis() );

        // Add the tracks
        treeRoot->Add(bestTrack);
        if (secondTrack) treeRoot->Add(secondTrack);

        // Keep relation between Event and Root TkrVecPoints
        TRef ref = treeRoot;
        m_common.m_tkrTreeMap[treeTds] = ref;

        // Ok, now add the TkrVecPoint to the right collection!
        recon->addTkrTree(treeRoot);
    }

    return;
}

void reconRootWriterAlg::fillTkrTreeCompressed(TkrRecon* recon, const Event::TkrTreeCol* treeColTds)
{
    // Purpose and Method:  This takes the TkrTree collection from the TDS,
    //                      creates root equivalents and stores them in the world of root

    // Iterate over the input TkrTree collection in the TDS
    for(Event::TkrTreeCol::const_iterator treeItr  = treeColTds->begin(); 
                                          treeItr != treeColTds->end(); 
                                          treeItr++)
    {
        // Recover pointer to the TkrVecPoint object in TDS
        Event::TkrTree* treeTds  = *treeItr;

        // Create the root equivalent object
        TkrTreeCompressed* treeRoot = new TkrTreeCompressed();

        // Before anything else, convert all the TkrVecNodes comprising this Tree
        TkrVecNodeCompressed* headNodeRoot = convertTkrVecNodeCompressed(recon, treeTds->getHeadNode());

        // Once the above is done, then we can look up references to the best/second branch nodes
        TkrVecNodeCompressed* bestLeafNodeRoot   = 0;
        TkrVecNodeCompressed* secondLeafNodeRoot = 0;

        if (treeTds->getBestLeaf())
        {
            TRef bestLeafNodeRef = m_common.m_tkrVecNodeMap[treeTds->getBestLeaf()];

            bestLeafNodeRoot = (TkrVecNodeCompressed*) bestLeafNodeRef.GetObject();
        }

        if (treeTds->getSecondLeaf())
        {
            TRef secondLeafNodeRef = m_common.m_tkrVecNodeMap[treeTds->getSecondLeaf()];

            secondLeafNodeRoot = (TkrVecNodeCompressed*) secondLeafNodeRef.GetObject();
        }

        // Don't forget to do the same for the tracks!
        TRef trackRootRef = m_common.m_tkrTrackMap[treeTds->getBestTrack()];

        TkrTrack* bestTrack   = (TkrTrack*) trackRootRef.GetObject();
        TkrTrack* secondTrack = 0;

        if (treeTds->size() > 1)
        {
            trackRootRef = m_common.m_tkrTrackMap[treeTds->back()];
            secondTrack  = (TkrTrack*) trackRootRef.GetObject();
        }

        // Get the "filter params" (Tree Axis information)
        TkrFilterParams* filterParamsRoot = convertTkrFilterParams(treeTds->getAxisParams());

        // Initialize the root TkrVecPoint
        treeRoot->initializeInfo(headNodeRoot, 
                                 bestLeafNodeRoot, 
                                 secondLeafNodeRoot, 
                                 filterParamsRoot,
                                 treeTds->getBestBranchAngleToAxis(),
                                 treeTds->getAxisSeededAngleToAxis() );

        // Add the tracks
        treeRoot->Add(bestTrack);
        if (secondTrack) treeRoot->Add(secondTrack);

        // Keep relation between Event and Root TkrVecPoints
        TRef ref = treeRoot;
        m_common.m_tkrTreeMap[treeTds] = ref;

        // Ok, now add the TkrVecPoint to the right collection!
        recon->addTkrTreeCompressed(treeRoot);
    }

    return;
}

void reconRootWriterAlg::fillVertices(TkrRecon* recon, Event::TkrVertexCol* verticesTds, Event::TkrTrackCol* /*tracksTds*/)
{
    // Purpose and Method:  This creates root vertex output objects from tds vertices 
    //                      and adds them to the list kept in TkrRecon
    
    // Loop over the candidate tracks in the TDS
    int                     vtxId  = 0;
    Event::TkrVertexConPtr  vtxPtr = verticesTds->begin();
    while(vtxPtr != verticesTds->end())
    {
        Event::TkrVertex*     vtxTds = *vtxPtr++;
        TkrVertex*            vtx    = new TkrVertex();
        TVector3              pos(vtxTds->getPosition().x(), vtxTds->getPosition().y(), vtxTds->getPosition().z());
        TVector3              dir(vtxTds->getDirection().x(),vtxTds->getDirection().y(),vtxTds->getDirection().z());
        idents::TkrId         hitTdsId  = vtxTds->getTkrId();
        commonRootData::TkrId hitId(hitTdsId.getTowerX(),hitTdsId.getTowerY(),hitTdsId.getTray(),
                                    hitTdsId.getBotTop(), hitTdsId.getView());

        vtx->setStatusBit(vtxTds->getStatusBits());
        vtx->setEnergy(vtxTds->getEnergy());
        vtx->setPosition(pos);
        vtx->setDirection(dir);
        vtx->setChiSquare(vtxTds->getChiSquare());
        vtx->setQuality(vtxTds->getQuality());
        vtx->setTkr1ArcLen(vtxTds->getTkr1ArcLen());
        vtx->setTkr2ArcLen(vtxTds->getTkr2ArcLen());
        vtx->setDOCA(vtxTds->getDOCA());
        vtx->setAddedRadLen(vtxTds->getAddedRadLen());
        vtx->setTkrID(hitId);

        TkrTrackParams params = convertTkrTrackParams(vtxTds->getVertexParams());

        vtx->setParams(params);
        
        // Now add the track ids 
        // This is pretty ugly because we don't store track ids in the TDS classes
        SmartRefVector<Event::TkrTrack>::const_iterator vtxTrkIter = vtxTds->getTrackIterBegin();
        while(vtxTrkIter != vtxTds->getTrackIterEnd())
        {
            // Basically... take the pointer to the track from the vertex list and loop
            // over the track pointers in the track list looking for the match. Assign 
            // the id according to the loop variable value.
            // This ALWAYS succeeds (and I'm also selling valuable swampland in Florida!)
            SmartRef<Event::TkrTrack> trackTds  = *vtxTrkIter++;
            TRef                      trackRef  = m_common.m_tkrTrackMap[trackTds];
            TkrTrack*                 trackRoot = (TkrTrack*) trackRef.GetObject();
            
            vtx->addTrack(trackRoot);
        }
        
        // Add the candidate to the list
        recon->addVertex(vtx);
        vtxId++;
    }
}

TkrTrackHit* reconRootWriterAlg::convertTkrTrackHit(const Event::TkrTrackHit* trackHitTds)
{
    TkrTrackHit* trackHit = new TkrTrackHit();
    idents::TkrId hitTdsId  = trackHitTds->getTkrId();

    // Fill cluster id only if it exists
    if (const Event::TkrCluster* clusTds   = trackHitTds->getClusterPtr())
    {
        TRef        clusRef  = m_common.m_tkrClusterMap[clusTds];
        TkrCluster* clusRoot = (TkrCluster*) clusRef.GetObject();
        trackHit->setClusterPtr(clusRoot);
    }
    if(hitTdsId.hasTowerX()&&hitTdsId.hasTowerY()&&hitTdsId.hasTray()&&hitTdsId.hasBotTop()) {
        if(hitTdsId.hasView()) {

            commonRootData::TkrId hitId(hitTdsId.getTowerX(), hitTdsId.getTowerY(),
                hitTdsId.getTray(), hitTdsId.getBotTop(), hitTdsId.getView());
            trackHit->setTkrID(hitId);

        } else { // no view if no cluster

            commonRootData::TkrId hitId(hitTdsId.getTowerX(), hitTdsId.getTowerY(),
                hitTdsId.getTray(), hitTdsId.getBotTop());
            trackHit->setTkrID(hitId);
        }
    } else {
        commonRootData::TkrId hitId = commonRootData::TkrId();
        trackHit->setTkrID(hitId);
    }
     
    // Fill in the "easy" stuff first 
    trackHit->setStatusBit(trackHitTds->getStatusBits());
    trackHit->setZPlane(trackHitTds->getZPlane());
    trackHit->setEnergy(trackHitTds->getEnergy());
    trackHit->setRadLen(trackHitTds->getRadLen());
    trackHit->setActiveDist(trackHitTds->getActiveDist());
    trackHit->setChiSquareFilter(trackHitTds->getChiSquareFilter());
    trackHit->setChiSquareSmooth(trackHitTds->getChiSquareSmooth());

    // Now set the track params
    if (trackHitTds->validMeasuredHit())  trackHit->getTrackParams(TkrTrackHit::MEASURED) = 
                    convertTkrTrackParams(trackHitTds->getTrackParams(Event::TkrTrackHit::MEASURED));
    if (trackHitTds->validPredictedHit()) trackHit->getTrackParams(TkrTrackHit::PREDICTED) = 
                    convertTkrTrackParams(trackHitTds->getTrackParams(Event::TkrTrackHit::PREDICTED));
    if (trackHitTds->validFilteredHit())  trackHit->getTrackParams(TkrTrackHit::FILTERED) = 
                    convertTkrTrackParams(trackHitTds->getTrackParams(Event::TkrTrackHit::FILTERED));
    if (trackHitTds->validSmoothedHit())  trackHit->getTrackParams(TkrTrackHit::SMOOTHED) = 
                    convertTkrTrackParams(trackHitTds->getTrackParams(Event::TkrTrackHit::SMOOTHED));
    if (trackHitTds->validRevFit())       trackHit->getTrackParams(TkrTrackHit::REVFIT) = 
                    convertTkrTrackParams(trackHitTds->getTrackParams(Event::TkrTrackHit::REVFIT));
    if (trackHitTds->validMaterial())     trackHit->getTrackParams(TkrTrackHit::QMATERIAL) = 
                    convertTkrTrackParams(trackHitTds->getTrackParams(Event::TkrTrackHit::QMATERIAL));
    
    return trackHit;
}

TkrTrackParams reconRootWriterAlg::convertTkrTrackParams(const Event::TkrTrackParams& trkParTds)
{
    TkrTrackParams params(trkParTds.getxPosition(), trkParTds.getxSlope(), 
                          trkParTds.getyPosition(), trkParTds.getySlope(),
                          trkParTds.getxPosxPos(),  trkParTds.getxPosxSlp(), 
                          trkParTds.getxPosyPos(),  trkParTds.getxPosySlp(),
                          trkParTds.getxSlpxSlp(),  trkParTds.getxSlpyPos(), 
                          trkParTds.getxSlpySlp(),  trkParTds.getyPosyPos(), 
                          trkParTds.getyPosySlp(),  trkParTds.getySlpySlp() );

    return params;
}
    
StatusCode reconRootWriterAlg::writeTkrCalAcdRelations()
{
    // Purpose and Method:  Retrieve the subsystem relation type information from the TDS
    //                      and then call the helper routines to convert to root versions
    
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    // Retrieve the tree/cluster relation collection  
    SmartDataPtr<Event::TreeClusterRelationCol> treeClusterColTds(eventSvc(), EventModel::Recon::TreeClusterRelationCol);   
    if (treeClusterColTds) fillTreeClusterRelations(treeClusterColTds);   
    
     return StatusCode::SUCCESS;
}

void reconRootWriterAlg::fillTreeClusterRelations(const Event::TreeClusterRelationCol* treeClusterRelationCol)
{
    // The Tree Cluster Relation objects are kept track of in the TkrRecon object... probably for no good reason...
    // So we need to recover it as a first step
    TkrRecon* recon = m_reconEvt->getTkrRecon();
    if (recon)
    {
        // We should probably not assume that its been initialized...
        // If it has this call will have no affect
        recon->initialize();

        // Ok, loop through the input collection
        for(Event::TreeClusterRelationCol::const_iterator relItr  = treeClusterRelationCol->begin();
                                                          relItr != treeClusterRelationCol->end();
                                                          relItr++)
        {
            Event::TreeClusterRelation* treeClusterTds = *relItr;

            // Obtain a shiny new root version
            TreeClusterRelation* treeClusterRoot = new TreeClusterRelation();

            // We need to look up the root pointers to the Tree and Cluster in this relation
            TkrTree*    treeRoot    = 0;
            CalCluster* clusterRoot = 0;

            if (treeClusterTds->getTree())
            {
                TRef treeRef = m_common.m_tkrTreeMap[treeClusterTds->getTree()];

                treeRoot = (TkrTree*) treeRef.GetObject();
            }

            if (treeClusterTds->getCluster())
            {
                TRef clusterRef = m_common.m_calClusterMap[treeClusterTds->getCluster()];

                clusterRoot = (CalCluster*) clusterRef.GetObject();
            }

            // Initialize the root object
            treeClusterRoot->initializeInfo(treeRoot,
                                            clusterRoot,
                                            treeClusterRoot->getTreeClusDoca(),
                                            treeClusterRoot->getTreeClusCosAngle(),
                                            treeClusterRoot->getTreeClusDistAtZ(),
                                            treeClusterRoot->getClusEnergy() );

            // Ok, add to our collection
            recon->addTreeClusterRelation(treeClusterRoot);
        }
    }

    return;
}


StatusCode reconRootWriterAlg::writeCalRecon() {
    // Purpose and Method:  Retrieve the Cal Recon data from the TDS 
    //  calls the two helper methods that do the work of filling the ROOT    
    //  version of the CalRecon object for this event. 
    
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    // Get pointer to CalRecon part of ReconEvent   
    CalRecon* calRec = m_reconEvt->getCalRecon();   
    if (!calRec) return StatusCode::FAILURE;   
    calRec->initialize();   

    // The way CalClusters are handled have changed from the ancient Roman times to the 
    // modern Tuscan dominated era. In the new school scheme, the CalClusterCol container serves
    // as the object owner in the TDS, but access is meant to be through the CalClusterMap. 
    // We need to be careful here about how to output in this era
    // Start by retrieving the CalClusterCol
    SmartDataPtr<Event::CalClusterCol> clusterColTds(eventSvc(), EventModel::CalRecon::CalClusterCol);

    // Check that the collection exists, if not there is nothing to do here
    if (clusterColTds)
    {
        // Start by filling the root CalClusterCol
        fillCalCluster(calRec, clusterColTds);   
    
    // Recover the CalClusterMap which will contain references to all CalClusters
    SmartDataPtr<Event::CalClusterMap> clusterMapTds(eventSvc(),EventModel::CalRecon::CalClusterMap); 

    // In modern times the cal clusters are stored in the CalClusterMap
    if (clusterMapTds) fillCalCluster(calRec, clusterMapTds);
    // However, we might be dealing with the peak of the Roman Empire era here, so leave hooks for that
    }
    
    // Retrieve the cal xtal collection   
    SmartDataPtr<Event::CalXtalRecCol> xtalColTds(eventSvc(), EventModel::CalRecon::CalXtalRecCol);   
    if (xtalColTds) fillCalXtalRec(calRec, xtalColTds); 
    
    //SG
    // Retrieve the calMipTrack collection   
    SmartDataPtr<Event::CalMipTrackCol> calMipTrackColTds(eventSvc(), EventModel::CalRecon::CalMipTrackCol);   
    if (calMipTrackColTds) fillCalMipTrack(calRec, calMipTrackColTds);   

    // Retrieve the CalEventEnergy collection
    SmartDataPtr<Event::CalEventEnergyCol> calEventEnergyColTds(eventSvc(), EventModel::CalRecon::CalEventEnergyCol) ;
    if (calEventEnergyColTds) fillCalEventEnergy(calRec, calEventEnergyColTds);

    // Retrieve the CalEventEnergyMap
    SmartDataPtr<Event::CalEventEnergyMap> calEventEnergyMapTds(eventSvc(), EventModel::CalRecon::CalEventEnergyMap) ;
    if (calEventEnergyMapTds) fillCalEventEnergyMap(calRec, calEventEnergyMapTds);

    //C.L: 08/22/06: Retrieve GcrXtal collection 
    SmartDataPtr<Event::GcrXtalCol> gcrXtalColTds(eventSvc(), EventModel::CalRecon::GcrXtalCol);       
    if (gcrXtalColTds) 
      fillGcrXtal(calRec, gcrXtalColTds); 

    //C.L: 08/22/06: Retrieve GcrTrack  
    SmartDataPtr<Event::GcrTrack> gcrTrackTds(eventSvc(), EventModel::CalRecon::GcrTrack);       
    if (gcrTrackTds) 
      fillGcrTrack(calRec, gcrTrackTds); 


    
    return sc;
}

void reconRootWriterAlg::fillCalCluster(CalRecon* calRec, Event::CalClusterMap* clusterMapTds)
{
    // Purpose and Method: Given an input CalClusterMap, we convert each entry to the root equivalent and
    // store on the PDS

    // We start by iterating over the keys in the CalClusterMap
    for(Event::CalClusterMap::iterator clusterMapTdsItr  = clusterMapTds->begin();
                                       clusterMapTdsItr != clusterMapTds->end();
                                       clusterMapTdsItr++)
    {
        // Recover the key to this entry
        const std::string& keyTds = clusterMapTdsItr->first;

        // Turn this into a TString...
        TObjString* keyRoot = new TObjString(keyTds.c_str());

        // Recover the collection
        const Event::CalClusterVec& clusterVecTds = clusterMapTdsItr->second;

        // Create the TObjArray to hold the clusters for this key
        TObjArray* clusterVecRoot = new TObjArray();

        clusterVecRoot->Clear();

        // Now we loop over the cluster vector
        for(Event::CalClusterVec::const_iterator clusterVecTdsItr  = clusterVecTds.begin();
                                                 clusterVecTdsItr != clusterVecTds.end();
                                                 clusterVecTdsItr++)
        {
            // Recover pointer to the actual Cal Cluster
            const Event::CalCluster* clusterTds = *clusterVecTdsItr;

            // Get a pointer to a root version
            if (m_common.m_calClusterMap.find(clusterTds) != m_common.m_calClusterMap.end()) 
            {
                TRef        clusterRootRef  = m_common.m_calClusterMap[clusterTds];
                CalCluster* clusterRoot     = (CalCluster*)clusterRootRef.GetObject();

            // Store in this key's TObjArray
            clusterVecRoot->Add(clusterRoot);

            // Keep track of in our TDS/root map
            }
            else
            {
                int thiscanthappen = 0;
            }
        }

        // Signal that we are not the owner of these values (for cleanup)
        clusterVecRoot->SetOwner(kFALSE);

        // Ok, now store this in our calRecon object
        calRec->getCalClusterMap()->Add(keyRoot, clusterVecRoot);
    }

    return;
}

void reconRootWriterAlg::fillCalEventEnergyMap(CalRecon* calRec, Event::CalEventEnergyMap* energyMapTds)
{
    // Purpose and Method: Given an input CalEventEnergyMap, we convert each entry to the root equivalent and
    // store on the PDS

    // Set ownership of objects
    // Key - we don't own
    // Value - we own the container (but not its objects)
    calRec->getCalEventEnergyMap()->SetOwnerKeyValue(kFALSE, kTRUE);

    // We start by iterating over the keys in the CalEventEnergyMap
    for(Event::CalEventEnergyMap::iterator energyMapTdsItr  = energyMapTds->begin();
                                           energyMapTdsItr != energyMapTds->end();
                                           energyMapTdsItr++)
    {
        // Recover the key to this entry
        const Event::CalCluster* clusterTds = energyMapTdsItr->first;

        // Get a pointer to a root version
        if (m_common.m_calClusterMap.find(clusterTds) != m_common.m_calClusterMap.end()) 
        {
            TRef        clusterRootRef  = m_common.m_calClusterMap[clusterTds];
            CalCluster* clusterRoot     = (CalCluster*)clusterRootRef.GetObject();

            // Recover the collection
            const Event::CalEventEnergyVec& energyVecTds = energyMapTdsItr->second;

            // Create the TObjArray to hold the CalEventEnergy objects for this cluster
            TObjArray* energyVecRoot = new TObjArray();

            energyVecRoot->Clear();

            // Now we loop over the cluster vector
            for(Event::CalEventEnergyVec::const_iterator energyVecTdsItr  = energyVecTds.begin();
                                                          energyVecTdsItr != energyVecTds.end();
                                                          energyVecTdsItr++)
            {
                // Recover pointer to the actual Cal Cluster
                const Event::CalEventEnergy* eventEnergyTds = *energyVecTdsItr;

                CalEventEnergy * eventEnergyRoot = new CalEventEnergy;
                RootPersistence::convert(*eventEnergyTds,*eventEnergyRoot) ;

                // Keep the relationship between the TDS and root objects
                TRef ref = eventEnergyRoot;
                m_common.m_calEventEnergyMap[eventEnergyTds] = ref;

                // Store in this key's TObjArray
                energyVecRoot->Add(eventEnergyRoot);
            }

            // Signal that we are not the owner of these values (for cleanup)
            energyVecRoot->SetOwner(kFALSE);

            // Ok, now store this in our calRecon object
            calRec->getCalEventEnergyMap()->Add(clusterRoot, energyVecRoot);
        }
    }

    return;
}

void reconRootWriterAlg::fillCalCluster(CalRecon *calRec, Event::CalClusterCol* clusterColTds) {   
    // Purpose and Method:  Given the CalCluster collection from the TDS, we fill the ROOT   
    //  CalCluster collection.   
    
    unsigned int numClusters = clusterColTds->size();   
    unsigned int iCluster;   
    for (iCluster = 0; iCluster < numClusters; iCluster++)
    {
        Event::CalCluster * clusterTds = (*clusterColTds)[iCluster] ;
        CalCluster * clusterRoot = new CalCluster ;
        RootPersistence::convert(*clusterTds,*clusterRoot) ;
        calRec->addCalCluster(clusterRoot) ;   

        // Keep the relationship between the TDS and root objects
        TRef ref = clusterRoot;
        m_common.m_calClusterMap[clusterTds] = ref;
    }   
    
    return;   
}   

void reconRootWriterAlg::fillCalXtalRec(CalRecon *calRec, Event::CalXtalRecCol* xtalColTds) {   
    // Purpose and Method:  Given the CalXtalRecData collection from the TDS,   
    //   this method fills the ROOT CalXtalRecData collection.   
    
    MsgStream log(msgSvc(), name());
    Event::CalXtalRecCol::const_iterator xtalTds;   
    
    for (xtalTds = xtalColTds->begin(); xtalTds != xtalColTds->end(); xtalTds++) {
 
        
        idents::CalXtalId::CalTrigMode modeTds = (*xtalTds)->getMode() ;   
        if (modeTds!=idents::CalXtalId::BESTRANGE) {   
            int range;   
            for ( range = idents::CalXtalId::LEX8 ;
                  range <= idents::CalXtalId::HEX1 ;
                  ++range) {   
                if (!(*xtalTds)->getRangeRecData(range)) {
                    log<<MSG::DEBUG ;
                    if (log.isActive())
                      log.stream()<<"xtal for range "<<range<<" does not exist." ;
                    log<<endreq ;
                }
            }
        }

        CalXtalRecData * xtalRoot = new CalXtalRecData() ;
        RootPersistence::convert(**xtalTds,*xtalRoot) ;        
        calRec->addXtalRecData(xtalRoot) ;   

        // Keep the relationship between the TDS and root objects
        TRef ref = xtalRoot;
        const Event::CalXtalRecData* xtalRecData = *xtalTds;
        m_common.m_calXtalRecDataMap[xtalRecData] = ref;
    }   
    
    return;   
} 

//CL: 08/22/06    @@@@@@@@@@@@@@@@@@@@@


void reconRootWriterAlg::fillGcrXtal(CalRecon *calRec, Event::GcrXtalCol* gcrXtalColTds) {   
    // Purpose and Method:  Given the GcrXtal collection from the TDS,   
    //   this method fills the ROOT CalXtalRecData collection.   
    
    MsgStream log(msgSvc(), name());
    Event::GcrXtalCol::const_iterator gcrXtalColIterTds;
    
    log << MSG::DEBUG << "reconRootWriterAlg::fillGcrXtal BEGIN, gcrXtalColTds->size()=" << gcrXtalColTds->size()<< endreq;   
    
    for (gcrXtalColIterTds = gcrXtalColTds->begin(); gcrXtalColIterTds != gcrXtalColTds->end(); gcrXtalColIterTds++) {
    
    //GcrXtal* gcrXtalRoot = new GcrXtal();
    GcrXtal* gcrXtalRoot = calRec->addGcrXtal();
    RootPersistence::convert(**gcrXtalColIterTds,*gcrXtalRoot) ; 
    
       // calRec->addGcrXtal(gcrXtalRoot) ;   
    }   
    

    
    return;   
}


void reconRootWriterAlg::fillGcrTrack(CalRecon *calRec, Event::GcrTrack* gcrTrackTds) {   
    // Purpose and Method:  Given the GcrXtal collection from the TDS,   
    //   this method fills the ROOT CalXtalRecData collection.   
    
    MsgStream log(msgSvc(), name());
    
   //log << MSG::INFO << "reconRootWriterAlg::fillGcrTrack BEGIN" << endreq;   
    
    
    GcrTrack* gcrTrackRoot = new GcrTrack();
    
    RootPersistence::convert(*gcrTrackTds,*gcrTrackRoot);
    
    calRec->addGcrTrack(gcrTrackRoot);
    
    TObject* toto = calRec->getGcrTrack();
    GcrTrack* totoCaste= (GcrTrack*)toto;
    
 
   //log << MSG::INFO << "reconRootWriterAlg::fillGcrTrack END" << endreq;   
    
    return;   
}




// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 


void reconRootWriterAlg::fillCalMipTrack(CalRecon *calRec, Event::CalMipTrackCol* calMipTrackColTds) 
{   
    // Purpose and Method:  Given the CalMipTrack collection from the TDS, we fill the ROOT   
    //  CalMiptrack collection.   
    //StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << " SG : fillCalMipTrack in reconRootWriterAlg.cxx" << endreq;

    // Loop over CalMipTracks in the TDS
    for(Event::CalMipTrackColItr tdsMipIter = calMipTrackColTds->begin(); tdsMipIter != calMipTrackColTds->end(); tdsMipIter++)
    {
        CalMipTrack* calMipTrackRoot = new CalMipTrack();
        Event::CalMipTrack* calMipTrackTds = *tdsMipIter;   
        RootPersistence::convert(*calMipTrackTds,*calMipTrackRoot) ;
        calRec->addCalMipTrack(calMipTrackRoot);
    }

    log << MSG::DEBUG << " SG : fillCalMipTrack - End" << endreq;
    return;
}   


void reconRootWriterAlg::fillCalEventEnergy(CalRecon *calRec, Event::CalEventEnergyCol* calEventEnergyCol) {

//    // David C. : currently, there is only one CalEventEnergy in the TDS,
//    // yet, I prefered to consider CalEventEnergy as a usual objet on
//    // the ROOT side, that is why in the root tree there is a collection
//    // of CalEventEnergy. This collection will always have a single element,
//    // until a change is made in the TDS.
    
    unsigned int numEnergy = calEventEnergyCol->size();
    unsigned int iEnergy;
    for (iEnergy = 0; iEnergy < numEnergy; iEnergy++)
    {
        Event::CalEventEnergy *eventEnergyTds = (*calEventEnergyCol)[iEnergy] ;
        CalEventEnergy * eventEnergyRoot = new CalEventEnergy;
        RootPersistence::convert(*eventEnergyTds,*eventEnergyRoot) ;
//        calRec->addCalEventEnergy(eventEnergyRoot) ;

        // Keep the relationship between the TDS and root objects
        TRef ref = eventEnergyRoot;
        m_common.m_calEventEnergyMap[eventEnergyTds] = ref;
    }
    return ;
}

StatusCode reconRootWriterAlg::writeAcdRecon() 
{
    AcdRecon* acdRec = m_reconEvt->getAcdRecon();

    // Standard ACD Recon
    if (acdRec)
    {
        SmartDataPtr<Event::AcdRecon> acdRecTds(eventSvc(), EventModel::AcdRecon::Event);  
        if (acdRecTds) RootPersistence::convert(*acdRecTds,*acdRec);
    }

    AcdReconV2* acdRecV2 = m_reconEvt->getAcdReconV2();
    if (acdRecV2)
    {
        SmartDataPtr<Event::AcdReconV2> acdRecTdsV2(eventSvc(), EventModel::AcdReconV2::Event);
        if (acdRecTdsV2) RootPersistence::convert(*acdRecTdsV2,*acdRecV2) ;
    }

    return StatusCode::SUCCESS;
}

StatusCode reconRootWriterAlg::writeAdfRecon()
{
    SmartDataPtr<AncillaryData::Recon> adfTds(eventSvc(), "/Event/AncillaryEvent/Recon");
    if (!adfTds) return StatusCode::SUCCESS;
    m_reconEvt->initAdf(new reconRootData::AdfRecon);
    reconRootData::AdfRecon *adfRoot = m_reconEvt->getAdfRecon();
    RootPersistence::convert(*adfTds, *adfRoot);
    return StatusCode::SUCCESS;

}

void reconRootWriterAlg::writeEvent() 
{
    // Purpose and Method:  Stores the DigiEvent data for this event in the ROOT
    //    tree.
    
    m_rootIoSvc->fillTree("recon");

    return;
}

void reconRootWriterAlg::close() 
{
    // Purpose and Method:  Writes the ROOT file at the end of the run.
    //    The TObject::kWriteDelete parameter is used in the Write method
   //    Used rather than TObject::kOverwrite - supposed to be safer but slower
    //    since ROOT will periodically write to the ROOT file when the bufSize
    //    is filled.  Writing would create 2 copies of the same tree to be
    //    stored in the ROOT file, if we did not specify kOverwrite.

    m_rootIoSvc->closeFile("recon");

    return;
}

StatusCode reconRootWriterAlg::finalize()
{
    MsgStream log(msgSvc(), name());
    
    // ADDED FOR THE FILE HEADERS DEMO
    m_headersTool->writeReconHeader(m_reconTree->GetCurrentFile()) ;
    
    close();
    
    StatusCode sc = StatusCode::SUCCESS;
    //setFinalized(); No longer available in Gaudi v21r7

    log << MSG::DEBUG;
    if( log.isActive()) log.stream() << "Finalized";
    log << endreq;

    return sc;
}


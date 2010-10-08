#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"

#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"

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
// HMK-v12r7 #include "RootConvert/Recon/TkrTruncationInfoConvert.h"

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
* $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/reconRootWriterAlg.cxx,v 1.90.14.1 2010/08/31 03:08:51 heather Exp $
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
    // HMK-v12r7 void fillTruncationInfo(TkrRecon* recon, Event::TkrTruncationInfo* truncationInfoTds);
    void fillFitTracks( TkrRecon* recon, Event::TkrTrackCol* tracksTds);
    void fillVertices( TkrRecon* recon, Event::TkrVertexCol*   verticesTds, Event::TkrTrackCol* tracksTds);
    
    TkrTrackHit*   convertTkrTrackHit(const Event::TkrTrackHit* trackHitTds);
    TkrTrackParams convertTkrTrackParams(const Event::TkrTrackParams& trackParamsTds);
    
    /// Retrieves the CAL reconstruction data from the TDS and fills the CalRecon
    /// ROOT object
    StatusCode writeCalRecon();
    
    /// These are the methods specific to filling the pieces of the CalRecon stuff
    void fillCalCluster(CalRecon *calRec, Event::CalClusterCol* clusterColTds);   
    void fillCalXtalRec(CalRecon *calRec, Event::CalXtalRecCol* xtalColTds); 
    void fillCalMipTrack(CalRecon *calRec, Event::CalMipTrackCol* calMipTrackColTds); 
    void fillCalEventEnergy(CalRecon *calRec, Event::CalEventEnergyCol* calEventEnergyCol);
    
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

    // HMK-v12e7 SmartDataPtr<Event::TkrTruncationInfo> truncationInfoTds(eventSvc(), EventModel::TkrRecon::TkrTruncationInfo);
    // if (truncationInfoTds &&  truncationInfoTds->isTruncated() ) fillTruncationInfo(recon, truncationInfoTds);

    // Retrieve the information on fit tracks
    SmartDataPtr<Event::TkrTrackCol> tracksTds(eventSvc(), EventModel::TkrRecon::TkrTrackCol);
    
    // Fill the fit tracks
    if (tracksTds) fillFitTracks(recon, tracksTds);
    
    // Retrieve the information on vertices
    SmartDataPtr<Event::TkrVertexCol> verticesTds(eventSvc(), EventModel::TkrRecon::TkrVertexCol);
    
    // Fill the vertices
    if (verticesTds && tracksTds) fillVertices(recon, verticesTds, tracksTds);

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

/* HMK-v12r7
void reconRootWriterAlg::fillTruncationInfo(TkrRecon* recon, Event::TkrTruncationInfo* truncationInfoTds)
{
//Purpose and Method : fill ROOT persistent version of the truncation info with the TDS class.
  Event::TkrTruncationInfo::TkrTruncationMap* truncMap = truncationInfoTds->getTruncationMap();
  Event::TkrTruncationInfo::TkrTruncationMap::const_iterator iter = truncMap->begin();
  for(; iter!=truncMap->end(); ++iter) {
    TkrTruncationData* truncationDataRoot = new TkrTruncationData();
    RootPersistence::convert(*iter, *truncationDataRoot);
    recon->addTruncationData(truncationDataRoot);

  }
}
*/

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

        // Ok, now add the track to the list!
        recon->addTrack(trackRoot);
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
    
    // Retrieve the cal cluster collection   
    SmartDataPtr<Event::CalClusterCol> clusterColTds(eventSvc(), EventModel::CalRecon::CalClusterCol);   
    if (clusterColTds) fillCalCluster(calRec, clusterColTds);   
    
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
      calRec->addCalEventEnergy(eventEnergyRoot) ;
    }
    return ;
}

StatusCode reconRootWriterAlg::writeAcdRecon() 
{
    AcdRecon* acdRec = m_reconEvt->getAcdRecon();   
    if (!acdRec) return StatusCode::FAILURE;
    SmartDataPtr<Event::AcdRecon> acdRecTds(eventSvc(), EventModel::AcdRecon::Event);  
    if (!acdRecTds) return StatusCode::SUCCESS;
    RootPersistence::convert(*acdRecTds,*acdRec) ;

    AcdReconV2* acdRecV2 = m_reconEvt->getAcdReconV2();
    if (!acdRecV2) return StatusCode::FAILURE;
    SmartDataPtr<Event::AcdReconV2> acdRecTdsV2(eventSvc(), EventModel::AcdReconV2::Event);
    if (!acdRecTdsV2) return StatusCode::SUCCESS;
    RootPersistence::convert(*acdRecTdsV2,*acdRecV2) ;

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


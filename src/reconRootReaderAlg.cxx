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
#include "Event/Recon/TkrRecon/TkrPatCand.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/CalRecon/CalCluster.h"   
#include "Event/Recon/CalRecon/CalXtalRecData.h"   

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TObjArray.h"
#include "TCollection.h"  // Declares TIter

#include "reconRootData/ReconEvent.h"

#include "facilities/Util.h"

#include <vector>
#include <map>

/** @class reconRootReaderAlg
 * @brief Reads Reconstruction data from a persistent ROOT file and stores the
 * the data in the TDS.
 *
 * @author Heather Kelly
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/reconRootReaderAlg.cxx,v 1.10 2002/10/16 02:02:06 heather Exp $
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

    StatusCode storeTkrCandidateTrackCol(TkrRecon *tkrRecRoot);

    StatusCode storeTrackAndVertexCol(TkrRecon *tkrRecRoot, bool vertexOnTdsFlag);

    /// convert a ROOT TkrFitHit to a TDS Event::TkrFitHit
    Event::TkrFitHit convertTkrFitHit(const TkrFitHit& hitRoot);
    /// convert a ROOT TkrCovMat to a CLHEP HepMatrix
    void convertMatrix(const TkrCovMat& matRoot, HepMatrix &matTds);
    /// convert a ROOT TkrCandHit::AXIS to a TDS Event::TkrCluster::view
    void convertCandHitView(TkrCandHit::AXIS viewRoot, Event::TkrCluster::view &viewTds);
    /// convert a ROOT TkrHitPlane::AXIS to a TDS Event::TkrCluster::view
    void convertHitPlaneView(TkrHitPlane::AXIS viewRoot, Event::TkrCluster::view &viewTds);

    /// Reads CAL recon data from ROOT and puts data on TDS
    StatusCode readCalRecon();

    /// read in CAL xtal recon data from ROOT and store on the TDS
    StatusCode storeCalXtalRecDataCol(CalRecon *calRecRoot);

    /// read CAL cluster data from ROOT and store on TDS
    StatusCode storeCalClusterCol(CalRecon *calRecRoot);

    /// Reads ACD recon data from ROOT and puts data on the TDS
    StatusCode readAcdRecon();

    /// Closes the ROOT file
    void close();
   
    /// ROOT file pointer
    TFile *m_reconFile;
    /// ROOT tree pointer
    TTree *m_reconTree;
    /// Top-level Monte Carlo ROOT object
    ReconEvent *m_reconEvt;
    /// name of the output ROOT file
    std::string m_fileName;
    /// name of the Recon TTree stored in the ROOT file
    std::string m_treeName;
    /// Number of events in the input ROOT TTree
    int m_numEvents;

};

static const AlgFactory<reconRootReaderAlg>  Factory;
const IAlgFactory& reconRootReaderAlgFactory = Factory;


reconRootReaderAlg::reconRootReaderAlg(const std::string& name, ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator)
{
    // Input pararmeters that may be set via the jobOptions file
    // Input ROOT file name
    declareProperty("reconRootFile",m_fileName="recon.root");
    // Input TTree name
    declareProperty("reconTreeName", m_treeName="Recon");

}

StatusCode reconRootReaderAlg::initialize()
{
    // Purpose and Method:  Called once before the run begins.  This method
    //    opens a new ROOT file and prepares for reading.

    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    // Use the Job options service to set the Algorithm's parameters
    // This will retrieve parameters set in the job options file
    setProperties();

    facilities::Util::expandEnvVar(&m_fileName);
    
    // Save the current directory for the ntuple writer service
    TDirectory *saveDir = gDirectory;   
    m_reconFile = new TFile(m_fileName.c_str(), "READ");
    if (!m_reconFile->IsOpen()) {
        log << MSG::ERROR << "ROOT file " << m_fileName 
            << " could not be opened for reading." << endreq;
        return StatusCode::FAILURE;
    }
    m_reconFile->cd();
    m_reconTree = (TTree*)m_reconFile->Get(m_treeName.c_str());
    if (!m_reconTree) {
        log << MSG::ERROR << "Could not load Tree " << m_treeName << 
            " from file " << m_fileName << endreq;
        return StatusCode::FAILURE;
    }

    m_reconEvt = 0;
    m_reconTree->SetBranchAddress("ReconEvent", &m_reconEvt);
    
    m_numEvents = m_reconTree->GetEntries();

    saveDir->cd();
    return sc;
    
}

StatusCode reconRootReaderAlg::execute()
{
    // Purpose and Method:  Called once per event.  This method calls
    //   the appropriate methods to read data from the ROOT file and store
    //   data on the TDS.

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    if (!m_reconFile->IsOpen()) {
        log << MSG::ERROR << "ROOT file " << m_fileName 
            << " could not be opened for reading." << endreq;
        return StatusCode::FAILURE;
    }

    static UInt_t evtId = 0;
    if (evtId >= m_numEvents) {
        log << MSG::ERROR << "ROOT file contains no more events" << endreq;
        return StatusCode::FAILURE;
    }

    m_reconTree->GetEvent(evtId);

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

    m_reconEvt->Clear();
    evtId++;
    
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
		log << MSG::INFO << "Tkr Cluster Collection is already on the TDS" << endreq;
	} else {
		sc = storeTkrClusterCol(tkrRecRoot);
		if (sc.isFailure()) {
			log << MSG::ERROR << "failed to store TKR cluster collection on TDS" << endreq;
			return sc;
		}
	}

    // check to see if track candidates exist on the TDS already
    // Store candidate track collection, only if it does not already exist on TDS
    SmartDataPtr<Event::TkrPatCandCol> candidatesColTds(eventSvc(), EventModel::TkrRecon::TkrPatCandCol);
    if (candidatesColTds) {
        log << MSG::INFO << "Tkr CandidateTracks Col is already on TDS" << endreq;
    } else {
        sc = storeTkrCandidateTrackCol(tkrRecRoot);
        if (sc.isFailure()) {
            log << MSG::ERROR << "failed to store candidate tracks on TDS" << endreq;
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
    SmartDataPtr<Event::TkrFitTrackCol> trackColTds(eventSvc(), EventModel::TkrRecon::TkrFitTrackCol);
    if (trackColTds) {
        log << MSG::INFO << "Tkr FitTrackCol is already on TDS" << endreq;
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

	while (clusterRoot = (TkrCluster*)clusterIter.Next()) {
		TkrCluster::view viewRoot = clusterRoot->getView();
                Event::TkrCluster::view viewTds;
		
		if (viewRoot == TkrCluster::X) viewTds = Event::TkrCluster::X;
		else viewTds = Event::TkrCluster::Y;

		TVector3 posRoot = clusterRoot->getPosition();
		Point posTds(posRoot.X(), posRoot.Y(), posRoot.Z());

		Event::TkrCluster *clusterTds = new Event::TkrCluster(clusterRoot->getId(),
			clusterRoot->getPlane(), viewTds, clusterRoot->getFirstStrip(),
			clusterRoot->getLastStrip(), posTds, 0.0, clusterRoot->getTower());

		clusterTdsCol->addCluster(clusterTds);
	}

	sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrClusterCol, clusterTdsCol);
    if (sc.isFailure()) {
        log << MSG::DEBUG << "Failed to register TkrClusterCol" << endreq;
        return StatusCode::FAILURE;
    }


	return sc;
}

StatusCode reconRootReaderAlg::storeTkrCandidateTrackCol(TkrRecon *tkrRecRoot) {
    // Purpose and Method: Read in candidate tracks from ROOT and store them on the
    //   TDS.

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Create TDS Candidate track collection
    Event::TkrPatCandCol *candTdsCol = new Event::TkrPatCandCol;

    // Retrieve ROOT candidate track collection
    const TObjArray *candTrackColRoot = tkrRecRoot->getTrackCandCol();
    TIter candTrackIter(candTrackColRoot);
    TkrCandTrack *candTrackRoot = 0;

    while(candTrackRoot = (TkrCandTrack*)candTrackIter.Next()) {
        TVector3 posRoot = candTrackRoot->getPosition();
        TVector3 dirRoot = candTrackRoot->getDirection();
        Ray rayTds( Point(posRoot.X(), posRoot.Y(), posRoot.Z()), 
            Vector(dirRoot.X(), dirRoot.Y(), dirRoot.Z()) );
       
        Event::TkrPatCand *track = new Event::TkrPatCand(
            candTrackRoot->getLayer(), candTrackRoot->getTower(),
            candTrackRoot->getEnergy(), candTrackRoot->getQuality(),
            rayTds); 

        std::vector<TkrCandHit>::const_iterator hitIt;
        for (hitIt = candTrackRoot->getHitIterBegin(); hitIt != candTrackRoot->getHitIterEnd(); hitIt++) {
            TVector3 posRoot = hitIt->getPosition();
            Point posTds(posRoot.X(), posRoot.Y(), posRoot.Z());
            Event::TkrCluster::view viewTds;
            convertCandHitView(hitIt->getView(), viewTds);
           
            Event::TkrPatCandHit hitTds(
                hitIt->getHitIndex(),
                posTds,
                hitIt->getTower(),
                hitIt->getPlane(),
                viewTds);

            track->addCandHit(hitTds);
        }
        candTdsCol->push_back(track);
    }

    sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrPatCandCol, candTdsCol);
    if (sc.isFailure()) {
        log << MSG::DEBUG << "Failed to register Tkr CandTrackCol" << endreq;
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
    std::map<int, Event::TkrFitTrack*> trackMap;
    trackMap.clear();

    // Create TDS version of track collection
    Event::TkrFitTrackCol *trackTdsCol = new Event::TkrFitTrackCol;

    // Retrieve ROOT version of track collection
    const TObjArray *trackRootCol = tkrRecRoot->getTrackCol();
    TIter trackIter(trackRootCol);
    TkrTrack *trackRoot = 0;

    while (trackRoot = (TkrTrack*)trackIter.Next()) {
        Event::TkrFitTrack *trackTds = new Event::TkrFitTrack();

        trackTds->initializeInfo(trackRoot->getXgaps(), trackRoot->getYgaps(),
            trackRoot->getXistGaps(), trackRoot->getYistGaps());
        trackTds->initializeQual(trackRoot->getChiSq(), trackRoot->getChiSqSmooth(),
            trackRoot->getRmsResid(), trackRoot->getQuality(), 
            trackRoot->getKalEnergy(), trackRoot->getKalThetaMS());

        // Store track id to properly associate tracks and vertices later
        trackMap[trackRoot->getId()] = trackTds;

        TkrHitPlaneIter hitItRoot;
        for (hitItRoot = trackRoot->getHitIterBegin(); hitItRoot != trackRoot->getHitIterEnd(); hitItRoot++) {
            Event::TkrFitPlane planeTds;
            Event::TkrCluster::view projTds, projPlusTds;
            convertHitPlaneView(hitItRoot->getProjection(), projTds);
            convertHitPlaneView(hitItRoot->getProjPlus(), projPlusTds);

            planeTds.initializeInfo(hitItRoot->getIdHit(), hitItRoot->getIdTower(), hitItRoot->getIdPlane(),
                projTds, projPlusTds, 
                hitItRoot->getZplane(), hitItRoot->getEnePlane(), hitItRoot->getRadLen(),
                hitItRoot->getActiveDist());

            Event::TkrFitHit measTds = convertTkrFitHit(hitItRoot->getHitMeas());
            Event::TkrFitHit predTds = convertTkrFitHit(hitItRoot->getHitPred());
            Event::TkrFitHit fitTds = convertTkrFitHit(hitItRoot->getHitFit());
            Event::TkrFitHit smoothTds = convertTkrFitHit(hitItRoot->getHitSmooth());
            TkrCovMat matRoot = hitItRoot->getQmaterial();
            HepMatrix hepMat(4, 4, 0);
            convertMatrix(matRoot, hepMat);
            Event::TkrFitMatrix matTds(hepMat);
            planeTds.initializeHits(measTds, predTds, fitTds, smoothTds, matTds);
        }

        trackTdsCol->push_back(trackTds);
    }

    sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrFitTrackCol, trackTdsCol);
    if (sc.isFailure()) {
        log << MSG::DEBUG << "Failed to register Tkr FitTrackCol" << endreq;
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

    while (vertexRoot = (TkrVertex*)vertexIter.Next()) {
        TVector3 posRoot = vertexRoot->getPosition();
        TVector3 dirRoot = vertexRoot->getDirection();
        Ray rayTds(Point(posRoot.X(), posRoot.Y(), posRoot.Z()), 
            Vector(dirRoot.X(), dirRoot.Y(), dirRoot.Z()));

        Event::TkrVertex *vertexTds = new Event::TkrVertex(
            vertexRoot->getLayer(), vertexRoot->getTower(),
            vertexRoot->getEnergy(), vertexRoot->getQuality(), rayTds);

        unsigned int numTracks = vertexRoot->getNumTracks();
        unsigned int iTrack;
        for (iTrack = 0; iTrack < numTracks; iTrack++) {
            unsigned int id = vertexRoot->getTrackId(iTrack);
            vertexTds->addTrack(trackMap[id]);
        }

        vertexColTds->push_back(vertexTds);
    }

    sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrVertexCol, vertexColTds);
    if (sc.isFailure()) {
        log << MSG::DEBUG << "Failed to register Tkr VertexCol" << endreq;
        return StatusCode::FAILURE;
    }

    return sc;
}


void reconRootReaderAlg::convertCandHitView(TkrCandHit::AXIS viewRoot, Event::TkrCluster::view &viewTds) {
    // Purpose and Method:  Converts from ROOT TkrCandHit::AXIS to a Event::TkrCluster::view

    viewTds = (viewRoot == TkrCandHit::X) ? Event::TkrCluster::X : Event::TkrCluster::Y;
    if (viewRoot == TkrCandHit::NONE) viewTds = Event::TkrCluster::XY;
    return;
}

void reconRootReaderAlg::convertHitPlaneView(TkrHitPlane::AXIS viewRoot, Event::TkrCluster::view &viewTds) {
    // Purpose and Method:  Convertes from ROOT TkrHitPlane::AXIS to an Event::TkrCluster::view

    viewTds = (viewRoot == TkrHitPlane::X) ? Event::TkrCluster::X : Event::TkrCluster::Y;
    if (viewRoot == TkrHitPlane::NONE) viewTds = Event::TkrCluster::XY;
    return;
}

Event::TkrFitHit reconRootReaderAlg::convertTkrFitHit(const TkrFitHit& hitRoot) {
    // Purpose and Method:  Converts from ROOT TkrFitHit to a TDS Event::TkrFitHit
    TkrParams paramRoot = hitRoot.getTkrParams();
    Event::TkrFitPar paramTds(paramRoot.getXPos(), paramRoot.getXSlope(), 
        paramRoot.getYPos(), paramRoot.getYSlope());
    
    TkrCovMat matRoot = hitRoot.getTkrCovMat();
    HepMatrix hepMat(4, 4, 0);
    convertMatrix(matRoot, hepMat);
    Event::TkrFitMatrix matTds(hepMat);
    TkrFitHit::TYPE typeRoot = hitRoot.getHitType();
    return Event::TkrFitHit(Event::TkrFitHit::TYPE(hitRoot.getHitType()), paramTds, matTds);
}


void reconRootReaderAlg::convertMatrix(const TkrCovMat& matRoot, HepMatrix &hepMat) {
    // Purpose and Method:  Converts from a ROOT TkrCovMat to CLHEP HepMatrix
    hepMat[0][0] = matRoot.getCovX0X0();
    hepMat[1][1] = matRoot.getCovSxSx();
    hepMat[0][1] = matRoot.getCovX0Sx();
    hepMat[2][2] = matRoot.getCovY0Y0();
    hepMat[3][3] = matRoot.getCovSySy();
    hepMat[2][2] = matRoot.getCovY0Sy();
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

    SmartDataPtr<Event::CalXtalRecCol> xtalRecColTds(eventSvc(),EventModel::CalRecon::CalClusterCol);
    if (xtalRecColTds){
        log << MSG::INFO << "XtalRecCol data is already on the TDS" << endreq;
    } else {
        sc = storeCalXtalRecDataCol(calRecRoot);
    }

    if (sc.isFailure()) {
        log << MSG::INFO << "Failed to store CalXtalRecCol on the TDS" << endreq;
        return sc;
    }


    SmartDataPtr<Event::CalClusterCol> checkCalClusterColTds(eventSvc(),EventModel::CalRecon::CalClusterCol);
    if (checkCalClusterColTds){
        log << MSG::INFO << "CalClusterCol data is already on the TDS" << endreq;
        return sc;
    } else {
        sc = storeCalClusterCol(calRecRoot);
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
    while (calXtalRecRoot = (CalXtalRecData*)calXtalIter.Next()) {
        CalXtalId idRoot = calXtalRecRoot->getPackedId();
        idents::CalXtalId idTds(idRoot.getTower(), idRoot.getLayer(), idRoot.getColumn());
        // create new object to store crystal reconstructed data 
        Event::CalXtalRecData* calXtalRecDataTds = 0;
        if (calXtalRecRoot->getMode() == CalXtalId::ALLRANGE) {
            calXtalRecDataTds = new Event::CalXtalRecData(idents::CalXtalId::ALLRANGE, idTds);
            unsigned int range;
            for (range = idents::CalXtalId::LEX8; range < idents::CalXtalId::HEX1; range++) {    
                const CalRangeRecData *xtalRangeRoot = calXtalRecRoot->getRangeRecData(range);
                Event::CalXtalRecData::CalRangeRecData *xtalRangeTds = 
                new Event::CalXtalRecData::CalRangeRecData(
                xtalRangeRoot->getRange(CalXtalId::POS), xtalRangeRoot->getEnergy(CalXtalId::POS),
                xtalRangeRoot->getRange(CalXtalId::NEG), xtalRangeRoot->getEnergy(CalXtalId::NEG));
                calXtalRecDataTds->addRangeRecData(*xtalRangeTds);
            }
        } else if (calXtalRecRoot->getMode() == CalXtalId::BESTRANGE) {

            calXtalRecDataTds = new Event::CalXtalRecData(idents::CalXtalId::BESTRANGE, idTds);
            const CalRangeRecData *xtalRangeRoot = calXtalRecRoot->getRangeRecData(0);   

            Event::CalXtalRecData::CalRangeRecData *xtalRangeTds = 
                new Event::CalXtalRecData::CalRangeRecData(
                xtalRangeRoot->getRange(CalXtalId::POS), xtalRangeRoot->getEnergy(CalXtalId::POS),
                xtalRangeRoot->getRange(CalXtalId::NEG), xtalRangeRoot->getEnergy(CalXtalId::NEG));
            calXtalRecDataTds->addRangeRecData(*xtalRangeTds);
        }

        calXtalRecColTds->push_back(calXtalRecDataTds);
    }

    //register output data collection as a TDS object
    sc = eventSvc()->registerObject(EventModel::CalRecon::CalXtalRecCol, calXtalRecColTds);
    if (sc.isFailure()) {
        log << MSG::DEBUG << "Failed to register CalXtalRecCol" << endreq;
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
    
    while (calClusterRoot = (CalCluster*)calClusterIter.Next()) {
        TVector3 calPosRoot = calClusterRoot->getPosition();
        Point pos(calPosRoot.X(), calPosRoot.Y(), calPosRoot.Z());
        Event::CalCluster *calClusterTds = new Event::CalCluster(
            calClusterRoot->getEnergySum(),
            pos);
        TVector3 dirRoot = calClusterRoot->getDirection();
        Vector dirTds(dirRoot.X(), dirRoot.Y(), dirRoot.Z());
        std::vector<Vector> posLayerTds;
        std::vector<TVector3> posLayerRoot = calClusterRoot->getPosLayer();
        std::vector<TVector3>::const_iterator posLayerIt;
        for (posLayerIt = posLayerRoot.begin(); posLayerIt != posLayerRoot.end(); posLayerIt++) {
            Vector vecTds(posLayerIt->X(), posLayerIt->Y(), posLayerIt->Z());
            posLayerTds.push_back(vecTds);
        }
        
        std::vector<Vector> rmsLayerTds;
        std::vector<TVector3> rmsLayerRoot = calClusterRoot->getRmsLayer();
        std::vector<TVector3>::const_iterator rmsLayerIt;
        for (rmsLayerIt = rmsLayerRoot.begin(); rmsLayerIt != rmsLayerRoot.end(); rmsLayerIt++) {
            Vector vecTds(rmsLayerIt->X(), rmsLayerIt->Y(), rmsLayerIt->Z());
            rmsLayerTds.push_back(vecTds);
        }

        calClusterTds->initialize(
            calClusterRoot->getEnergyLeak(),
            calClusterRoot->getEneLayer(),
            posLayerTds,
            rmsLayerTds,
            calClusterRoot->getRmsLong(),
            calClusterRoot->getRmsTrans(),
            dirTds,
            calClusterRoot->getTransvOffset());
        calClusterTds->initProfile(
            calClusterRoot->getFitEnergy(),
            calClusterRoot->getProfChisq(),
            calClusterRoot->getCsiStart(),
            calClusterRoot->getCsiAlpha(),
            calClusterRoot->getCsiLambda());

        calClusterColTds->push_back(calClusterTds);
    }
    
    sc = eventSvc()->registerObject(EventModel::CalRecon::CalClusterCol, calClusterColTds);

    return sc;
}



StatusCode reconRootReaderAlg::readAcdRecon() {
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    const AcdRecon *acdRecRoot = m_reconEvt->getAcdRecon();
    if (!acdRecRoot) {
        log << MSG::DEBUG << "No AcdRecon found in ROOT file" << endreq;
        return StatusCode::SUCCESS;
    }

    SmartDataPtr<Event::AcdRecon> checkAcdRecTds(eventSvc(), EventModel::AcdRecon::Event);  
    if (checkAcdRecTds) {
        log << MSG::INFO << "AcdRecon data already on TDS!" << endreq;
        return StatusCode::SUCCESS;
    }

    // create the TDS location for the AcdRecon
    const AcdId acdIdRoot = acdRecRoot->getMinDocaId();
    const idents::AcdId acdIdTds(acdIdRoot.getLayer(), acdIdRoot.getFace(), 
        acdIdRoot.getRow(), acdIdRoot.getColumn());
	std::vector<idents::AcdId> idColTds;
	std::vector<AcdId>::const_iterator idRootIt;
	for (idRootIt = acdRecRoot->getIdCol().begin(); idRootIt != acdRecRoot->getIdCol().end(); idRootIt++) {
		idColTds.push_back(idents::AcdId(idRootIt->getLayer(), idRootIt->getFace(),
			idRootIt->getRow(), idRootIt->getColumn()));
	}
    std::vector<double> energyColTds = acdRecRoot->getEnergyCol();
	Event::AcdRecon *acdRecTds = new Event::AcdRecon(acdRecRoot->getEnergy(), acdRecRoot->getTileCount(),
        acdRecRoot->getGammaDoca(), acdRecRoot->getDoca(), 
        acdRecRoot->getActiveDist(), acdIdTds, 
        acdRecRoot->getRowDocaCol(), acdRecRoot->getRowActDistCol(), idColTds, energyColTds);
    
    sc = eventSvc()->registerObject(EventModel::AcdRecon::Event, acdRecTds);
    if (sc.isFailure()) {
        log << MSG::DEBUG << "Failed to register AcdRecon" << endreq;
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

    TDirectory *saveDir = gDirectory;
    m_reconFile->cd();
    m_reconFile->Close();
    saveDir->cd();
}

StatusCode reconRootReaderAlg::finalize()
{
    close();
    
    StatusCode sc = StatusCode::SUCCESS;
    return sc;
}
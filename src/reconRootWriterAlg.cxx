#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"

#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "Event/Recon/TkrRecon/TkrPatCandCol.h"
#include "Event/Recon/TkrRecon/TkrFitTrackCol.h"
#include "Event/Recon/TkrRecon/TkrVertexCol.h"

#include "idents/CalXtalId.h"

#include "facilities/Util.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"

#include "reconRootData/ReconEvent.h"

/** @class reconRootWriterAlg
 * @brief Writes Recon TDS data to a persistent ROOT file.
 *
 * @author Heather Kelly
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/reconRootWriterAlg.cxx,v 1.1 2002/05/15 22:30:54 heather Exp $
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
    void       fillFitTracks(TkrRecon* recon, Event::TkrFitTrackCol* tracksTds);

    /// Retrieves the CAL reconstruction data from the TDS and fills the CalRecon
    /// ROOT object
    StatusCode writeCalRecon();

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
    
};


static const AlgFactory<reconRootWriterAlg>  Factory;
const IAlgFactory& reconRootWriterAlgFactory = Factory;

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

}

StatusCode reconRootWriterAlg::initialize()
{
    // Purpose and Method:  Called once before the run begins.  This method
    //    opens a new ROOT file and prepares for writing.

    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    // Use the Job options service to set the Algorithm's parameters
    // This will retrieve parameters set in the job options file
    setProperties();
    
    facilities::Util::expandEnvVar(&m_fileName);

    // Save the current directory for the ntuple writer service
    TDirectory *saveDir = gDirectory;   
    // Create the new ROOT file
    m_reconFile = new TFile(m_fileName.c_str(), "RECREATE");
    if (!m_reconFile->IsOpen()) sc = StatusCode::FAILURE;
    m_reconFile->cd();
    m_reconFile->SetCompressionLevel(m_compressionLevel);
    m_reconTree = new TTree(m_treeName.c_str(), "GLAST Reconstruction Data");
    m_reconEvt = new ReconEvent();
    m_reconTree->Branch("Recon","Recon", &m_reconEvt, m_bufSize, m_splitMode);
    
    saveDir->cd();
    return sc;
    
}

StatusCode reconRootWriterAlg::execute()
{
    // Purpose and Method:  Called once per event.  This method calls
    //   the appropriate methods to read data from the TDS and write data
    //   to the ROOT file.

    MsgStream log(msgSvc(), name());

    StatusCode sc = StatusCode::SUCCESS;

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
    evtTds->fillStream(log.stream());
    log << endreq;

    m_reconEvt->initialize(evtId, runId);

    return sc;
}

StatusCode reconRootWriterAlg::writeTkrRecon() {
    // Purpose and Method:  Retrieve the Tkr Recon data from the TDS

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Get pointer to TkrRecon part of ReconEvent
    TkrRecon* recon = m_reconEvt->getTkrRecon();

    // Retrieve the information on fit tracks
    SmartDataPtr<Event::TkrFitTrackCol> tracksTds(eventSvc(), EventModel::TkrRecon::TkrFitTrackCol);

    // Fill the fit tracks
    fillFitTracks(recon, tracksTds);

    return sc;
}

void reconRootWriterAlg::fillFitTracks(TkrRecon* recon, Event::TkrFitTrackCol* tracksTds)
{
    // Purpose and Method:  This creates root tracks from tds tracks 
    //                      and adds them to the list kept in TkrRecon

    // Iterate over the tracks in the TDS
    int                 trkId  = 0;
    Event::TkrFitColPtr trkPtr = tracksTds->getTrackIterBegin();
    while(trkPtr != tracksTds->getTrackIterEnd())
    {
        Event::TkrFitTrack* trackTds  = *trkPtr++;       // The TDS track
        TkrTrack*           track     = new TkrTrack();  // Create a new Root Track

        // Initialize the track
        track->initializeInfo(trkId++,
                              trackTds->getNumXGaps(),
                              trackTds->getNumYGaps(),
                              trackTds->getNumXFirstGaps(),
                              trackTds->getNumYFirstGaps() );

        track->initializeQaul(trackTds->getChiSquare(),
                              trackTds->getChiSquareSmooth(),
                              trackTds->getScatter(),
                              trackTds->getQuality(),
                              trackTds->getKalEnergy(),
                              trackTds->getKalThetaMS() );

        // Now loop over the hit planes and fill that information
        Event::TkrFitPlaneConPtr planePtr = trackTds->getHitIterBegin();
        while(planePtr != trackTds->getHitIterEnd())
        {
            Event::TkrFitPlane  planeTds = *planePtr++;
            TkrHitPlane         plane;

            // Wouldn't it be great if there was a STANDARD definition for these!
            TkrHitPlane::AXIS   proj     = planeTds.getProjection() == Event::TkrCluster::view::X ? TkrHitPlane::AXIS::X : TkrHitPlane::AXIS::Y;
            TkrHitPlane::AXIS   projPlus = planeTds.getNextProj()   == Event::TkrCluster::view::X ? TkrHitPlane::AXIS::X : TkrHitPlane::AXIS::Y;

            plane.initializeInfo(planeTds.getIDHit(),
                                 planeTds.getIDTower(),
                                 planeTds.getIDPlane(),
                                 proj,
                                 projPlus,
                                 planeTds.getZPlane(),
                                 planeTds.getEnergy(),
                                 planeTds.getRadLen() );

            // Here we build the hit info (one at a time) starting with measured
            Event::TkrFitHit    hitTds = planeTds.getHit(Event::TkrFitHit::TYPE::MEAS);
            Event::TkrFitPar    parTds = hitTds.getPar();
            Event::TkrFitMatrix covTds = hitTds.getCov();

            TkrFitHit           measHit(TkrFitHit::TYPE::MEAS,
                                        TkrParams(parTds.getXPosition(),
                                                  parTds.getXSlope(),
                                                  parTds.getYPosition(),
                                                  parTds.getYSlope() ),
                                        TkrCovMat(covTds.getcovX0X0(),
                                                  covTds.getcovSxSx(),
                                                  covTds.getcovX0Sx(),
                                                  covTds.getcovY0Y0(),
                                                  covTds.getcovSySy(),
                                                  covTds.getcovY0Sy() ));

            // Now the predicted hit
            hitTds = planeTds.getHit(Event::TkrFitHit::TYPE::PRED);
            parTds = hitTds.getPar();
            covTds = hitTds.getCov();

            TkrFitHit           predHit(TkrFitHit::TYPE::PRED,
                                        TkrParams(parTds.getXPosition(),
                                                  parTds.getXSlope(),
                                                  parTds.getYPosition(),
                                                  parTds.getYSlope() ),
                                        TkrCovMat(covTds.getcovX0X0(),
                                                  covTds.getcovSxSx(),
                                                  covTds.getcovX0Sx(),
                                                  covTds.getcovY0Y0(),
                                                  covTds.getcovSySy(),
                                                  covTds.getcovY0Sy() ));

            // Now the fit (filtered) hit
            hitTds = planeTds.getHit(Event::TkrFitHit::TYPE::FIT);
            parTds = hitTds.getPar();
            covTds = hitTds.getCov();

            TkrFitHit           filtHit(TkrFitHit::TYPE::FIT,
                                        TkrParams(parTds.getXPosition(),
                                                  parTds.getXSlope(),
                                                  parTds.getYPosition(),
                                                  parTds.getYSlope() ),
                                        TkrCovMat(covTds.getcovX0X0(),
                                                  covTds.getcovSxSx(),
                                                  covTds.getcovX0Sx(),
                                                  covTds.getcovY0Y0(),
                                                  covTds.getcovSySy(),
                                                  covTds.getcovY0Sy() ));

            // Now the smoothed hit
            hitTds = planeTds.getHit(Event::TkrFitHit::TYPE::SMOOTH);
            parTds = hitTds.getPar();
            covTds = hitTds.getCov();

            TkrFitHit           smooHit(TkrFitHit::TYPE::SMOOTH,
                                        TkrParams(parTds.getXPosition(),
                                                  parTds.getXSlope(),
                                                  parTds.getYPosition(),
                                                  parTds.getYSlope() ),
                                        TkrCovMat(covTds.getcovX0X0(),
                                                  covTds.getcovSxSx(),
                                                  covTds.getcovX0Sx(),
                                                  covTds.getcovY0Y0(),
                                                  covTds.getcovSySy(),
                                                  covTds.getcovY0Sy() ));

            // Here retrieve the scattering matrix
            Event::TkrFitMatrix scatTds = planeTds.getQmaterial();
            TkrCovMat           scatCov(scatTds.getcovX0X0(),
                                        scatTds.getcovSxSx(),
                                        scatTds.getcovX0Sx(),
                                        scatTds.getcovY0Y0(),
                                        scatTds.getcovSySy(),
                                        scatTds.getcovY0Sy() );

            plane.initializeHits(measHit, predHit, filtHit, smooHit, scatCov);

            // Finally! Add this root track to the list
            track->addHit(plane);
        }

        // Ok, now add the track to the list!
        recon->addTrack(track);
    }

    return;
}


StatusCode reconRootWriterAlg::writeCalRecon() {
    // Purpose and Method:  Retrieve the Cal Recon data from the TDS

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;


    return sc;
}

void reconRootWriterAlg::writeEvent() 
{
    // Purpose and Method:  Stores the DigiEvent data for this event in the ROOT
    //    tree.  The m_digiEvt object is cleared for the next event.

    TDirectory *saveDir = gDirectory;
    m_reconFile->cd();
    m_reconTree->Fill();
    m_reconEvt->Clear();
    saveDir->cd();
    return;
}

void reconRootWriterAlg::close() 
{
    // Purpose and Method:  Writes the ROOT file at the end of the run.
    //    The TObject::kOverWrite parameter is used in the Write method
    //    since ROOT will periodically write to the ROOT file when the bufSize
    //    is filled.  Writing would create 2 copies of the same tree to be
    //    stored in the ROOT file, if we did not specify kOverwrite.

    TDirectory *saveDir = gDirectory;
    m_reconFile->cd();
    m_reconFile->Write(0, TObject::kOverwrite);
    m_reconFile->Close();
    saveDir->cd();
    return;
}

StatusCode reconRootWriterAlg::finalize()
{
    close();
    
    StatusCode sc = StatusCode::SUCCESS;
    return sc;
}


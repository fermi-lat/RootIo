#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"

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
 * $Header$
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


    return sc;
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


#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/DigiEvent.h"
#include "Event/Digi/CalDigi.h"

#include "idents/CalXtalId.h"

#include "facilities/Util.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"

#include "digiRootData/DigiEvent.h"

/** @class digiRootWriterAlg
 * @brief Writes Digi TDS data to a persistent ROOT file.
 *
 * @author Heather Kelly
 * $Header$
 */

class digiRootWriterAlg : public Algorithm
{	
public:
    
    digiRootWriterAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    /// Handles setup by opening ROOT file in write mode and creating a new TTree
    StatusCode initialize();
   
    /// Orchastrates reading from TDS and writing to ROOT for each event
    StatusCode execute();
    
    /// Closes the ROOT file and cleans up
    StatusCode finalize();

private:

    /// Retrieves event Id and run Id from TDS and fills the McEvent ROOT object
    StatusCode writeDigiEvent();

    /// Retrieves the CAL digitization data from the TDS and fills the CalDigi
    /// ROOT object
    StatusCode writeCalDigi();

    /// Calls TTree::Fill for each event and clears m_mcEvt
    void writeEvent();

    /// Performs the final write to the ROOT file and closes
    void close();
   
    /// ROOT file pointer
    TFile *m_digiFile;
    /// ROOT tree pointer
    TTree *m_digiTree;
    /// Top-level Monte Carlo ROOT object
    DigiEvent *m_digiEvt;
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


static const AlgFactory<digiRootWriterAlg>  Factory;
const IAlgFactory& digiRootWriterAlgFactory = Factory;

digiRootWriterAlg::digiRootWriterAlg(const std::string& name, 
                                 ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator)
{
    // Input parameters available to be set via the jobOptions file
    declareProperty("digiRootFile",m_fileName="digi.root");
    declareProperty("splitMode", m_splitMode=1);
    declareProperty("bufferSize", m_bufSize=64000);
    // ROOT default compression
    declareProperty("compressionLevel", m_compressionLevel=1);
    declareProperty("treeName", m_treeName="Digi");

}

StatusCode digiRootWriterAlg::initialize()
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
    m_digiFile = new TFile(m_fileName.c_str(), "RECREATE");
    if (!m_digiFile->IsOpen()) sc = StatusCode::FAILURE;
    m_digiFile->cd();
    m_digiFile->SetCompressionLevel(m_compressionLevel);
    m_digiTree = new TTree(m_treeName.c_str(), "GLAST Digitization Data");
    m_digiEvt = new DigiEvent();
    m_digiTree->Branch("DigiEvent","DigiEvent", &m_digiEvt, m_bufSize, m_splitMode);
    
    saveDir->cd();
    return sc;
    
}

StatusCode digiRootWriterAlg::execute()
{
    // Purpose and Method:  Called once per event.  This method calls
    //   the appropriate methods to read data from the TDS and write data
    //   to the ROOT file.

    MsgStream log(msgSvc(), name());

    StatusCode sc = StatusCode::SUCCESS;

    sc = writeDigiEvent();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to write DigiEvent" << endreq;
        return sc;
    }

    sc = writeCalDigi();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to write Cal Digi Collection" << endreq;
        return sc;
    }
   
    writeEvent();
    return sc;
}


StatusCode digiRootWriterAlg::writeDigiEvent() {
    // Purpose and Method:  Retrieve the Event object from the TDS and set the
    //    event and run numbers in the DigiEvent ROOT object

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Retrieve the Event data for this event
    SmartDataPtr<Event::EventHeader> evtTds(eventSvc(), EventModel::EventHeader);

    if (!evtTds) return sc;

    UInt_t evtId = evtTds->event();
    UInt_t runId = evtTds->run();

    Bool_t fromMc = true;

    log << MSG::DEBUG;
    evtTds->fillStream(log.stream());
    log << endreq;

    m_digiEvt->initialize(evtId, runId, fromMc);

    return sc;
}


StatusCode digiRootWriterAlg::writeCalDigi() {
    // Purpose and Method:  Retrieve the Event object from the TDS and set the
    //    event and run numbers in the DigiEvent ROOT object

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    SmartDataPtr<Event::CalDigiCol> calDigiColTds(eventSvc(), EventModel::Digi::CalDigiCol);
    if (!calDigiColTds) return sc;
    Event::CalDigiCol::const_iterator calDigiTds;

    for (calDigiTds = calDigiColTds->begin(); calDigiTds != calDigiColTds->end(); calDigiTds++) {
        CalDigi *calDigiRoot = new CalDigi();
        idents::CalXtalId::CalTrigMode modeTds = (*calDigiTds)->getMode();
        idents::CalXtalId idTds = (*calDigiTds)->getPackedId();
        CalXtalId idRoot(idTds.getTower(), idTds.getLayer(), idTds.getColumn());
        CalXtalId::CalTrigMode modeRoot;
        if (modeTds == idents::CalXtalId::BESTRANGE) {
            modeRoot = CalXtalId::BESTRANGE;
            calDigiRoot->initialize(modeRoot, idRoot);
            const Event::CalDigi::CalXtalReadout *readoutTds = (*calDigiTds)->getXtalReadout(0);
            Char_t rangePlusRoot = readoutTds->getRange(idents::CalXtalId::POS);
            UInt_t adcPlusRoot = readoutTds->getAdc(idents::CalXtalId::POS);
            Char_t rangeMinRoot = readoutTds->getRange(idents::CalXtalId::NEG);
            UInt_t adcMinRoot = readoutTds->getAdc(idents::CalXtalId::NEG);
            calDigiRoot->addReadout(rangePlusRoot, adcPlusRoot, rangeMinRoot, adcMinRoot);
        } else {
            modeRoot = CalXtalId::ALLRANGE;
            calDigiRoot->initialize(modeRoot, idRoot);
            int range;
            for (range = idents::CalXtalId::LEX8; range <= idents::CalXtalId::HEX1; range++) {
                const Event::CalDigi::CalXtalReadout *readoutTds = (*calDigiTds)->getXtalReadout(range);
                Char_t rangePlusRoot = readoutTds->getRange(idents::CalXtalId::POS);
                UInt_t adcPlusRoot = readoutTds->getAdc(idents::CalXtalId::POS);
                Char_t rangeMinRoot = readoutTds->getRange(idents::CalXtalId::NEG);
                UInt_t adcMinRoot = readoutTds->getAdc(idents::CalXtalId::NEG);
                calDigiRoot->addReadout(rangePlusRoot, adcPlusRoot, rangeMinRoot, adcMinRoot);
            }
        }

        m_digiEvt->addCalDigi(calDigiRoot);
    }

    return sc;
}

void digiRootWriterAlg::writeEvent() 
{
    // Purpose and Method:  Stores the DigiEvent data for this event in the ROOT
    //    tree.  The m_digiEvt object is cleared for the next event.

    TDirectory *saveDir = gDirectory;
    m_digiFile->cd();
    m_digiTree->Fill();
    m_digiEvt->Clear();
    saveDir->cd();
    return;
}

void digiRootWriterAlg::close() 
{
    // Purpose and Method:  Writes the ROOT file at the end of the run.
    //    The TObject::kOverWrite parameter is used in the Write method
    //    since ROOT will periodically write to the ROOT file when the bufSize
    //    is filled.  Writing would create 2 copies of the same tree to be
    //    stored in the ROOT file, if we did not specify kOverwrite.

    TDirectory *saveDir = gDirectory;
    m_digiFile->cd();
    m_digiFile->Write(0, TObject::kOverwrite);
    m_digiFile->Close();
    saveDir->cd();
    return;
}

StatusCode digiRootWriterAlg::finalize()
{
    close();
    
    StatusCode sc = StatusCode::SUCCESS;
    return sc;
}


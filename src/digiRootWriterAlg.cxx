#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/DigiEvent.h"
#include "Event/Digi/AcdDigi.h"
#include "Event/Digi/CalDigi.h"
#include "Event/Digi/TkrDigi.h"


#include "idents/CalXtalId.h"
#include "idents/TowerId.h"

#include "facilities/Util.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"

#include "digiRootData/DigiEvent.h"

#include "commonData.h"

#include "RootIo/IRootIoSvc.h"

/** @class digiRootWriterAlg
 * @brief Writes Digi TDS data to a persistent ROOT file.
 *
 * @author Heather Kelly
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/digiRootWriterAlg.cxx,v 1.18 2003/09/28 23:54:15 heather Exp $
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

    /// Retrieves ACD digitization data from the TDS and fill the AcdDigi
    /// ROOT collection
    StatusCode writeAcdDigi();

    /// Retrieves the CAL digitization data from the TDS and fills the CalDigi
    /// ROOT collection
    StatusCode writeCalDigi();

    /// Retrieves TKR digitization data from the TDS and fills the TkrDigi
    /// ROOT collection
    StatusCode writeTkrDigi();

    /// Calls TTree::Fill for each event and clears m_digiEvt
    void writeEvent();
    /// Converts from TDS VolId to ROOT VolId
    void convertVolumeId(idents::VolumeIdentifier tdsVolId, VolumeIdentifier& rootVolId) ;
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
    /// auto save every events
    int m_autoSaveEvents;
    
    commonData m_common;
    IRootIoSvc* m_rootIoSvc;

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
    declareProperty("autoSave", m_autoSaveEvents=1000);

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

    if ( service("RootIoSvc", m_rootIoSvc, true).isFailure() ){
        log << MSG::INFO << "Couldn't find the RootIoSvc!" << endreq;
        log << MSG::INFO << "No Auto Saving" << endreq;
        m_rootIoSvc = 0;
    } 

    facilities::Util::expandEnvVar(&m_fileName);

    // Save the current directory for the ntuple writer service
    TDirectory *saveDir = gDirectory;   
    // Create the new ROOT file
    m_digiFile = new TFile(m_fileName.c_str(), "RECREATE");
    if (!m_digiFile->IsOpen()) {
        log << MSG::ERROR << "ROOT file " << m_fileName 
            << " could not be opened for writing." << endreq;
        return StatusCode::FAILURE;
    }
    m_digiFile->cd();
    m_digiFile->SetCompressionLevel(m_compressionLevel);
    m_digiTree = new TTree(m_treeName.c_str(), "GLAST Digitization Data");
    m_digiEvt = new DigiEvent();
    m_common.m_digiEvt = m_digiEvt;
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
    
    if (!m_digiFile->IsOpen()) {
        log << MSG::ERROR << "ROOT file " << m_fileName 
            << " could not be opened for writing." << endreq;
        return StatusCode::FAILURE;
    }
    
    m_digiEvt->Clear();

    sc = writeDigiEvent();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to write DigiEvent" << endreq;
        return sc;
    }

    sc = writeAcdDigi();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to write Acd Digi Collection" << endreq;
        return sc;
    }

    sc = writeCalDigi();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to write Cal Digi Collection" << endreq;
        return sc;
    }

    sc = writeTkrDigi();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to write Tkr Digi Collection" << endreq;
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
    TimeStamp timeObj = evtTds->time();

    Bool_t fromMc = true;

    log << MSG::DEBUG;
    if( log.isActive()) evtTds->fillStream(log.stream());
    log << endreq;

    L1T levelOne(evtTds->trigger());
    m_digiEvt->initialize(evtId, runId, timeObj.time(), levelOne, fromMc);

    return sc;
}

StatusCode digiRootWriterAlg::writeAcdDigi() {
    // Purpose and Method:  Retrieve the AcdDigi collection from the TDS and 
    //    set the AcdDigi ROOT collection

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    SmartDataPtr<Event::AcdDigiCol> acdDigiColTds(eventSvc(), EventModel::Digi::AcdDigiCol);
    if (!acdDigiColTds) return sc;
    Event::AcdDigiCol::const_iterator acdDigiTds;
    
    for (acdDigiTds = acdDigiColTds->begin(); acdDigiTds != acdDigiColTds->end(); acdDigiTds++) {
        log << MSG::DEBUG ;
        if( log.isActive()) {
            (*acdDigiTds)->fillStream(log.stream());
        }
        log << endreq;
        Float_t energyRoot = (*acdDigiTds)->getEnergy();
        UShort_t phaRoot[2] = { (*acdDigiTds)->getPulseHeight(Event::AcdDigi::A),
            (*acdDigiTds)->getPulseHeight(Event::AcdDigi::B) };
        Bool_t vetoRoot[2] = { (*acdDigiTds)->getVeto(Event::AcdDigi::A),
            (*acdDigiTds)->getVeto(Event::AcdDigi::B) };
        Bool_t lowRoot[2] = { (*acdDigiTds)->getLowDiscrim(Event::AcdDigi::A),
            (*acdDigiTds)->getLowDiscrim(Event::AcdDigi::B) };
        Bool_t highRoot[2] = { (*acdDigiTds)->getHighDiscrim(Event::AcdDigi::A),
            (*acdDigiTds)->getHighDiscrim(Event::AcdDigi::B) };
        idents::AcdId idTds = (*acdDigiTds)->getId();
        AcdId idRoot;
        if (idTds.tile())
            idRoot.initialize(idTds.layer(), idTds.face(), idTds.row(), idTds.column());
        else
            idRoot.initialize(idTds.ribbonOrientation(), idTds.ribbonNum());

        const idents::VolumeIdentifier volIdTds = (*acdDigiTds)->getVolId();
        VolumeIdentifier volIdRoot;
        convertVolumeId(volIdTds, volIdRoot);

        m_digiEvt->addAcdDigi(idRoot, volIdRoot, energyRoot, phaRoot, 
            vetoRoot, lowRoot, highRoot);
    }

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
        log << MSG::DEBUG;
        if( log.isActive()) (*calDigiTds)->fillStream(log.stream());
        log << endreq;
        
        CalDigi *calDigiRoot = m_digiEvt->addCalDigi();
        m_common.m_calDigiMap[(*calDigiTds)] = calDigiRoot;

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

    }

    return sc;
}

StatusCode digiRootWriterAlg::writeTkrDigi() {
    // Purpose and Method:  Retrieve the TkrDigi collection from the TDS and set the
    //    TkrDigi ROOT collection

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    SmartDataPtr<Event::TkrDigiCol> tkrDigiColTds(eventSvc(), EventModel::Digi::TkrDigiCol);
    if (!tkrDigiColTds) return sc;
    Event::TkrDigiCol::const_iterator tkrDigiTds;

    for (tkrDigiTds = tkrDigiColTds->begin(); tkrDigiTds != tkrDigiColTds->end(); tkrDigiTds++) {
        idents::GlastAxis::axis axisTds = (*tkrDigiTds)->getView();
        GlastAxis::axis axisRoot = (axisTds == idents::GlastAxis::X) ? GlastAxis::X : GlastAxis::Y;
        idents::TowerId idTds = (*tkrDigiTds)->getTower();
        TowerId towerRoot(idTds.ix(), idTds.iy());
        Int_t totRoot[2] = {(*tkrDigiTds)->getToT(0), (*tkrDigiTds)->getToT(1)};
        Int_t lastController0Strip = (*tkrDigiTds)->getLastController0Strip();
       
        TkrDigi *tkrDigiRoot = new TkrDigi();
        m_common.m_tkrDigiMap[(*tkrDigiTds)] = tkrDigiRoot;

        tkrDigiRoot->initialize((*tkrDigiTds)->getBilayer(), axisRoot, towerRoot, totRoot);
        UInt_t numHits = (*tkrDigiTds)->getNumHits();
        unsigned int iHit;
        for (iHit = 0; iHit < numHits; iHit++) {
            Int_t strip = (*tkrDigiTds)->getHit(iHit);
            if (strip <= lastController0Strip) {
                tkrDigiRoot->addC0Hit(strip);
            } else {
                tkrDigiRoot->addC1Hit(strip);
            }
        }
        m_digiEvt->addTkrDigi(tkrDigiRoot);
    }

    return sc;
}

void digiRootWriterAlg::writeEvent() 
{
    // Purpose and Method:  Stores the DigiEvent data for this event in the ROOT
    //    tree.  The m_digiEvt object is cleared for the next event.
    static int eventCounter = 0;
    TDirectory *saveDir = gDirectory;
    m_digiFile->cd();
    m_digiTree->Fill();
    //m_digiEvt->Clear();
    saveDir->cd();
    ++eventCounter;
    if (m_rootIoSvc)
        if (eventCounter % m_rootIoSvc->getAutoSaveInterval() == 0) m_digiTree->AutoSave();

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

void digiRootWriterAlg::convertVolumeId(idents::VolumeIdentifier tdsVolId, 
                     VolumeIdentifier& rootVolId) 
{
    // Purpose and Method:  We must store the volume ids as two 32 bit UInt_t
    //     in the ROOT class.  Hence, we must convert the 64 bit representation
    //     used in the idents::VolumeIdentifier class into two 32 bit UInt_t.
    //     To perform the conversion, we iterate over all the ids in the TDS
    //     version of the idents::VolumeIdentifier and append each to the ROOT
    //     VolumeIdentifier
    
    int index;
    rootVolId.Clear();
    for (index = 0; index < tdsVolId.size(); index++) {
        rootVolId.append(tdsVolId.operator [](index));
    }
}

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/DigiEvent.h"
#include "Event/Digi/AcdDigi.h"
#include "Event/Digi/CalDigi.h"
#include "Event/Digi/TkrDigi.h"

#include "Event/RelTable/Relation.h"

#include "idents/CalXtalId.h"
#include "idents/TowerId.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TObjArray.h"
#include "TCollection.h"  // Declares TIter

#include "mcRootData/McEvent.h"
#include "digiRootData/DigiEvent.h"
#include "commonRootData/RelTable.h"

#include "facilities/Util.h"

#include "commonData.h"

/** @class relationRootReaderAlg
 * @brief Reads Digitization data from a persistent ROOT file and stores the
 * the data in the TDS.
 *
 * @author Heather Kelly
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/relationRootReaderAlg.cxx,v 1.7 2002/06/20 21:36:31 heather Exp $
 */

class relationRootReaderAlg : public Algorithm
{	
public:
    
    relationRootReaderAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    /// Handles setup by opening ROOT file in read mode and creating a new TTree
    StatusCode initialize();
   
    /// Orchastrates reading from ROOT file and storing the data on the TDS for each event
    StatusCode execute();
    
    /// Closes the ROOT file and cleans up
    StatusCode finalize();
            
private:

    /// Reads top-level DigiEvent
    StatusCode readRelations();
    
    StatusCode readTkrDigiRelations();
  
    StatusCode readCalDigiRelations();

    /// Closes the ROOT file
    void close();
   
    /// ROOT file pointer
    TFile *m_relFile;
    /// ROOT tree pointer
    TTree *m_relTree;
    /// Top-level Monte Carlo ROOT object
    RelTable *m_relTab;
    /// name of the output ROOT file
    std::string m_fileName;
    /// name of the Monte Carlo TTree stored in the ROOT file
    std::string m_treeName;
    /// Stores number of events available in the input ROOT TTree
    int m_numEvents;

};

static const AlgFactory<relationRootReaderAlg>  Factory;
const IAlgFactory& relationRootReaderAlgFactory = Factory;


relationRootReaderAlg::relationRootReaderAlg(const std::string& name, ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator)
{
    // Input pararmeters that may be set via the jobOptions file
    // Input ROOT file name
    declareProperty("rootFile",m_fileName="relations.root");
    // Input TTree name
    declareProperty("treeName", m_treeName="Relations");

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

    facilities::Util::expandEnvVar(&m_fileName);
    
    // Save the current directory for the ntuple writer service
    TDirectory *saveDir = gDirectory;   
    m_relFile = new TFile(m_fileName.c_str(), "READ");
    if (!m_relFile->IsOpen()) {
        log << MSG::ERROR << "ROOT file " << m_fileName 
            << " could not be opened for reading." << endreq;
        return StatusCode::FAILURE;
    }
    m_relFile->cd();
    m_relTree = (TTree*)m_relFile->Get(m_treeName.c_str());
    if (!m_relTree) {
        log << MSG::ERROR << "Could not load Tree " << m_treeName << 
            " from file " << m_fileName << endreq;
        return StatusCode::FAILURE;
    }
    m_relTab = 0;
    m_relTree->SetBranchAddress("RelTable", &m_relTab);

    m_numEvents = m_relTree->GetEntries();
    
    saveDir->cd();
    return sc;
    
}

StatusCode relationRootReaderAlg::execute()
{
    // Purpose and Method:  Called once per event.  This method calls
    //   the appropriate methods to read data from the ROOT file and store
    //   data on the TDS.

    MsgStream log(msgSvc(), name());

    StatusCode sc = StatusCode::SUCCESS;
    
    if (!m_relFile->IsOpen()) {
        log << MSG::ERROR << "ROOT file " << m_fileName 
            << " could not be opened for reading." << endreq;
        return StatusCode::FAILURE;
    }

    static UInt_t evtId = 0;

    if (evtId >= m_numEvents) {
        log << MSG::ERROR << "ROOT file contains no more events" << endreq;
        return StatusCode::FAILURE;
    }

    m_relTree->GetEvent(evtId);

    sc = readRelations();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to read Relational Table" << endreq;
        return sc;
    }


    m_relTab->Clear();
    evtId++;
    
    return sc;
}


StatusCode relationRootReaderAlg::readRelations() {

    MsgStream log(msgSvc(), name());

    StatusCode sc = StatusCode::SUCCESS;
    
    return sc;
}

StatusCode relationRootReaderAlg::readTkrDigiRelations() {

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    

    return sc;
}

StatusCode relationRootReaderAlg::readCalDigiRelations() {

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    

    return sc;
}


void relationRootReaderAlg::close() 
{
    // Purpose and Method:  Writes the ROOT file at the end of the run.
    //    The TObject::kOverWrite parameter is used in the Write method
    //    since ROOT will periodically write to the ROOT file when the bufSize
    //    is filled.  Writing would create 2 copies of the same tree to be
    //    stored in the ROOT file, if we did not specify kOverwrite.

    TDirectory *saveDir = gDirectory;
    m_relFile->cd();
    m_relFile->Close();
    saveDir->cd();
}

StatusCode relationRootReaderAlg::finalize()
{
    close();
    
    StatusCode sc = StatusCode::SUCCESS;
    return sc;
}


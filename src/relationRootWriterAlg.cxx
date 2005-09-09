#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

//#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/DigiEvent.h"
#include "Event/Digi/AcdDigi.h"
#include "Event/Digi/CalDigi.h"
#include "Event/Digi/TkrDigi.h"

#include "Event/RelTable/RelTable.h"

#include "idents/CalXtalId.h"
#include "idents/TowerId.h"

#include "facilities/Util.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TRef.h"
#include "TRefArray.h"

#include "commonRootData/RelTable.h"

#include "commonData.h"
#include "RootIo/IRootIoSvc.h"

#include "digiRootData/DigiEvent.h"
#include "mcRootData/McEvent.h"

#include <map>

/** @class relationRootWriterAlg
 * @brief Writes relational table TDS data to a persistent ROOT file.
 *
 * @author Heather Kelly
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/relationRootWriterAlg.cxx,v 1.12 2004/12/16 00:24:29 usher Exp $
 */

class relationRootWriterAlg : public Algorithm
{	
public:
    
    relationRootWriterAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    /// Handles setup by opening ROOT file in write mode and creating a new TTree
    StatusCode initialize();
   
    /// Orchastrates reading from TDS and writing to ROOT for each event
    StatusCode execute();
    
    /// Closes the ROOT file and cleans up
    StatusCode finalize();

private:
    /// Fills the TkrDigi/McPositionHits relational table
    StatusCode writeTkrDigiRelations();
    /// Fills the CalDigi/McIntegratingHits relational table
    StatusCode writeCalDigiRelations();

    /// Typedefs for map between root objects and TRefArrays
    typedef std::map<TObject*, TRefArray> relationMap;
    typedef relationMap::const_iterator relationMapIt;

    /// Standard method for filling RelTable
    void fillRelTable(const relationMap& relMap);

    /// Calls TTree::Fill for each event and clears m_digiEvt
    void writeEvent();
    /// Performs the final write to the ROOT file and closes
    void close();
   
    /// ROOT file pointer
    TFile *m_relFile;
    /// ROOT tree pointer
    TTree *m_relTree;
    /// Top-level ROOT object
    RelTable *m_relTable;
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
    
    commonData m_common;
    IRootIoSvc* m_rootIoSvc;
};


static const AlgFactory<relationRootWriterAlg>  Factory;
const IAlgFactory& relationRootWriterAlgFactory = Factory;

relationRootWriterAlg::relationRootWriterAlg(const std::string& name, 
                                 ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator)
{
    // Input parameters available to be set via the jobOptions file
    declareProperty("rootFile",m_fileName="relations.root");
    declareProperty("splitMode", m_splitMode=1);
    declareProperty("bufferSize", m_bufSize=64000);
    // ROOT default compression
    declareProperty("compressionLevel", m_compressionLevel=1);
    declareProperty("treeName", m_treeName="Relations");
}

StatusCode relationRootWriterAlg::initialize()
{
    // Purpose and Method:  Called once before the run begins.  This method
    //    opens a new ROOT file and prepares for writing.

    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    // Use the Job options service to set the Algorithm's parameters
    // This will retrieve parameters set in the job options file
    setProperties();

    m_rootIoSvc = 0 ;
    if ( service("RootIoSvc", m_rootIoSvc, true).isFailure() ){
        log << MSG::INFO << "Couldn't find the RootIoSvc!" << endreq;
        log << MSG::INFO << "No Auto Saving" << endreq;
        m_rootIoSvc = 0;
    } 
    
    facilities::Util::expandEnvVar(&m_fileName);

    // Save the current directory for the ntuple writer service
    TDirectory *saveDir = gDirectory;   
    // Create the new ROOT file
    m_relFile = new TFile(m_fileName.c_str(), "RECREATE");
    if (!m_relFile->IsOpen()) {
        log << MSG::ERROR << "ROOT file " << m_fileName 
            << " could not be opened for writing." << endreq;
        return StatusCode::FAILURE;
    }
    m_relFile->cd();
    m_relFile->SetCompressionLevel(m_compressionLevel);
    m_relTree = new TTree(m_treeName.c_str(), "GLAST Relational Table");
    m_relTable = new RelTable();
    m_relTree->Branch("RelTable","RelTable", &m_relTable, m_bufSize, m_splitMode);
    
    saveDir->cd();

    return sc;
    
}

StatusCode relationRootWriterAlg::execute()
{
    // Purpose and Method:  Called once per event.  This method calls
    //   the appropriate methods to read data from the TDS and write data
    //   to the ROOT file.

    MsgStream log(msgSvc(), name());

    StatusCode sc = StatusCode::SUCCESS;
    
    if (!m_relTree->GetCurrentFile()->IsOpen()) {
        log << MSG::ERROR << "ROOT file " << m_fileName 
            << " could not be opened for writing." << endreq;
        return StatusCode::FAILURE;
    }

    sc = writeTkrDigiRelations();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to write TkrDigiRelations" << endreq;
        return sc;
    }

    sc = writeCalDigiRelations();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to write CalDigiRelations" << endreq;
        return sc;
    }

    writeEvent();
    return sc;
}


StatusCode relationRootWriterAlg::writeTkrDigiRelations() {
    // Purpose and Method:  Retrieve the relations concerning TkrDigi from TDS
    //    Write the relations to ROOT

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    typedef Event::Relation<Event::TkrDigi,Event::McPositionHit> relType;
    typedef ObjectList<relType> tabType;

    SmartDataPtr<tabType> tkrTable(eventSvc(), EventModel::Digi::TkrDigiHitTab);

    relationMap tkrRelationMap;

    tabType::const_iterator relation;
    for (relation = tkrTable->begin(); relation != tkrTable->end(); relation++) {
        Event::TkrDigi *tkrDigiTds = (*relation)->getFirst();
        Event::McPositionHit *posHitTds = (*relation)->getSecond();

        TkrDigi *tkrDigiRoot = 0;
        McPositionHit *posHitRoot = 0;

        if (m_common.m_mcPosHitMap.find(posHitTds) != m_common.m_mcPosHitMap.end()) {
            TRef ref = m_common.m_mcPosHitMap[posHitTds];
            posHitRoot = (McPositionHit*)ref.GetObject();
        } else {
            log << MSG::WARNING << "Could not located McPositionHit TDS/ROOT pair" << endreq;
        }

        if (m_common.m_tkrDigiMap.find(tkrDigiTds) != m_common.m_tkrDigiMap.end()) {
            TRef ref = m_common.m_tkrDigiMap[tkrDigiTds];
            tkrDigiRoot = (TkrDigi*)ref.GetObject();
        } else {
            log << MSG::WARNING << "Could not located TkrDigi TDS/ROOT pair" << endreq;
        }

        if (tkrDigiRoot && posHitRoot) tkrRelationMap[tkrDigiRoot].Add(posHitRoot);


    }

   fillRelTable(tkrRelationMap);

    return sc;
}

StatusCode relationRootWriterAlg::writeCalDigiRelations() {
    // Purpose and Method:  Retrieve the relations concerning CalDigi from TDS
    //    Write the relations to ROOT

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    typedef Event::Relation<Event::CalDigi,Event::McIntegratingHit> relType;
    typedef ObjectList<relType> tabType;

    SmartDataPtr<tabType> calTable(eventSvc(), EventModel::Digi::CalDigiHitTab);
    tabType::const_iterator relation;
    relationMap calRelationMap;

    for (relation = calTable->begin(); relation != calTable->end(); relation++) {
        Event::CalDigi *calDigiTds = (*relation)->getFirst();
        Event::McIntegratingHit *intHitTds = (*relation)->getSecond();

        CalDigi *calDigiRoot = 0;
        McIntegratingHit *intHitRoot = 0;
        if (m_common.m_mcIntHitMap.find(intHitTds) != m_common.m_mcIntHitMap.end()) {
            TRef ref = m_common.m_mcIntHitMap[intHitTds];
            intHitRoot = (McIntegratingHit*)ref.GetObject();
        } else {
            log << MSG::WARNING << "Could not located McIntegratingHit TDS/ROOT pair" << endreq;
        }

        if (m_common.m_calDigiMap.find(calDigiTds) != m_common.m_calDigiMap.end()) {
            TRef ref = m_common.m_calDigiMap[calDigiTds];
            calDigiRoot = (CalDigi*)ref.GetObject();
        } else {
            log << MSG::WARNING << "Could not located CalDigi TDS/ROOT pair" << endreq;
        }

        if (calDigiRoot && intHitRoot) calRelationMap[calDigiRoot].Add(intHitRoot);

    }

    fillRelTable(calRelationMap);

    return sc;
}

void relationRootWriterAlg::fillRelTable(const relationMap& relMap)
{
    // Purpose and Method:  Centralizes the filling of root relations into the table
    //    Write the relations to ROOT
    relationMapIt mapIt;
    for (mapIt = relMap.begin(); mapIt != relMap.end(); mapIt++) 
    {
        TRef       tkrObj = (*mapIt).first;
        TRefArray  tkrArr = (*mapIt).second;
        Relation*  rel    = new Relation(tkrObj, tkrArr);
        m_relTable->addRelation(rel);
    }

    return;
}

void relationRootWriterAlg::writeEvent() 
{
    // Purpose and Method:  Stores the Relations data for this event in the ROOT
    //    tree.  The m_common object is cleared for the next event.

    static int eventCounter = 0;
    try {
        TDirectory *saveDir = gDirectory;
        m_relTree->GetCurrentFile()->cd();
        m_relTree->Fill();
        m_relTable->Clear();
        m_common.clear();
        ++eventCounter;
        if (m_rootIoSvc)
            if (eventCounter % m_rootIoSvc->getAutoSaveInterval() == 0) 
                m_relTree->AutoSave();
        saveDir->cd();
      } catch(...) { 
          std::cerr << "Failed to write the event to file" << std::endl; 
          std::cerr << "Exiting..." << std::endl; 
          std::cerr.flush(); 
          exit(1); 
      } 


    return;
}

void relationRootWriterAlg::close() 
{
    // Purpose and Method:  Writes the ROOT file at the end of the run.
    //    The TObject::kWriteDelete parameter is used in the Write method
    //    Used rather than TObject::kOverwrite - supposed to be safer but slower
    //    since ROOT will periodically write to the ROOT file when the bufSize
    //    is filled.  Writing would create 2 copies of the same tree to be
    //    stored in the ROOT file, if we did not specify kOverwrite.

    try {
        TDirectory *saveDir = gDirectory;
        TFile *f = m_relTree->GetCurrentFile();
        f->cd();
        m_relTree->BuildIndex("m_runId", "m_eventId");
        f->Write(0, TObject::kWriteDelete);
        f->Close();
        saveDir->cd();
     } catch(...) { 
          std::cerr << "Failed final write to RELATION file" << std::endl; 
          std::cerr << "Exiting..." << std::endl; 
          std::cerr.flush(); 
          exit(1); 
     } 

    return;
}

StatusCode relationRootWriterAlg::finalize()
{
    close();
    
    StatusCode sc = StatusCode::SUCCESS;
    setFinalized();
    return sc;
}


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
#include "Event/Recon/TkrRecon/TkrTrackTab.h"
#include "Event/Recon/TkrRecon/TkrVertexTab.h"

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
#include "reconRootData/ReconEvent.h"
#include "commonRootData/RelTable.h"

#include "facilities/Util.h"

#include "commonData.h"

/** @class relationRootReaderAlg
 * @brief Reads Digitization data from a persistent ROOT file and stores the
 * the data in the TDS.
 *
 * @author Heather Kelly
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/relationRootReaderAlg.cxx,v 1.1 2002/12/02 21:54:19 heather Exp $
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

    /// Creates TDS tables
    StatusCode createTDSTables();

    /// Reads top-level DigiEvent
    StatusCode readRelations();
    

    /// Handles converting TkrDigi relations back to TDS
    StatusCode readTkrDigiRelations(Relation*);
    /// Handles converting CalDigi relations back to TDS
    StatusCode readCalDigiRelations(Relation*);
    /// Handles converting track relations back to TDS
    StatusCode readTkrTrackRelations(Relation*);
    /// Handles converting vertex relations back to TDS
    StatusCode readTkrVertexRelations(Relation*);

    //StatusCode readTkrDigiRelations();  
    //StatusCode readCalDigiRelations();

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

    commonData m_common;

    /// typedefs for tables
    typedef Event::RelTable<Event::TkrDigi,Event::McPositionHit>    TkrDigiRelTab;
    typedef Event::Relation<Event::TkrDigi,Event::McPositionHit>    TkrDigiRelType;
    typedef Event::RelTable<Event::CalDigi,Event::McIntegratingHit> CalDigiRelTab;
    typedef Event::Relation<Event::CalDigi,Event::McIntegratingHit> CalDigiRelType;

    /// Internal pointer to TDS TkrDigi table
    TkrDigiRelTab*         m_tkrDigiRelTab;
    /// Internal pointer to tDS CalDigi table
    CalDigiRelTab*         m_calDigiRelTab;
    /// Internal pointer to TDS TkrTrackTab relational table
    Event::TkrFitTrackTab* m_trackRelTab;
    /// Internal pointer to TDS TkrVertexTab relational table
    Event::TkrVertexTab*   m_vertexRelTab;


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


    sc = createTDSTables();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to create TDS relational tables" << endreq;
        return sc;
    }

    sc = readRelations();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to read Relational Table" << endreq;
        return sc;
    }


    /// Register our relational tables in the TDS (and turn over ownership)
    sc = eventSvc()->registerObject(EventModel::Digi::TkrDigiHitTab,    m_tkrDigiRelTab->getAllRelations() );
    sc = eventSvc()->registerObject(EventModel::Digi::CalDigiHitTab,    m_calDigiRelTab->getAllRelations() );
    sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrTrackTab,  m_trackRelTab->getAllRelations()   );
    sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrVertexTab, m_vertexRelTab->getAllRelations()  );


    m_relTab->Clear();
    evtId++;

    m_common.clear();
    
    return sc;
}

StatusCode relationRootReaderAlg::createTDSTables()
{
    // Purpose and Method:  Creates the TDS relational tables for output
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    m_tkrDigiRelTab = new TkrDigiRelTab();
    m_tkrDigiRelTab->init();

    m_calDigiRelTab = new CalDigiRelTab();
    m_calDigiRelTab->init();

    m_trackRelTab = new Event::TkrFitTrackTab();
    m_trackRelTab->init();

    m_vertexRelTab = new Event::TkrVertexTab();
    m_vertexRelTab->init();

    return sc;
}

StatusCode relationRootReaderAlg::readRelations() {
    // Purpose and Method:  "Reads" the individual relations and, depending upon type,
    //                      calls the appropriate method for converting back to TDS relation
    MsgStream log(msgSvc(), name());

    StatusCode sc = StatusCode::SUCCESS;

    // Get pointer to the root relations
    const TObjArray* relations = m_relTab->getRelationTable();

    // Loop over these relations
    for(int idx = 0; idx < relations->GetEntries(); idx++)
    {
        TObject*  relObj   = relations->At(idx);
        Relation* relation = dynamic_cast<Relation*>(relObj);

        const TObject* key = relation->getKey();

        // Dynamic cast the key back to its original object and decide how to read the relation
        if      (const TkrDigi*      digi = dynamic_cast<const TkrDigi*>(key))      readTkrDigiRelations(relation);
        else if (const CalDigi*      digi = dynamic_cast<const CalDigi*>(key))      readCalDigiRelations(relation);
        else if (const TkrCandTrack* cand = dynamic_cast<const TkrCandTrack*>(key)) readTkrTrackRelations(relation);
        else if (const TkrVertex*    vert = dynamic_cast<const TkrVertex*>(key))    readTkrVertexRelations(relation);
        else
        {
            log << MSG::WARNING << "Unrecognized Relation key" << endreq;
        }
    }
    
    return sc;
}

StatusCode relationRootReaderAlg::readTkrDigiRelations(Relation* relation) 
{
    // Purpose and Method:  Converts a root TkrDigi/McPositionHit relation to TDS version
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    const TkrDigi*        digiRoot  = dynamic_cast<const TkrDigi*>(relation->getKey());
    TObject*              mcHitRoot = relation->getValueCol().At(0);

    Event::TkrDigi*       digiTds   = 0;
    Event::McPositionHit* mcHitTds  = 0;

    // Look up TDS to root relation for digis
    if (m_common.m_rootTkrDigiMap.find(digiRoot) != m_common.m_rootTkrDigiMap.end()) {
        digiTds = const_cast<Event::TkrDigi*>(m_common.m_rootTkrDigiMap[digiRoot]);
    } else {
        log << MSG::WARNING << "Could not located TkrDigi TDS/ROOT pair" << endreq;
    }

    // Look up TDS to root relation for position hits
    if (m_common.m_rootMcPosHitMap.find(mcHitRoot) != m_common.m_rootMcPosHitMap.end()) {
        mcHitTds = const_cast<Event::McPositionHit*>(m_common.m_rootMcPosHitMap[mcHitRoot]);
    } else {
        log << MSG::WARNING << "Could not located McPositionHit TDS/ROOT pair" << endreq;
    }

    // Make the event relation here
    typedef Event::Relation<Event::TkrDigi,Event::McPositionHit> relType;
    relType* rel = new relType(digiTds, mcHitTds);
    m_tkrDigiRelTab->addRelation(rel);

    return sc;
}


StatusCode relationRootReaderAlg::readCalDigiRelations(Relation* relation) 
{
    // Purpose and Method:  Converts a root CalDigi/McIntegratingHit relation to TDS version

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    const CalDigi*           digiRoot  = dynamic_cast<const CalDigi*>(relation->getKey());
    TObject*                 mcHitRoot = relation->getValueCol().At(0);

    Event::CalDigi*          digiTds   = 0;
    Event::McIntegratingHit* mcHitTds  = 0;

    // Look up TDS to root relation for digis
    if (m_common.m_rootCalDigiMap.find(digiRoot) != m_common.m_rootCalDigiMap.end()) {
        digiTds = const_cast<Event::CalDigi*>(m_common.m_rootCalDigiMap[digiRoot]);
    } else {
        log << MSG::WARNING << "Could not located CalDigi TDS/ROOT pair" << endreq;
    }

    // Look up TDS to root relation for position hits
    if (m_common.m_rootMcIntHitMap.find(mcHitRoot) != m_common.m_rootMcIntHitMap.end()) {
        mcHitTds = const_cast<Event::McIntegratingHit*>(m_common.m_rootMcIntHitMap[mcHitRoot]);
    } else {
        log << MSG::WARNING << "Could not located McIntegratingHit TDS/ROOT pair" << endreq;
    }

    // Make the event relation here
    typedef Event::Relation<Event::CalDigi,Event::McIntegratingHit> relType;
    relType* rel = new relType(digiTds, mcHitTds);
    m_calDigiRelTab->addRelation(rel);

    return sc;
}


StatusCode relationRootReaderAlg::readTkrTrackRelations(Relation* relation) 
{
    // Purpose and Method:  Converts a root TkrPatCand/TkrFitTrack relation to TDS version

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    const TkrCandTrack*     tkrCandRoot  = dynamic_cast<const TkrCandTrack*>(relation->getKey());

    const TRefArray&        refArray     = relation->getValueCol();
    TObject*                tkrTrackRoot = refArray.At(0);
    Event::TkrPatCand*      tkrCandTds   = 0;
    Event::TkrFitTrackBase* tkrTrackTds  = 0;

    // Look up TDS to root relation for candidate tracks
    if (m_common.m_rootTkrCandMap.find(tkrCandRoot) != m_common.m_rootTkrCandMap.end()) {
        tkrCandTds = const_cast<Event::TkrPatCand*>(m_common.m_rootTkrCandMap[tkrCandRoot]);
    } else {
        log << MSG::WARNING << "Could not located TkrPatCand TDS/ROOT pair" << endreq;
    }

    // Look up TDS to root relation for candidate tracks
    if (m_common.m_rootTkrTrackMap.find(tkrTrackRoot) != m_common.m_rootTkrTrackMap.end()) {
        tkrTrackTds = const_cast<Event::TkrFitTrackBase*>(m_common.m_rootTkrTrackMap[tkrTrackRoot]);
    } else {
        log << MSG::WARNING << "Could not located TkrFitTrackBase TDS/ROOT pair" << endreq;
    }

    // Make the event relation here
    Event::TkrFitTrackRel* rel = new Event::TkrFitTrackRel(tkrCandTds, tkrTrackTds);
    m_trackRelTab->addRelation(rel);
    
    return sc;
}

StatusCode relationRootReaderAlg::readTkrVertexRelations(Relation* relation) 
{
    // Purpose and Method:  Converts a root TkrVertex/TkrFitTrack relation to TDS version

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    const TkrVertex*       vertexRoot    = dynamic_cast<const TkrVertex*>(relation->getKey());

    const TRefArray&        refArray     = relation->getValueCol();
    TObject*                tkrTrackRoot = refArray.At(0);
    Event::TkrVertex*       vertexTds    = 0;
    Event::TkrFitTrackBase* tkrTrackTds  = 0;

    // Look up TDS to root relation for vertices
    if (m_common.m_rootTkrVertexMap.find(vertexRoot) != m_common.m_rootTkrVertexMap.end()) {
        vertexTds = const_cast<Event::TkrVertex*>(m_common.m_rootTkrVertexMap[vertexRoot]);
    } else {
        log << MSG::WARNING << "Could not located TkrVertex TDS/ROOT pair" << endreq;
    }

    // Look up TDS to root relation for candidate tracks
    if (m_common.m_rootTkrTrackMap.find(tkrTrackRoot) != m_common.m_rootTkrTrackMap.end()) {
        tkrTrackTds = const_cast<Event::TkrFitTrackBase*>(m_common.m_rootTkrTrackMap[tkrTrackRoot]);
    } else {
        log << MSG::WARNING << "Could not located TkrFitTrackBase TDS/ROOT pair" << endreq;
    }

    // Make the event relation here
    Event::TkrVertexRel* rel = new Event::TkrVertexRel(vertexTds, tkrTrackTds);
    m_vertexRelTab->addRelation(rel);

    return sc;
}



void relationRootReaderAlg::close() 
{
    // Purpose and Method:  Closes the ROOT file at the end of the run.

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


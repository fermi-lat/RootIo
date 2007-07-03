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
#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McTrajectory.h"
#include "Event/MonteCarlo/McRelTableDefs.h"

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
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/relationRootWriterAlg.cxx,v 1.16 2007/02/21 01:57:32 usher Exp $
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
    /// Fills the Monte Carlo specific relational tables
    StatusCode writeMcRelations();

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

    // Get the run and event info first (used for event indexing)
    SmartDataPtr<Event::EventHeader> evt(eventSvc(), EventModel::EventHeader);
    
    if (!evt) return sc;

    UInt_t evtId = evt->event();
    UInt_t runId = evt->run();

    m_relTable->initialize(evtId, runId);


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

    sc = writeMcRelations();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to write McRelations" << endreq;
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
        // Note that noise hits will not have any McPositionHit associated with them
        } else {
        //    log << MSG::WARNING << "Could not located McPositionHit TDS/ROOT pair" << endreq;
            continue;
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
        // It can happen that a noise hit will not have an McIntegratingHit associated
        } else {
            //log << MSG::WARNING << "Could not located McIntegratingHit TDS/ROOT pair" << endreq;
            continue;
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

StatusCode relationRootWriterAlg::writeMcRelations() {
    // Purpose and Method:  Retrieve the relations concerning the Monte Carlo
    //                      from the TDS, convert them and write to ROOT output
    // 
    // There are three sets of MC relations we are concerned with here:
    // 1) McParticle <--> McTrajectory
    // 2) McTrajectoryPoint <--> McPositionHit
    // 3) McTrajectoryPoint <--> McIntegratingHit
    // 

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Root map for intermediate conversion results
    relationMap mcRelationMap;
    relationMap mcPosHitMap;
    relationMap mcIntHitMap;

    // Retrieve the list of relations between McParticles and McTrajectory's
    // which should be one to one
    SmartDataPtr<Event::McPartToTrajectoryTabList> 
        mcPartToTraj(eventSvc(), "/Event/MC/McPartToTrajectory");

    // Recover the McPositionHits to trajectory points relational tables
    SmartDataPtr<Event::McPointToPosHitTabList> 
        mcPointToPosHitList(eventSvc(), "/Event/MC/McPointToPosHit");
    Event::McPointToPosHitTab mcPosHitTab(mcPointToPosHitList);

    // Recover the McPositionHits to trajectory points relational tables
    SmartDataPtr<Event::McPointToIntHitTabList> 
        mcPointToIntHitList(eventSvc(), "/Event/MC/McPointToIntHit");
    Event::McPointToIntHitTab mcIntHitTab(mcPointToIntHitList);


    // Look up McParticle <--> McTrajectory relational table
    if (mcPartToTraj != 0)
    {
        // Loop through this list at the top level
        for(Event::McPartToTrajectoryTabList::const_iterator trajIter = mcPartToTraj->begin(); 
            trajIter != mcPartToTraj->end(); trajIter++) 
        {
            Event::McPartToTrajectoryRel* partToTraj = *trajIter;

            const Event::McParticle*   mcPartTds = partToTraj->getFirst();
            const Event::McTrajectory* mcTrajTds = partToTraj->getSecond();

            McParticle*   mcPartRoot = 0;
            McTrajectory* mcTrajRoot = 0;
            if (m_common.m_mcPartMap.find(mcPartTds) != m_common.m_mcPartMap.end()) {
                TRef ref = m_common.m_mcPartMap[mcPartTds];
                mcPartRoot = (McParticle*)ref.GetObject();
            } else {
                log << MSG::WARNING << "Could not located McParticle TDS/ROOT pair" << endreq;
            }

            if (m_common.m_mcTrajectoryMap.find(mcTrajTds) != m_common.m_mcTrajectoryMap.end()) {
                TRef ref = m_common.m_mcTrajectoryMap[mcTrajTds];
                mcTrajRoot = (McTrajectory*)ref.GetObject();
            } else {
                log << MSG::WARNING << "Could not located McTrajectory TDS/ROOT pair" << endreq;
            }

            // Pointers must be found at this level to proceed
            if (mcPartRoot && mcTrajRoot)
            {
                mcRelationMap[mcPartRoot].Add(mcTrajRoot);

                // Set up to loop through Trajectory points
                std::vector<Event::McTrajectoryPoint*> points = mcTrajTds->getPoints();
                std::vector<Event::McTrajectoryPoint*>::const_iterator pointIter;

                for(pointIter = points.begin(); pointIter != points.end(); pointIter++) 
                {
                    // De-reference the McTrajectoryPoint pointer
                    Event::McTrajectoryPoint* mcPointTds  = *pointIter;
                    McTrajectoryPoint*        mcPointRoot = 0;

                    if (m_common.m_mcTrajectoryPointMap.find(mcPointTds) != m_common.m_mcTrajectoryPointMap.end()) {
                        TRef ref = m_common.m_mcTrajectoryPointMap[mcPointTds];
                        mcPointRoot = (McTrajectoryPoint*)ref.GetObject();
                    } else {
                        log << MSG::WARNING << "Could not located McTrajectoryPoint TDS/ROOT pair" << endreq;
                        continue;
                    }

                    // First look for McTrajectoryPoint <--> McPositionHit
                    Event::McPointToPosHitVec posRelVec = mcPosHitTab.getRelByFirst(mcPointTds);

                    if (!posRelVec.empty())
                    {
                        Event::McPointToPosHitRel* hitRel = *(posRelVec.begin());

                        Event::McPositionHit* mcPosHitTds  = hitRel->getSecond();
                        McPositionHit*        mcPosHitRoot = 0;

                        if (m_common.m_mcPosHitMap.find(mcPosHitTds) != m_common.m_mcPosHitMap.end()) {
                            TRef ref = m_common.m_mcPosHitMap[mcPosHitTds];
                            mcPosHitRoot = (McPositionHit*)ref.GetObject();
                        } else {
                            log << MSG::WARNING << "Could not located McPositionHit TDS/ROOT pair" << endreq;
                        }

                        if (mcPosHitRoot) mcPosHitMap[mcPointRoot].Add(mcPosHitRoot);
                    }

                    // Now look for McTrajectoryPoint <--> McIntegratingHit
                    Event::McPointToIntHitVec intRelVec = mcIntHitTab.getRelByFirst(mcPointTds);

                    if (!intRelVec.empty())
                    {
                        Event::McPointToIntHitRel* hitRel = *(intRelVec.begin());

                        Event::McIntegratingHit* mcIntHitTds  = hitRel->getSecond();
                        McIntegratingHit*        mcIntHitRoot = 0;

                        if (m_common.m_mcIntHitMap.find(mcIntHitTds) != m_common.m_mcIntHitMap.end()) {
                            TRef ref = m_common.m_mcIntHitMap[mcIntHitTds];
                            mcIntHitRoot = (McIntegratingHit*)ref.GetObject();
                        } else {
                            log << MSG::WARNING << "Could not located McIntegratingHit TDS/ROOT pair" << endreq;
                        }

                        if (mcIntHitRoot) mcIntHitMap[mcPointRoot].Add(mcIntHitRoot);
                    }
                }
            }
        }
    }

    fillRelTable(mcRelationMap);
    fillRelTable(mcPosHitMap);
    fillRelTable(mcIntHitMap);

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
    if (m_relTree->GetCurrentFile()->TestBits(TFile::kWriteError)) {
        throw;
    }
    m_relTree->Fill();
    const_cast<TObjArray*>(m_relTable->getRelationTable())->SetOwner(true);
    m_relTable->Clear();
    m_common.clear();
    ++eventCounter;
    if (m_rootIoSvc)
        if (eventCounter % m_rootIoSvc->getAutoSaveInterval() == 0) 
            if (m_relTree->AutoSave() == 0) throw;
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


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
#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McTrajectory.h"
#include "Event/MonteCarlo/McRelTableDefs.h"

#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalClusterTab.h"

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
#include "reconRootData/ReconEvent.h"

#include <map>

/** @class relationRootWriterAlg
 * @brief Writes relational table TDS data to a persistent ROOT file.
 * Note that this algorithm should be run AFTER all other ROOT writing algorithms, so that the objects pointed to in
 * the relation table exist when the table is written.
 *
 * @author Heather Kelly
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/relationRootWriterAlg.cxx,v 1.27 2010/11/03 22:13:42 usher Exp $
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
    /// Fills the CalRecon specific relational tables
    StatusCode writeCalReconRelations();

    /// Standard method for filling a root relation
    void addRootRelation(TRef& first, TRef& second, std::vector<std::string> infos);

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
        log << MSG::WARNING << "Couldn't find the RootIoSvc!" << endreq;
        m_rootIoSvc = 0;
        // RootIoSvc is now required for writing/reading, we cannot continue without it
        return StatusCode::FAILURE;
    }

    m_relTree = m_rootIoSvc->prepareRootOutput("rel", m_fileName, m_treeName, 
        m_compressionLevel, "GLAST Digitization Data");

    m_relTable = new RelTable();
    m_rootIoSvc->setupBranch("rel", "RelTable", "RelTable", &m_relTable, m_bufSize, m_splitMode);
    
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

    sc = writeCalReconRelations();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to write CalReconRelations" << endreq;
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
    // When processing data, no MC is available
    if (!tkrTable) return sc;
    
    tabType::const_iterator relation;
    for (relation = tkrTable->begin(); relation != tkrTable->end(); relation++) {
        Event::TkrDigi *tkrDigiTds = (*relation)->getFirst();
        Event::McPositionHit *posHitTds = (*relation)->getSecond();

        TkrDigi *tkrDigiRoot = 0;
        TRef     tkrDigiRef  = 0;
        TRef     posHitRef   = 0;
        McPositionHit *posHitRoot = 0;

        // Look up the key first
        if (m_common.m_tkrDigiMap.find(tkrDigiTds) != m_common.m_tkrDigiMap.end()) {
            tkrDigiRef  = m_common.m_tkrDigiMap[tkrDigiTds];
            tkrDigiRoot = (TkrDigi*)tkrDigiRef.GetObject();
        } else {
            log << MSG::WARNING << "Could not located TkrDigi TDS/ROOT pair" << endreq;
            continue;
        }

        // Now the thing it depends on
        if (m_common.m_mcPosHitMap.find(posHitTds) != m_common.m_mcPosHitMap.end()) {
            posHitRef  = m_common.m_mcPosHitMap[posHitTds];
            posHitRoot = (McPositionHit*)posHitRef.GetObject();
        // Note that noise hits will not have any McPositionHit associated with them
        } else {
        //    log << MSG::WARNING << "Could not located McPositionHit TDS/ROOT pair" << endreq;
            continue;
        }

        if (tkrDigiRoot && posHitRoot) addRootRelation(tkrDigiRef, posHitRef, (*relation)->getInfos());
    }

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
    // When processing data, no MC is available
    if (!calTable) return sc;

    tabType::const_iterator relation;

    for (relation = calTable->begin(); relation != calTable->end(); relation++) {
        Event::CalDigi*          calDigiTds = (*relation)->getFirst();
        Event::McIntegratingHit* intHitTds  = (*relation)->getSecond();

        CalDigi *calDigiRoot = 0;
        TRef     calDigiRef  = 0;
        TRef     intHitRef   = 0;
        McIntegratingHit *intHitRoot = 0;
        if (m_common.m_mcIntHitMap.find(intHitTds) != m_common.m_mcIntHitMap.end()) {
            intHitRef  = m_common.m_mcIntHitMap[intHitTds];
            intHitRoot = (McIntegratingHit*)intHitRef.GetObject();
        // It can happen that a noise hit will not have an McIntegratingHit associated
        } else {
            //log << MSG::WARNING << "Could not located McIntegratingHit TDS/ROOT pair" << endreq;
            continue;
        }

        if (m_common.m_calDigiMap.find(calDigiTds) != m_common.m_calDigiMap.end()) {
            calDigiRef  = m_common.m_calDigiMap[calDigiTds];
            calDigiRoot = (CalDigi*)calDigiRef.GetObject();
        } else {
            log << MSG::WARNING << "Could not located CalDigi TDS/ROOT pair" << endreq;
        }

        if (calDigiRoot && intHitRoot) addRootRelation(calDigiRef, intHitRef, (*relation)->getInfos());
    }

    return sc;
}

StatusCode relationRootWriterAlg::writeCalReconRelations() {
    // Purpose and Method:  Retrieve the relations concerning CalDigi from TDS
    //    Write the relations to ROOT

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    Event::CalClusterHitTabList* xTal2ClusTabList = 
        SmartDataPtr<Event::CalClusterHitTabList>(eventSvc(), EventModel::CalRecon::CalClusterHitTab);

    if (!xTal2ClusTabList) return sc;

    for (Event::CalClusterHitTabList::iterator calRelIter = xTal2ClusTabList->begin(); calRelIter != xTal2ClusTabList->end(); calRelIter++)
    {
        Event::CalXtalRecData* xtalRecDataTds = (*calRelIter)->getFirst();
        Event::CalCluster*     calClusterTds  = (*calRelIter)->getSecond();

        CalXtalRecData* calXtalRecDataRoot = 0;
        TRef            calXtalRecDataRef  = 0;
        CalCluster*     calClusterRoot     = 0;
        TRef            calClusterRef      = 0;
        if (m_common.m_calXtalRecDataMap.find(xtalRecDataTds) != m_common.m_calXtalRecDataMap.end()) 
        {
            calXtalRecDataRef  = m_common.m_calXtalRecDataMap[xtalRecDataTds];
            calXtalRecDataRoot = (CalXtalRecData*)calXtalRecDataRef.GetObject();
        // It can happen that a noise hit will not have an McIntegratingHit associated
        } else {
            //log << MSG::WARNING << "Could not located McIntegratingHit TDS/ROOT pair" << endreq;
            continue;
        }

        if (m_common.m_calClusterMap.find(calClusterTds) != m_common.m_calClusterMap.end()) {
            calClusterRef  = m_common.m_calClusterMap[calClusterTds];
            calClusterRoot = (CalCluster*)calClusterRef.GetObject();
        } else {
            log << MSG::WARNING << "Could not located CalDigi TDS/ROOT pair" << endreq;
        }

        if (calClusterRoot && calXtalRecDataRoot) addRootRelation(calXtalRecDataRef, calClusterRef, (*calRelIter)->getInfos());
    }

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

    // Retrieve the list of relations between McParticles and McTrajectory's
    // which should be one to one
    SmartDataPtr<Event::McPartToTrajectoryTabList> 
        mcPartToTraj(eventSvc(), "/Event/MC/McPartToTrajectory");

    // Look up McParticle <--> McTrajectory relational table
    if (!mcPartToTraj) return sc;

    // Recover the McPositionHits to trajectory points relational tables
    SmartDataPtr<Event::McPointToPosHitTabList> 
        mcPointToPosHitList(eventSvc(), "/Event/MC/McPointToPosHit");
    Event::McPointToPosHitTab mcPosHitTab(mcPointToPosHitList);

    // Recover the McPositionHits to trajectory points relational tables
    SmartDataPtr<Event::McPointToIntHitTabList> 
        mcPointToIntHitList(eventSvc(), "/Event/MC/McPointToIntHit");
    Event::McPointToIntHitTab mcIntHitTab(mcPointToIntHitList);

    // Loop through this list at the top level
    for(Event::McPartToTrajectoryTabList::const_iterator trajIter = mcPartToTraj->begin(); 
        trajIter != mcPartToTraj->end(); trajIter++) 
    {
        Event::McPartToTrajectoryRel* partToTraj = *trajIter;

        const Event::McParticle*   mcPartTds = partToTraj->getFirst();
        const Event::McTrajectory* mcTrajTds = partToTraj->getSecond();

        McParticle*   mcPartRoot = 0;
        TRef          mcPartRef  = 0;
        McTrajectory* mcTrajRoot = 0;
        TRef          mcTrajRef  = 0;

        if (m_common.m_mcPartMap.find(mcPartTds) != m_common.m_mcPartMap.end()) {
            mcPartRef  = m_common.m_mcPartMap[mcPartTds];
            mcPartRoot = (McParticle*)mcPartRef.GetObject();
        } else {
            log << MSG::WARNING << "Could not located McParticle TDS/ROOT pair" << endreq;
        }

        if (m_common.m_mcTrajectoryMap.find(mcTrajTds) != m_common.m_mcTrajectoryMap.end()) {
            mcTrajRef  = m_common.m_mcTrajectoryMap[mcTrajTds];
            mcTrajRoot = (McTrajectory*)mcTrajRef.GetObject();
        } else {
            log << MSG::WARNING << "Could not located McTrajectory TDS/ROOT pair" << endreq;
        }

        // Pointers must be found at this level to proceed
        if (mcPartRoot && mcTrajRoot)
        {
            addRootRelation(mcPartRef, mcTrajRef, (*trajIter)->getInfos());

            // Set up to loop through Trajectory points
            std::vector<Event::McTrajectoryPoint*> points = mcTrajTds->getPoints();
            std::vector<Event::McTrajectoryPoint*>::const_iterator pointIter;

            for(pointIter = points.begin(); pointIter != points.end(); pointIter++) 
            {
                // De-reference the McTrajectoryPoint pointer
                Event::McTrajectoryPoint* mcPointTds  = *pointIter;
                McTrajectoryPoint*        mcPointRoot = 0;
                TRef                      mcPointRef  = 0;

                if (m_common.m_mcTrajectoryPointMap.find(mcPointTds) != m_common.m_mcTrajectoryPointMap.end()) {
                    mcPointRef  = m_common.m_mcTrajectoryPointMap[mcPointTds];
                    mcPointRoot = (McTrajectoryPoint*)mcPointRef.GetObject();
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
                    TRef                  mcPosHitRef  = 0;

                    if (m_common.m_mcPosHitMap.find(mcPosHitTds) != m_common.m_mcPosHitMap.end()) {
                        mcPosHitRef  = m_common.m_mcPosHitMap[mcPosHitTds];
                        mcPosHitRoot = (McPositionHit*)mcPosHitRef.GetObject();
                    } else {
                        log << MSG::WARNING << "Could not located McPositionHit TDS/ROOT pair" << endreq;
                    }

                    if (mcPosHitRoot) addRootRelation(mcPointRef, mcPosHitRef, hitRel->getInfos());
                }

                // Now look for McTrajectoryPoint <--> McIntegratingHit
                Event::McPointToIntHitVec intRelVec = mcIntHitTab.getRelByFirst(mcPointTds);

                if (!intRelVec.empty())
                {
                    Event::McPointToIntHitRel* hitRel = *(intRelVec.begin());

                    Event::McIntegratingHit* mcIntHitTds  = hitRel->getSecond();
                    McIntegratingHit*        mcIntHitRoot = 0;
                    TRef                     mcIntHitRef  = 0;

                    if (m_common.m_mcIntHitMap.find(mcIntHitTds) != m_common.m_mcIntHitMap.end()) {
                        mcIntHitRef  = m_common.m_mcIntHitMap[mcIntHitTds];
                        mcIntHitRoot = (McIntegratingHit*)mcIntHitRef.GetObject();
                    } else {
                        log << MSG::WARNING << "Could not located McIntegratingHit TDS/ROOT pair" << endreq;
                    }

                    if (mcIntHitRoot) addRootRelation(mcPointRef, mcIntHitRef, hitRel->getInfos());
                }
            }
        }
    }

    return sc;
}

void relationRootWriterAlg::addRootRelation(TRef& first, TRef& second, std::vector<std::string> infos)
{
    // Our first job is to copy any string information in the "infos" vector
    TObjArray infosRoot;

    // Skip the real work if the vector is empty
    if (!infos.empty())
    {
        for(std::vector<std::string>::iterator infosItr = infos.begin(); infosItr != infos.end(); infosItr++)
        {
            TObjString* tString = new TObjString((*infosItr).data());
            infosRoot.Add(tString);
        }
    }

    // make sure we "own" the infos objects
    infosRoot.SetOwner();

    // Create the root version of the relation
    Relation*  rel    = new Relation(first, second, infosRoot);

    // Add to our table (note that duplicate relations can't exist here)
    m_relTable->addRelation(rel);
    
    return;
}

void relationRootWriterAlg::writeEvent() 
{
    // Purpose and Method:  Stores the Relations data for this event in the ROOT
    //    tree.

    m_rootIoSvc->fillTree("rel");

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

    m_rootIoSvc->closeFile("rel");

    return;
}

StatusCode relationRootWriterAlg::finalize()
{
    close();
    
    StatusCode sc = StatusCode::SUCCESS;
    setFinalized();
    return sc;
}


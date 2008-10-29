#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"

//#include "Trigger/TriRowBits.h"

//#include "OnboardFilterTds/FilterStatus.h"

#include "facilities/Util.h"

#include "RootIo/IRootIoSvc.h"
//#include "RootConvert/Utilities/RootReaderUtil.h"
#include "ntupleWriterSvc/INTupleWriterSvc.h"

#include <map>
#include <string>

#include "TChain.h"
#include "TLeaf.h"
#include "TROOT.h"

// ADDED FOR THE FILE HEADERS DEMO
//#include "RootIo/FhTool.h"

/** @class ntupleRootReaderAlg
 * @brief Reads ntuple data from a persistent ROOT file and stores the
 * the data in the TDS in /Event 
 *
 * @author Heather Kelly
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/ntupleRootReaderAlg.cxx,v 1.97 2008/10/05 05:23:04 heather Exp $
 */

class ntupleRootReaderAlg : public Algorithm
{	
public:
    
    ntupleRootReaderAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    /// Handles setup by opening ROOT file in read mode and creating a new TTree
    StatusCode initialize();
   
    /// Orchastrates reading from ROOT file and storing the data on the TDS for each event
    StatusCode execute();
    
    /// Closes the ROOT file and cleans up
    StatusCode finalize();

    std::string ntupleRootReaderAlg::getItem( const std::string& itemName, void*& pval);

            
private:

    /// Closes the ROOT file
    void close();

    /// Array of input file names
    StringArrayProperty m_fileListProp;
    std::vector<std::string> m_fileList;
    /// name of the merit TTree stored in the ROOT file
    std::string m_treeName;

    TChain *m_chain;

    Long64_t m_numEvents;
  
    IRootIoSvc* m_rootIoSvc;
    INTupleWriterSvc *m_rootTupleSvc;

   /// collection of leaf addresses for items that we have to create an object for.
    /// This occurs in reprocessing, when not all AnaTup Tools are executed - so not all branches have a corresponding
    /// variable to TChain::SetBranchAddress for, so that we have a stable location to provide via the getItem call.
    std::map<std::string, void*> m_itemPool;


};

static const AlgFactory<ntupleRootReaderAlg>  Factory;
const IAlgFactory& ntupleRootReaderAlgFactory = Factory;


ntupleRootReaderAlg::ntupleRootReaderAlg(const std::string& name, ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator)
{
    // Input pararmeters that may be set via the jobOptions file
    // Input ROOT file name  provided for backward compatibility, digiRootFileList is preferred
    // The files will be ignored if RootIoSvc is provided a meta (event collection) ROOT file
    StringArrayProperty initList;
    std::vector<std::string> initVec;
    initVec.clear();
    initList.setValue(initVec);
    //declareProperty("tupleRootFileList",m_fileListProp=initList);
    // Input TTree name
    declareProperty("tupleTreeName", m_treeName="MeritTuple");
}

StatusCode ntupleRootReaderAlg::initialize()
{
    // Purpose and Method:  Called once before the run begins.  This method
    //    opens a new ROOT file and prepares for reading.

    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    // Use the Job options service to set the Algorithm's parameters
    // This will retrieve parameters set in the job options file
    setProperties();

    m_itemPool.clear();
   
    m_rootIoSvc = 0 ;
    if ( service("RootIoSvc", m_rootIoSvc, true).isFailure() ){
        log << MSG::ERROR << "Couldn't find the RootIoSvc!" << endreq;
        m_rootIoSvc = 0;
        // We truly depend on the RootIoSvc now to read/write files, so we cannot just
        // continue if this service is not found.
        return StatusCode::FAILURE;
    } 

    // Retrieve tuple info an set up
    // get a pointer to RootTupleSvc
    sc  = service("RootTupleSvc", m_rootTupleSvc, true);
    if( sc.isFailure() )
    {
        MsgStream log( msgSvc(), name() );
        log << MSG::WARNING << "failed to get the RootTupleSvc" << endreq;
        m_rootTupleSvc = 0;
    }


    //if ( m_fileList.value().size() <= 0)  {
       if (m_rootTupleSvc) {
         log << MSG::INFO << "looking to RootTupleSvc for input file list" 
             << endreq;
         bool stat = m_rootTupleSvc->getInputFileList(m_fileList);
         if (m_fileList.size() <= 0) {
             log << MSG::ERROR << "No input files specified" << endreq;
             return StatusCode::FAILURE; 
         }
       } else {
         log << MSG::ERROR << "No input tuple specified" << endreq;
         return StatusCode::FAILURE;
       }
    //} 



   // Assuming all input ROOT files have the same tree name,
   // we could adjust the code to pick out the TTree and provide 
   // its name in the TChain::Add call
   m_chain = new TChain(m_treeName.c_str());

   // Add all files in the JO list
   std::vector<std::string>::const_iterator fileListItr;
   for (fileListItr = m_fileList.begin();
        fileListItr != m_fileList.end();
        fileListItr++ ) {
        std::string fileName = *fileListItr;
        facilities::Util::expandEnvVar(&fileName);
        int stat = m_chain->Add(fileName.c_str());
        if (stat <= 0) {
            log << MSG::WARNING << "Failed to TChain::Add " 
                << fileName << " return code: " << stat << endreq;
            return StatusCode::FAILURE;
        }
        else
            log << MSG::DEBUG << "Added File: " << fileName << endreq;

    }  // end for loop TChain initialized
    // call GetEntries to load the headers of the TFiles
    m_numEvents = m_chain->GetEntries();
    log << MSG::INFO << "Number of events in input files = " 
        << m_numEvents << " StartingIndex: " << m_rootIoSvc->index() << endreq;
    if ((m_rootIoSvc->index() > m_numEvents-1) || (m_rootIoSvc->index() < 0)) {
        log << MSG::WARNING << "StartingIndex invalid" << endreq;
    }
    int numbytes = m_chain->GetEntry(m_rootIoSvc->index());
    if (numbytes <= 0) 
       log << MSG::WARNING << "Unable to read tuple event, "
           << m_rootIoSvc->index() << endreq;
    // Here is our chance to set up branch pointers for the whole input TChain,
    //  so that no elements are missed
    TObjArray *brCol = m_chain->GetListOfBranches();
    int numBranches = brCol->GetEntries();
    int iBranch;
    for (iBranch=0;iBranch<numBranches;iBranch++) {
        std::string branchName(((TBranch*)(brCol->At(iBranch)))->GetName());
        std::string leafName = branchName;
        int index = leafName.find("[");
        if(index!=std::string::npos) leafName = leafName.substr(0,index); 
        log << MSG::DEBUG << "setting branch: " << branchName
            << " and Leaf: " << leafName << endreq;
        TLeaf *leaf = ((TBranch*)brCol->At(iBranch))->GetLeaf(leafName.c_str());
        if (!leaf) {
           log << MSG::WARNING << "Leaf: " << leafName << " not found" << endreq;
           continue;
         }
         std::string type_name = leaf->GetTypeName();
         if (type_name == "Float_t")
         {
            m_itemPool[branchName] = new Float_t[leaf->GetNdata()];
         }
         else if (type_name == "Int_t")
         {
            m_itemPool[branchName]= new Int_t[leaf->GetNdata()];
         }
         else if (type_name == "UInt_t")
         {
            m_itemPool[branchName]= new UInt_t[leaf->GetNdata()];
         }        
         else if (type_name == "ULong64_t")
         {
            m_itemPool[branchName]= new ULong64_t[leaf->GetNdata()];
         }
         else if (type_name == "Double_t")
         {
            m_itemPool[branchName] = new Double_t[leaf->GetNdata()];
         }
         else if (type_name == "Char_t")
         {
            m_itemPool[branchName] = new Char_t[leaf->GetNdata()];
         }
         else
         {
           log << MSG::WARNING << "type: " << type_name <<" not found" << endreq;
         }
         m_chain->SetBranchAddress(branchName.c_str(), m_itemPool[branchName]);
     }


    return sc;
    
}

std::string ntupleRootReaderAlg::getItem( const std::string& itemName, void*& pval) {
 
    TDirectory *saveDir = gDirectory;

    TLeaf *leaf = m_chain->GetLeaf(itemName.c_str());
    std::string type_name = "";
    if (!leaf) {
        pval = 0;
        return(type_name);
    }

    pval = leaf->GetValuePointer();

    type_name = leaf->GetTypeName();

    saveDir->cd();
    return type_name;


}

StatusCode ntupleRootReaderAlg::execute()
{
    // Purpose and Method:  Called once per event.  This method calls
    //   the appropriate methods to read data from the ROOT file and store
    //   data on the TDS.

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Assuming we are always reading via index for now
    Int_t numBytes = m_chain->GetEvent(m_rootIoSvc->index());

    if (numBytes <= 0) {
         log << MSG::ERROR << "Failed to load event " <<  m_rootIoSvc->index()
             << endreq;
         return StatusCode::FAILURE;
    }

    // Retrieve the Event data for this event
    SmartDataPtr<Event::EventHeader> evt(eventSvc(), EventModel::EventHeader);
    if (!evt) {
        log << MSG::ERROR << "Failed to retrieve /Event on TDS" << endreq;
        return StatusCode::FAILURE;
    }

    void *runRoot;
    void *evtIdRoot;
    std::string typeRun = getItem("EvtRun", runRoot);
    std::string typeEventId = getItem("EvtEventId", evtIdRoot);

    evt->setEvent(*((UInt_t*)evtIdRoot));
    evt->setRun(*((UInt_t*)runRoot));

    log << MSG::DEBUG << std::dec << "Reading Event (run, event): (" << *((UInt_t*)runRoot)
        << ", " << *((UInt_t*)evtIdRoot) << ")" << endreq;

    void* liveTimeRoot;
    std::string typeLiveTime = getItem("EvtLiveTime", liveTimeRoot);
    evt->setLivetime(*((Double_t*)liveTimeRoot));

    void* elapsedTimeRoot;
    std::string typeElaspsedTime = getItem("EvtElapsedTime", elapsedTimeRoot);
    evt->setTime(*((Double_t*)elapsedTimeRoot));

    //evt->setTrigger(m_digiEvt->getL1T().getTriggerWord());
    //evt->setTriggerWordTwo(m_digiEvt->getL1T().getTriggerWordTwo());
    //evt->setGemPrescale(m_digiEvt->getL1T().getGemPrescale());
    //evt->setGltPrescale(m_digiEvt->getL1T().getGltPrescale());
    //evt->setPrescaleExpired(m_digiEvt->getL1T().getPrescaleExpired());

    return sc;
}

void ntupleRootReaderAlg::close() 
{

    if (m_chain) {
        delete m_chain;
        m_chain = 0;
    }

}

StatusCode ntupleRootReaderAlg::finalize()
{
    close();
    
    StatusCode sc = StatusCode::SUCCESS;
    setFinalized();
    return sc;
}


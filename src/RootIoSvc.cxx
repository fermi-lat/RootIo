
/** 
* @file RootIoSvc.cxx
* @brief definition of the class RootIoSvc
*
*  $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/RootIo/src/RootIoSvc.cxx,v 1.70 2013/05/17 13:47:43 heather Exp $
*  Original author: Heather Kelly heather@lheapop.gsfc.nasa.gov
*/

#include "commonData.h"
#include "RootInputDesc.h"
#include "RootOutputDesc.h"
#include "GleamMessageHandler.h"
#include "RootConvert/Utilities/RootReaderUtil.h"
#include "rootUtil/CelManager.h"
#include "rootUtil/CompositeEventList.h"

#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/GaudiException.h"
//#include "GaudiKernel/IObjManager.h"
//#include "GaudiKernel/IToolFactory.h"
#include "GaudiKernel/IAlgManager.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/IAppMgrUI.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/IRunable.h"
#include "GaudiKernel/Property.h"
#include "RootIo/IRootIoSvc.h"
#include "ntupleWriterSvc/INTupleWriterSvc.h"
#include "RootIo/FhTool.h"
#include "GlastSvc/GlastRandomSvc/IRandomAccess.h"

#include "TSystem.h"
#include "TFile.h"
#include "facilities/Util.h"
#include "facilities/Timestamp.h"

#include <vector>
#include <algorithm>

//forward declarations
template <class TYPE> class SvcFactory;
class IAppMgrUI;



//==============================================================================
/**
* \class RootIoSvc
*
* \brief Service that implements the IRunable interface, to control the event loop.
* \author Heather Kelly heather@slac.stanford.edu
*/
//==============================================================================


class RootIoSvc : 
  virtual public Service, 
  virtual public IIncidentListener,
  virtual public IRootIoSvc,
  virtual public IRunable
 {  
  public:  

    //====================
    // Gaudi machinery
    //====================

    /// for the IRunnable interfce
    virtual StatusCode run() ;

    //------------------------------
    //  stuff required by a Service
    //------------------------------
    
    /// perform initializations for this service. 
    virtual StatusCode initialize() ;
    
    /// perform the finalization, as required for a service.
    virtual StatusCode finalize() ;
    
    /// Query interface
    virtual StatusCode queryInterface( const InterfaceID & riid, void** ppvUnknown ) ;

    /// Handles incidents, implementing IIncidentListener interface
    virtual void handle( const Incident & ) ;    
	
    
    //====================
    // For HepRepSvc
    //====================
    
    virtual bool setRootFile
     ( const char * mc, const char * digi, 
       const char * rec, const char* relation, const char * gcr ) ;

    virtual bool setIndex( Long64_t i ) ;
    virtual Long64_t index() { return m_index ; }

    virtual bool setRunEventPair(std::pair<int,int> ids) ;
    virtual std::pair<int,int> runEventPair() { return m_runEventPair ; }

    virtual Long64_t getIndexByEventID(int run, int event);

    
    //====================
    // For readers
    //====================
    
    // file list
    virtual bool setFileList( const std::string & type, const StringArrayProperty & fileList ) ;
    virtual StringArrayProperty getFileList( const std::string & type) const ;
    virtual bool appendFileList( StringArrayProperty & fileList, const std::string & fileName ) ;

    virtual StatusCode prepareRootInput
     ( const std::string&         type, 
       const std::string&         tree,
       const std::string&         branch,
       TObject**                  branchPtr,
       const StringArrayProperty& fileList ) ;

    virtual bool setBranchStatus
     ( const std::string&         type,
       const std::string&         branch,
       int                        status );

 
    virtual StatusCode closeInput(const std::string& type);
       
    virtual TObject * getNextEvent( const std::string & type, long long inputIndex) ;

    /// Max Event for event loop for read jobs
    virtual long long getEvtMax() const {return m_evtMax;}

    virtual long long getRootEvtMax(const std::string& type) const;

    /// Allows reader algorithms to register the number of events in their
    /// input ROOT files
    virtual void setEvtMax(Long64_t max) ;
    virtual void setRootEvtMax(Long64_t max) {setEvtMax(max);}

    virtual bool terminateOnReadError() { return m_terminateOnReadError; }


    //====================
    // For writers
    //====================

    virtual TTree * prepareRootOutput
      ( const std::string & type,
        const std::string & fileName,
        const std::string & treeName,
        int compressionLevel,
        const std::string & treeTitle ) ;

    //[David] APPARENTLY NEVER USED ?!?
    virtual TTree * getTree( const std::string & type) ;

    virtual StatusCode setupBranch(const std::string & type,
        const std::string & branchName, 
        const std::string & classname, void * branchAddr,
        int bufSize =64000, int splitLevel =1 ) ;

    virtual StatusCode closeFile( const std::string & type) ;

    virtual StatusCode fillTree( const std::string & type) ;
    
    virtual int getAutoSaveInterval() { return m_autoSaveInterval ; }

    virtual CompositeEventList* getCel() { return m_outputCel; }

  protected : 
    
    RootIoSvc ( const std::string & name, ISvcLocator * al ) ;
    virtual ~RootIoSvc() ;
    
    
  private :

    const RootInputDesc * getRootInputDesc( const std::string & type ) const ;
    RootInputDesc * getRootInputDesc( const std::string & type ) ; 
    StatusCode deleteRootInputDesc( const std::string& type );
    unsigned setSingleFile( const std::string & type, const char * fileName ) ;
    void registerRootTree( TChain * ch ) ;

    const RootOutputDesc * getRootOutputDesc( const std::string & type ) const;
    RootOutputDesc * getRootOutputDesc( const std::string & type ) ;
    void setOutputUpdate( bool flag ) ;
    bool checkOutputUpdate( const std::string & type ) ;
    bool checkOutputUpdate() ;

    void beginEvent() ;
    void endEvent() ;
    void loopStatus( long long eventNumber, double currentTime, MsgStream & log) ;	
    /// Allow SvcFactory to instantiate the service.
    friend class SvcFactory<RootIoSvc> ;
   
    /// Reference to application manager UI
    IAppMgrUI * m_appMgrUI ;
    IRandomAccess *m_randTool; // setup RandomNum Tool
    //IntegerProperty m_evtMax ;
    long long m_evtMax ;
    IntegerProperty m_autoSaveInterval ;
    UnsignedLongProperty m_treeSize ;
    StringProperty m_tupleName;
    IntegerProperty m_noFailure;
    bool m_rebuildIndex;
    bool m_abortOnRootError;
    bool m_autoFlush;
    int m_readRetries;
    bool m_terminateOnReadError;
    std::string m_compressAlg;
    int m_compressionAlgFlag;

    // starting and ending times for orbital simulation
    DoubleProperty m_startTime ;
    DoubleProperty m_endTime ;

    unsigned int m_rootEvtMax ;
    Long64_t m_index ;
    UnsignedLongProperty m_startIndex ;
    std::pair< int, int > m_runEventPair ;
    
    commonData m_common ;
    UInt_t m_objectNumber ;

    bool m_useIndex ;

    std::map< std::string, RootInputDesc * > m_rootIoMap ;
    std::map< std::string, RootOutputDesc *> m_rootOutputMap ;

   
    CelManager m_celManager ; // so to read
    CompositeEventList * m_outputCel ; // so to write
    std::vector<TTree*> m_celTreeCol ; // so to write
    std::string m_celFileNameWrite, m_celFileNameRead ;

    /// access the ntupleWriter service to access the output Ntuple
    INTupleWriterSvc *   m_rootTupleSvc;
    TTree *m_outputTuple;

    GleamMessageHandler* m_rootFileMessageHandler;

    IToolSvc * m_gaudiToolSvc;
    //IFhTool * m_headersTool;
    
 } ;



//==============================================================================
//
// RootIoSvc definition
//
//==============================================================================


// declare the service factories for the RootIoSvc
//static SvcFactory<RootIoSvc> a_factory ;
//const ISvcFactory& RootIoSvcFactory = a_factory ;
DECLARE_SERVICE_FACTORY(RootIoSvc);

/// Standard Constructor
RootIoSvc::RootIoSvc(const std::string& name,ISvcLocator* svc)
: Service(name,svc), m_outputCel(0), m_rootTupleSvc(0), m_outputTuple(0)
{
    
    //declareProperty("EvtMax"     , m_evtMax=0);

    // disable for now
    //declareProperty("StartTime"   , m_startTime=0);
    //declareProperty("EndTime",      m_endTime=0);
    m_endTime = 0;
    m_startTime = 0;

    declareProperty("AutoSaveInterval", m_autoSaveInterval=1000);
    declareProperty("StartingIndex", m_startIndex=0);
    // limited by the size of an unsigned int, expecting MB
    declareProperty("MaxTreeSize", m_treeSize=0);
    declareProperty("TupleName", m_tupleName="MeritTuple");
    declareProperty("NoFailure", m_noFailure=0);
    declareProperty("RebuildIndex", m_rebuildIndex=false);
    
    declareProperty("CompressAlg", m_compressAlg="LZMA");   // Was ZLIB (TU 7/1/13)
    // 
    declareProperty("CelRootFileWrite", m_celFileNameWrite="");
    declareProperty("CelRootFileRead", m_celFileNameRead="");

    declareProperty("AbortOnRootError", m_abortOnRootError=true);
    declareProperty("AllowAutoFlush", m_autoFlush=false);

    declareProperty("ReadRetries", m_readRetries=3);

    declareProperty("TerminateOnReadError", m_terminateOnReadError = false);


    m_index = 0;
    m_rootEvtMax = 0;
    m_runEventPair = std::pair<int,int>(-1,-1);
    m_useIndex = false;

    //m_headersTool = 0;
    m_gaudiToolSvc = 0;

#ifdef WIN32
    gSystem->Load("libTreePlayer.dll");
#endif

}


/// Standard Destructor
RootIoSvc::~RootIoSvc()  
 {
  m_rootIoMap.clear() ;
  m_rootOutputMap.clear();
  m_celTreeCol.clear();
  delete m_outputCel ;
 }

StatusCode RootIoSvc::initialize () 
{   
    StatusCode  status =  Service::initialize();
    
    // bind all of the properties for this service
    setProperties ();
    
    m_appMgrUI = 0 ;
    status = serviceLocator()->queryInterface(IAppMgrUI::interfaceID(), (void**)&m_appMgrUI);
    
    // use the incident service to register begin, end events
    IIncidentSvc* incsvc = 0;
    status = service ("IncidentSvc", incsvc, true);

    if( status.isFailure() ) return status;

    incsvc->addListener(this, "BeginEvent", 100);
    incsvc->addListener(this, "EndEvent", 0);

    // Tell ROOT to reset signals to their default behavior
    gSystem->ResetSignal(kSigBus); 
    gSystem->ResetSignal(kSigSegmentationViolation); 
    gSystem->ResetSignal(kSigIllegalInstruction); 
    gSystem->ResetSignal(kSigFloatingException);  

    // Create a local MessageHandler to check for ROOT errors
    m_rootFileMessageHandler = new GleamMessageHandler("TFile", kTRUE);

    if (m_treeSize > 0) {
        // JO parameter in MB - need to convert to bytes
        Long64_t maxTreeSize = m_treeSize * 1000000;
        TTree::SetMaxTreeSize(maxTreeSize);
    } else if (m_treeSize == 0) {
        // 500 GB default
        //Long64_t maxTreeSize = 25000000000;
        Long64_t maxTreeSize = 500000000000LL;
        TTree::SetMaxTreeSize(maxTreeSize);
    }

    MsgStream log2( msgSvc(), name() );
    m_compressionAlgFlag = (m_compressAlg.compare("LZMA")) ? ROOT::kZLIB : ROOT::kLZMA;
    log2 << MSG::INFO << "Setting compression alg " << m_compressAlg << " "
        << m_compressionAlgFlag << endreq;

    // DC : I guess the original author
    // was meaning to test m_index, or the
    // test has no sense
    m_index = m_startIndex ;
    if (m_index >= 0) m_useIndex = true;


    // Retrieve tuple info an set up
    // get a pointer to RootTupleSvc
    status = service("RootTupleSvc", m_rootTupleSvc, true);       
    if( status.isFailure() ) 
    {
        MsgStream log( msgSvc(), name() );
        log << MSG::WARNING << "failed to get the RootTupleSvc" << endreq;
        m_rootTupleSvc = 0;
    }

    if (!m_celFileNameWrite.empty()) 
      {
       //m_celManager.initWrite(m_celFileNameWrite, "RECREATE");
       m_outputCel = new CompositeEventList(m_celFileNameWrite,"RECREATE") ;
      }
    
    if (!m_celFileNameRead.empty()) {
        MsgStream log( msgSvc(), name() );
        log << MSG::INFO << "Cel ROOT file Reading" << endreq;
        m_celManager.initRead(m_celFileNameRead);
    }


    StatusCode toolSvcSc = service("ToolSvc", m_gaudiToolSvc, true);
    if (toolSvcSc.isSuccess())
        if (m_gaudiToolSvc->retrieveTool("RootIoRandom", m_randTool).isFailure() ) {
            MsgStream log(msgSvc(), name() );
            log << MSG::WARNING << "Failed to create RootIoRandom" << endreq;
   }
/*
    if  ( toolSvcSc.isSuccess() ) {
        StatusCode headersSc = m_gaudiToolSvc->retrieveTool("FhTool", m_headersTool);
        if (headersSc.isFailure() ) {
            MsgStream log( msgSvc(), name() );
            m_headersTool = 0;
            log << MSG::WARNING << "Failed to set up FhTool" << endreq;
        }
     } else {
        MsgStream log( msgSvc(), name() );
        log << MSG::WARNING << "Failed to set up FhTool" << endreq;
     }

*/
    return StatusCode::SUCCESS;
}



//============================================================================
//
// RootInputDesc
//
//============================================================================


const RootInputDesc * RootIoSvc::getRootInputDesc( const std::string & type ) const
 {
  const RootInputDesc * rootInputDesc = 0 ;
  std::map< std::string, RootInputDesc * >::const_iterator typeItr = m_rootIoMap.find(type) ;
  if (typeItr!=m_rootIoMap.end()) rootInputDesc = typeItr->second ;
  return rootInputDesc ;
 }

RootInputDesc * RootIoSvc::getRootInputDesc( const std::string & type )
 {
  return const_cast<RootInputDesc *>
   ( static_cast<const RootIoSvc *>(this)->getRootInputDesc(type) ) ;
 }

StatusCode RootIoSvc::deleteRootInputDesc(const std::string& type)
{
    StatusCode sc = StatusCode::SUCCESS;

    std::map<std::string,RootInputDesc*>::iterator descItr = m_rootIoMap.find(type);

    if (descItr != m_rootIoMap.end())
    {
        delete descItr->second;
        m_rootIoMap.erase(descItr);
    }
    else
    {
        MsgStream log( msgSvc(), name() );

        log << MSG::ERROR << "Failed to find " << type << " in the input descriptor map" << endreq;

        sc = StatusCode::FAILURE;
    }
    return sc;
}


//============================================================================
//
// FileList
//
//============================================================================


bool RootIoSvc::setFileList( const std::string & type, const StringArrayProperty & fileNames )
 {
  bool fileSet = false ;
  RootInputDesc * rootInputDesc = getRootInputDesc(type) ;
  if (rootInputDesc)
   {
    Long64_t retVal = rootInputDesc->setFileList(fileNames) ;
    if (retVal > 0) fileSet = true ;
   }
  else
   {
    fileSet = false ;
   }
  return fileSet ;
 }

StringArrayProperty RootIoSvc::getFileList( const std::string & type) const
 {
  StringArrayProperty result ;
  const RootInputDesc * rootInputDesc = getRootInputDesc(type) ;
  if (rootInputDesc) result = rootInputDesc->getFileList() ;
  return result ;
 }

bool RootIoSvc::appendFileList(StringArrayProperty &fileList, const std::string &fileName) {
    bool status = false;
    if (fileName.empty())  return status;

    std::vector<std::string> initVec;
    typedef std::vector<std::string>::const_iterator StringVecIter ;
    StringVecIter fileListItr ;
    for ( fileListItr = fileList.value().begin() ; 
          fileListItr != fileList.value().end() ;
          fileListItr++ )
    {
          std::string file = *fileListItr ;
          initVec.push_back(file);
    }
    initVec.push_back(fileName);
    fileList.setValue(initVec);
    status = true;
    return status;
    
}

unsigned RootIoSvc::setSingleFile( const std::string & type, const char * fileName )
 {
  // return value :
  // 0 : all ok
  // 1 : empty filename
  // 2 : unknown type
  // 3 : unknown file
  
  MsgStream log (msgSvc(), name() );
  log << MSG::DEBUG << "setSingleFile type: " << type << " file: " 
      << fileName << endreq;
  std::string myFileName = fileName ;
  facilities::Util::expandEnvVar(&myFileName) ;
  if (myFileName.empty())
   { 
    closeInput(type); // remove the corresponding RootInputDesc object
    return 1 ; 
   }
   
  RootInputDesc * rootInputDesc = getRootInputDesc(type) ;
  if (rootInputDesc==0)
   { 
     log << MSG::DEBUG 
         << "setSingleFile rootInputDesc does not exist for type " 
         << type << endreq;

      // For now we only allow MC files to be added on the fly (for WIRED)
      if (type == "mc") {
          StringArrayProperty fileList;
          std::vector<std::string> initVec;
          initVec.push_back(myFileName);
          fileList.setValue(initVec);
          prepareRootInput("mc", "Mc", "McEvent", 0, fileList);
          rootInputDesc = getRootInputDesc(type);
          if (rootInputDesc == 0) return 2;
      } else
         return 2 ; 
   }
  
  if (!RootPersistence::fileExists(myFileName))
   { 
       log << MSG::DEBUG << "setSingleFile file doesn't exist for type " 
           << type << endreq;
       return 3 ; 
   }
  
  std::vector<std::string> fileVec ;
  fileVec.push_back(myFileName) ;
  StringArrayProperty fileList ;
  fileList.setValue(fileVec) ;
  rootInputDesc->setFileList(fileList,true) ;

  // fixAcdStreamer
  if (type=="RECON")
   {
    TFile * f = TFile::Open(myFileName.c_str()) ;
    if (f->IsOpen())
     {
      AcdRecon::fixAcdStreamer(f->GetVersion()) ;
      f->Close() ;
     }
   }
          
  return 0 ;
 }

bool RootIoSvc::setRootFile( const char * mc, const char * digi, const char * rec, const char* relation, const char * gcr ) 
 {
  bool success = false ;

  unsigned mcResult = setSingleFile("mc",mc) ;
  unsigned digiResult = setSingleFile("digi",digi) ;
  unsigned reconResult = setSingleFile("recon",rec) ;
  unsigned relationResult = setSingleFile("rel", relation);
  unsigned gcrResult = setSingleFile("gcr",gcr) ;

  // at least one string must be non-null
  if ( (mcResult==1) &&
       (digiResult==1) &&
       (reconResult==1) &&
       (relationResult==1) &&
       (gcrResult==1) )
   { return success ; }

  // any other problem is an error
  if ( (mcResult>1) ||
       (digiResult>1) ||
       (reconResult>1) ||
       (relationResult>1) ||
       (gcrResult>1) )
   { return success ; }

  success = true ;
  return success ;
 }



//============================================================================
//
// 
//
//============================================================================


StatusCode RootIoSvc::prepareRootInput
 ( const std::string&         type, 
   const std::string&         tree, 
   const std::string&         branch,
   TObject**                  branchPtr,
   const StringArrayProperty& fileList )
{
    MsgStream log( msgSvc(), name() );
    StatusCode sc = StatusCode::SUCCESS;
   
    RootInputDesc *rootInputDesc = 0;
    if (m_celFileNameRead.empty()) {

        // Create a new RootInputDesc object for this algorithm
        rootInputDesc = new RootInputDesc(fileList, tree, branch, m_readRetries,
                                     branchPtr, 
                                     m_rebuildIndex, log.level()<=MSG::DEBUG);
        // Register with RootIoSvc
        if (rootInputDesc->getNumEvents() <= 0) {
            log << MSG::WARNING << "Failed to setup ROOT input for " 
                << type << endreq;
            if(!m_noFailure) return StatusCode::FAILURE;
        }
        setEvtMax(rootInputDesc->getNumEvents());

    } else {  // Event Collection
        MsgStream log( msgSvc(), name() );
        log << MSG::INFO << "Cel ROOT file Reading" << endreq;   

        TChain *chain = m_celManager.getChainByType(type);
        if (!chain) {
            MsgStream log(msgSvc(), name());
            log << MSG::WARNING << "TChain " << type 
                << " could not be constructed from event collection" << endreq;
            if (!m_noFailure) return StatusCode::FAILURE;
        }
        rootInputDesc = new RootInputDesc(chain, tree, branch, m_readRetries,
                          branchPtr, m_rebuildIndex, log.level()<=MSG::DEBUG);
        // Register number of events
        if (m_celManager.getNumEvents() <= 0) {
            log << MSG::WARNING << "Failed to setup ROOT input for " 
                << type << endreq;
            if (!m_noFailure) return StatusCode::FAILURE;
        }
        setEvtMax(m_celManager.getNumEvents());
    }

    // Store the pointer to this in our map
    m_rootIoMap[type] = rootInputDesc;



    return sc;
}

bool RootIoSvc::setBranchStatus(const std::string& type,
  const std::string& branch, int status) {

    MsgStream log( msgSvc(), name() );

    RootInputDesc* rootInputDesc = getRootInputDesc(type);
    if (rootInputDesc) {
        return rootInputDesc->setBranchStatus(branch,status);
    }
    log << MSG::WARNING << "Did not find RootInputDesc type " << type << endreq;
    return false;

}

StatusCode RootIoSvc::closeInput(const std::string& type)
{
    StatusCode sc = StatusCode::SUCCESS;

    RootInputDesc* inputDesc = getRootInputDesc(type);

    if (inputDesc) 
    {
        sc = deleteRootInputDesc(type);
    } 
    else 
    {
        MsgStream log (msgSvc(), name() );
        log << MSG::WARNING << "RootInputDesc of type " << type << " not found" << endreq;
        sc = StatusCode::FAILURE;
    }
    return sc;
}

long long RootIoSvc::getRootEvtMax(const std::string& type) const
{ 
    long long nEventsMax = -1;

    // Standard mode is to return the cached number
    if (type == "")
    {
        nEventsMax = m_rootEvtMax;
    }
    // Otherwise, we want the number of events from a specific file (useful for Overlays)
    else
    {
        if (const RootInputDesc* rootInputDesc = getRootInputDesc(type))
        {
            nEventsMax = rootInputDesc->getNumEvents();
        }
    }

    return nEventsMax; 
}

TObject* RootIoSvc::getNextEvent(const std::string& type, long long inputIndex)
{
    MsgStream log( msgSvc(), name() );

    RootInputDesc* rootInputDesc = getRootInputDesc(type);
    TObject*    pData  = 0;

    if (rootInputDesc)
    {
        // Clear the event occurs in the xxxRootReader algorithm execute methods
        //rootInputDesc->clearEvent();
        
        // Read the event depending upon the mode
        if (!m_celFileNameRead.empty()) {  // access via event collection

            Long64_t ind = m_celManager.getEventIndexInTree(type, index());
            pData = rootInputDesc->getEvent(ind);

        } else if (m_useIndex)  // typical serial reading
        {
            if (inputIndex < 0) pData = rootInputDesc->getEvent(index());
            else                pData = rootInputDesc->getEvent(inputIndex);
        }
        // Otherwise we get by run/event number
        else
        {
            std::pair<int,int> runEvtPair = runEventPair();
            int                runNum     = runEvtPair.first;
            int                evtNum     = runEvtPair.second;
            pData = rootInputDesc->getEvent(runNum,evtNum) ;
            if (pData == 0) log << MSG::DEBUG << "getEvent failed" << endreq;
        }
    }

    return pData;
}

Long64_t RootIoSvc::getIndexByEventID(int run, int event)
{
    std::map<std::string,RootInputDesc *>::iterator typeItr ;
    for ( typeItr = m_rootIoMap.begin() ; typeItr != m_rootIoMap.end() ; ++typeItr )
    {
        RootInputDesc * rootInputDesc = typeItr->second ;
        if(rootInputDesc->checkEventAvailability(run, event)) {
            unsigned index = rootInputDesc->getIndexByEventID(run, event);
            return index;
        }
    }
    return (unsigned)-1;
}

//============================================================================
//
// RootOutputDesc
//
//============================================================================



 const RootOutputDesc * RootIoSvc::getRootOutputDesc( const std::string &type) const
 {
     const RootOutputDesc * rootOutputDesc = 0;
     std::map<std::string, RootOutputDesc *>::const_iterator typeItr = m_rootOutputMap.find(type);
     if (typeItr!=m_rootOutputMap.end()) rootOutputDesc = typeItr->second;
     return rootOutputDesc;
 }

 RootOutputDesc* RootIoSvc::getRootOutputDesc( const std::string &type)
 { 
     return const_cast<RootOutputDesc *>
         (static_cast<const RootIoSvc*>(this)->getRootOutputDesc(type) );
 }

 void RootIoSvc::setOutputUpdate(bool flag) {
     std::map<std::string, RootOutputDesc*>::iterator it;
     for (it = m_rootOutputMap.begin(); it != m_rootOutputMap.end(); it++)
         it->second->setUpdated(flag);
 }

 bool RootIoSvc::checkOutputUpdate(const std::string &type) {
     RootOutputDesc* outputDesc = getRootOutputDesc(type);
     if (!outputDesc) return false;
     return outputDesc->getUpdated();
 }

 bool RootIoSvc::checkOutputUpdate() {
     std::map<std::string, RootOutputDesc*>::const_iterator it;
     for (it = m_rootOutputMap.begin(); it != m_rootOutputMap.end(); it++) {
         if (it->second->getUpdated()) return true;
     }
     return false;
 }
     

TTree* RootIoSvc::prepareRootOutput
        ( const std::string &type,
        const std::string &fileName,
        const std::string &treeName,
        int compressionLevel,
        const std::string &treeTitle) {

    MsgStream log(msgSvc(),name());
    if (getRootOutputDesc(type) != 0) {
        log << MSG::WARNING << "Already found RootOutputDesc entry for type " << type << endreq;
        return 0;
    }
    // Create a new RootOutputDesc object for this algorithm
    RootOutputDesc* outputDesc = new RootOutputDesc(fileName, treeName,
compressionLevel, treeTitle, log.level()<=MSG::DEBUG, m_compressionAlgFlag);

    if (!m_autoFlush) outputDesc->turnOffAutoFlush();

    // Store the pointer to this in our map
    m_rootOutputMap[type] = outputDesc;

    // If a composite event list file is requested, add this component to the CompositeEventList
    if (!m_celFileNameWrite.empty())
     {
      m_outputCel->declareComponent(type) ;
      m_celTreeCol.push_back(outputDesc->getTree());
      //m_celManager.addComponent(treeName, outputDesc->getTree()) ;
     }

    return outputDesc->getTree();

}

//[David] NEVER USED ?!?
 TTree* RootIoSvc::getTree(const std::string &type) { 
    RootOutputDesc* outputDesc = getRootOutputDesc(type);

    if (outputDesc) {
        return outputDesc->getTree();
    } 
    MsgStream log (msgSvc(), name() );
    log << MSG::WARNING << "RootOutputDesc of type " << type << " not found" << endreq;
    return 0;
 }

 StatusCode RootIoSvc::setupBranch(const std::string& type, const std::string &bname, 
     const std::string &classname, void* branchAddr,
     int bufSize, int splitLevel) {
    RootOutputDesc* outputDesc = getRootOutputDesc(type);

    if (outputDesc) {
        bool stat = outputDesc->setupBranch(bname, classname, branchAddr, bufSize, splitLevel);
        return ((stat == true) ? StatusCode::SUCCESS : StatusCode::FAILURE);
    } else {
        MsgStream log (msgSvc(), name() );
        log << MSG::WARNING << "RootOutputDesc of type " << type << " not found" << endreq;
    }
    return StatusCode::FAILURE;

}

StatusCode RootIoSvc::closeFile( const std::string & type ) {
    RootOutputDesc* outputDesc = getRootOutputDesc(type);

    if (outputDesc) {
        bool stat = outputDesc->closeFile();
        return (stat == true) ? StatusCode::SUCCESS : StatusCode::FAILURE;
    } else {
        MsgStream log (msgSvc(), name() );
        log << MSG::WARNING << "RootOutputDesc of type " << type << " not found" << endreq;
    }
    return StatusCode::FAILURE;
}
     
StatusCode RootIoSvc::fillTree(const std::string &type) {
    RootOutputDesc* outputDesc = getRootOutputDesc(type);

    if (outputDesc) {
        /*bool status =*/ outputDesc->fillTree(this->getAutoSaveInterval());
        outputDesc->setUpdated(true);        
        return StatusCode::SUCCESS;
    } 
    MsgStream log (msgSvc(), name() );
    log << MSG::WARNING << "RootOutputDesc of type " << type << " not found" << endreq;
    return StatusCode::FAILURE;

}


// finalize
StatusCode RootIoSvc::finalize ()
{
    StatusCode  status = StatusCode::SUCCESS;
/*

    if (m_rootTupleSvc) {
        void *treePtr;
        m_rootTupleSvc->getOutputTreePtr(treePtr,"MeritTuple");
        TTree *theTree = (TTree*)treePtr;
        if (theTree) { 
            FileHeader *meritHeader = m_headersTool->meritHeader();
            if (meritHeader) {
                meritHeader->setInteger("MeritVersion", m_rootTupleSvc->getMeritVersion());
                m_headersTool->writeMeritHeader(theTree->GetCurrentFile());
            } else {
                MsgStream log(msgSvc(),name()) ;
                log << MSG::WARNING << "MeritHeader object is NULL, not writing"
                    << " Merit Header, even though an output file exists"
                    << endreq;
            }
        }

    }
*/

    if (!m_celFileNameWrite.empty())
     {
      //m_celManager.fillFileAndTreeSet();
      TDirectory * saveDir = gDirectory ;
      long long retVal = m_outputCel->fillFileAndTreeSet() ;
      if (retVal < 0) {
          MsgStream log(msgSvc(),name()) ;
          log << MSG::WARNING << "Error Condition when closing CompositeEventList via fillFileAndTreeSet" << endreq;
      }
      m_outputCel->writeAndClose() ;
      delete m_outputCel ;
      m_outputCel = 0 ;
      saveDir->cd() ;
     }
    
    return status;
}

/// Query interface
StatusCode RootIoSvc::queryInterface(const InterfaceID& riid, void** ppvInterface)  {
    if ( IRootIoSvc::interfaceID() == riid )  {
        *ppvInterface = (IRootIoSvc*)this;
    }else if (IRunable::interfaceID()==riid ) {
      *ppvInterface = (IRunable*)this;
	} else if (IIncidentListener::interfaceID()==riid ) {
		*ppvInterface = (IIncidentListener*)this;
	} else  {
        return Service::queryInterface(riid, ppvInterface);
    }

    addRef();
    return SUCCESS;
}


void RootIoSvc::setEvtMax(Long64_t max)
 {
  // Purpose and Method:  Allow users of the RootIoSvc to specify the number
  //  of events found in their ROOT files
  if ( m_rootEvtMax == 0 )
   {
    m_rootEvtMax = max;
    return ;
   }   
  if ( m_rootEvtMax > max )
   { m_rootEvtMax = max ; }
 }

 bool RootIoSvc::setIndex( Long64_t i )
 {
     // Purpose and Method:  Interface routine for the event display to allow GUI to check on the
     // availability of an index into the ROOT files being read in.
     // As long as any input ROOT file has the index, we will report true.

     bool status = false;
     // In the context of an event collection, the index i will refer to the cel ROOT file.
     if (i<0) return false ;

     std::map<std::string,RootInputDesc *>::iterator typeItr ;
     for ( typeItr = m_rootIoMap.begin() ; typeItr != m_rootIoMap.end() ; ++typeItr )
     {
         RootInputDesc * rootInputDesc = typeItr->second ;
         if (!m_celFileNameRead.empty()) {  // access via event collection
             Long64_t ind = m_celManager.getEventIndexInTree(typeItr->first, i);
             if (ind < 0) 
                 status |= false;
             else
                 status |= true;
         } else {
             status |= rootInputDesc->checkEventAvailability(i);
         }
     }

     if (!status) return false;

     m_index = i ;
     m_runEventPair = std::pair<int,int>(-1,-1) ;
     m_useIndex = true ;
     return true ;
 }


 bool RootIoSvc::setRunEventPair( std::pair<int, int> ids )
 {
     // Purpose and Method:  Interface routine for the event display.  Allows the client to check
     // if a run/event id pair is available for loading.  If at least one of the TTrees has the event
     // available we will return true, otherwise return false.

     MsgStream log( msgSvc(), name() );
     bool status = false;
     std::map<std::string,RootInputDesc *>::iterator typeItr ;
     for ( typeItr = m_rootIoMap.begin() ; typeItr != m_rootIoMap.end() ; ++typeItr )
     {
         RootInputDesc * rootInputDesc = typeItr->second ;
         status |= rootInputDesc->checkEventAvailability(ids.first, ids.second);
     }

     if (!status) {
        return false;
     }
     m_runEventPair = ids;
     m_useIndex = false;

     return true;
 }


//==============================================================
//
// Gaudi machinery
//
//==============================================================

// handle "incidents"
void RootIoSvc::handle(const Incident &inc)
{
    if( inc.type()=="BeginEvent")beginEvent();
    else if(inc.type()=="EndEvent")endEvent();
}


void RootIoSvc::beginEvent() // should be called at the beginning of an event
{ 
    m_objectNumber = TProcessID::GetObjectCount();
    //static bool gotTuple = false;

    //if ((!m_celFileNameWrite.empty()) && (!gotTuple) && (m_rootTupleSvc)) {
    //    try {
        // Use this to recover a pointer to the output tree
   //     void * ptr= 0;
   //     m_rootTupleSvc->getItem(m_tupleName.value().c_str(),"", ptr);
   //     if( ptr==0) {
   //         MsgStream log( msgSvc(), name() );
   //         log << MSG::WARNING << "Could not find the tuple: " +  m_tupleName.value() << endreq;
    //    } else {
    //        m_outputTuple = static_cast<TTree*>(ptr);
            // Include the merit ntuple in the output CEL
   //         m_outputCel->declareComponent("merit") ;
    //        m_celTreeCol.push_back(m_outputTuple);
    //    }
    //    gotTuple = true;
    //    } catch(...) {
    //        m_outputTuple = 0;
    //    }
 //   }


    if (!m_celFileNameWrite.empty())
     {
      setOutputUpdate(false) ;
     }
}

void RootIoSvc::endEvent()  // must be called at the end of an event to update, allow pause
 {        
    m_useIndex = false;
    // temp
    m_useIndex = true ;
    m_index++;

    m_runEventPair = std::pair<int, int>(-1,-1) ;

    // clear out the static maps
    m_common.clear();

    // reset object in order to avoid memleak
    TProcessID::SetObjectCount(m_objectNumber) ;


    // Fill Composite Event List if requested
    if ( (!m_celFileNameWrite.empty()) && (checkOutputUpdate()) )
     {
//    [David] The difference between trees in read mode and write mode
//    is now taken in charge bt the CEL, which will use "GetEntries()-1"
//    instead of "GetReadEntry()" when the TDirectory associated to a tree
//    "IsWritable()".
//      // [David] During a writing job, filling a tree apparently
//      // does not upgrade the value which is returned bt GetReadEntry(),
//      // and that we use to generate the CEL. That's why we make
//      // below the call to LoadTree.
//      std::vector<TTree*>::iterator treeItr ;
//      Long64_t numBytes ;
//      for ( treeItr=m_celTreeCol.begin() ; treeItr != m_celTreeCol.end() ; treeItr++ )
//       { numBytes = (*treeItr)->LoadTree((*treeItr)->GetEntries()-1) ; }
      m_outputCel->fillEvent(m_celTreeCol) ;
     }

 }

// Purpose and Method:  Control the event loop
StatusCode RootIoSvc::run()
 {
  StatusCode status = StatusCode::FAILURE ;
  MsgStream log(msgSvc(),name()) ;

    if ( 0 == m_appMgrUI )  return status; 

    // Check for an input merit tuple, and add its number of events to our
    // m_rootEvtMax
    if (m_rootTupleSvc) {
        void *treePtr;
        Long64_t entries = m_rootTupleSvc->getInputTreePtr(treePtr);
        if (entries > 0) setEvtMax(entries);
    }

    IProperty* propMgr=0;
    status = serviceLocator()->service("ApplicationMgr", propMgr );
    if( status.isFailure()) {
        log << MSG::ERROR << "Unable to locate PropertyManager Service" << endreq;
        return status;
    }

    // Pick up ApplicationMgr.EvtMax
    IntegerProperty evtMax("EvtMax",0);
    status = propMgr->getProperty( &evtMax );
    if (status.isFailure()) return status;
 
    // Determine if the min number of ROOT events is less than the
    // requested number of events in the jobOptions file
    IntegerProperty rootEvtMax("EvtMax", m_rootEvtMax);
    if ( (static_cast<int>(rootEvtMax-m_startIndex) < evtMax) || (evtMax==0) ) {
       evtMax = rootEvtMax - m_startIndex;
       setProperty(evtMax);
    } else setProperty(evtMax);

    m_evtMax = evtMax;

    // now find the top alg so we can monitor its error count
    IAlgManager* theAlgMgr =0 ;
    status = serviceLocator( )->getService( "ApplicationMgr",
        IAlgManager::interfaceID(),
        (IInterface*&)theAlgMgr );
    IAlgorithm* theIAlg;
    Algorithm*  theAlgorithm=0;
    IntegerProperty errorProperty("ErrorCount",0);
    
    status = theAlgMgr->getAlgorithm( "Top", theIAlg );
    if ( status.isSuccess( ) ) {
        try{
            theAlgorithm = dynamic_cast<Algorithm*>(theIAlg);
        } catch(...){
            status = StatusCode::FAILURE;
        }
    }
    if ( status.isFailure( ) ) {
        log << MSG::WARNING << "Could not find algorithm 'Top'; will not monitor errors" << endreq;
    }
    
    
    // loop over the events

    long long eventNumber= 0;
    double currentTime=m_startTime;
    
    { bool noend=true;
    log << MSG::INFO << "Runable interface starting event loop as :" ; 
    if( m_evtMax>0)  { log << " MaxEvt = " << m_evtMax; noend=false;  }
    if( m_endTime>0) { log << " EndTime= " << m_endTime; noend=false; }
    log << endreq;
    
    if(noend) { 
        log << MSG::WARNING << "No valid end condition specified: will not "
            << "process any events! MaxEvt = " << m_evtMax << endreq; 
    }
    }


     // Check for ROOT errors once before event loop in case there are 	 
     // errors to report 	 
     if ( (m_rootFileMessageHandler->getError()) || 	 
          (m_rootFileMessageHandler->getSysError()) || 	 
          (m_rootFileMessageHandler->getFatal()) ) { 	 
         m_rootFileMessageHandler->Print(); 	 
         if (m_abortOnRootError) { 	 
             // Bailing out and returning error code 	 
             log << MSG::ERROR << "Terminating due to ROOT error" << endreq; 	
             exit(-1); 	 
 	  } 	       
          m_rootFileMessageHandler->resetFlags();
     }


    // Not yet using time as a control on the event loop for ROOT
    while( m_evtMax>0 && eventNumber < m_evtMax
        || m_endTime>0 && currentTime< m_endTime ) {
        
        loopStatus( eventNumber, currentTime, log);

        status =  m_appMgrUI->nextEvent(1); // currently, always success
        
        // the single event may have created a failure. 
        // Check the ErrorCount propery of the Top alg.
        if( theAlgorithm !=0) theAlgorithm->getProperty(&errorProperty);
        if( status.isFailure() || errorProperty.value() > 0){
            status = StatusCode::FAILURE;
        }
        
        if( status.isFailure()) break;
        //if(flux!=0){
         //   currentTime = flux->gpsTime();
       // }


        // Check for ROOT errors
        if ( (m_rootFileMessageHandler->getError()) || 
             (m_rootFileMessageHandler->getSysError()) ||
             (m_rootFileMessageHandler->getFatal()) ) {
             m_rootFileMessageHandler->Print();
             if (m_abortOnRootError) {
                 // Bailing out and returning error code
                 log << MSG::ERROR << "Terminating due to ROOT error" << endreq;
                 exit(-1);
             }
             m_rootFileMessageHandler->resetFlags();
         }
   
        eventNumber++;
    }
    if( status.isFailure()){
        log << MSG::ERROR << "Terminating RootIoSvc loop due to error" << endreq;
        
    }else if( m_endTime>0 && currentTime >= m_endTime ) {
        log << MSG::INFO << "Loop terminated by time " << endreq;
    }else {
        log << MSG::INFO << "Processing loop terminated by event count" << endreq;
    }
    m_rootFileMessageHandler->Print();
    return status;
}

// purpose: print periodic messages
void RootIoSvc::loopStatus( long long eventNumber, double currentTime, MsgStream & log )
 {
  static int last_fraction = 0 ;

  double efrac = ( (m_evtMax>0) ? (100.*eventNumber/m_evtMax) : (0.0) ) ;
  double tfrac = ( (m_endTime>0) ? (100.*(currentTime-m_startTime)/(m_endTime-m_startTime)) : (0.0) ) ;

  int percent_complete= static_cast<int>(std::max(efrac,tfrac)) ;
  if (percent_complete!=last_fraction)
   {
    last_fraction = percent_complete ;
    if ((percent_complete<10)||((percent_complete%10)==0))
     {
      facilities::Timestamp tstamp;
      ios_base::fmtflags ioFlags = 0;
      ioFlags = 4096;
      log<<MSG::INFO
        << " [" << tstamp.getString() << "] "
        << std::setprecision(12) << std::resetiosflags(ioFlags)
        <<percent_complete<<"% complete: "
        <<" event "<<eventNumber<<",  time= "<<currentTime<<endreq ;
     }
   }
 }


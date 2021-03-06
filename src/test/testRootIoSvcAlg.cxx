#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"
#include "RootIo/IRootIoSvc.h"

/** @class testRootIoSvcAlg
 * @brief Takes data from the TDS to test reading from ROOT files
 *
 * @author Heather Kelly
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/RootIo/src/test/testRootIoSvcAlg.cxx,v 1.8.246.1 2010/10/08 16:40:46 heather Exp $
 */

class testRootIoSvcAlg : public Algorithm
{	
public:
    
    testRootIoSvcAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    StatusCode initialize();
   
    StatusCode execute();
    
    StatusCode finalize();   
        
private:
    IRootIoSvc*   m_rootIoSvc;
 
};

//static const AlgFactory<testRootIoSvcAlg>  Factory;
//const IAlgFactory& testRootIoSvcAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(testRootIoSvcAlg);

testRootIoSvcAlg::testRootIoSvcAlg(const std::string& name, ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator)
{
}

StatusCode testRootIoSvcAlg::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    setProperties();

    if ( service("RootIoSvc", m_rootIoSvc).isFailure() ){
        log << MSG::ERROR << "Couldn't find the RootIoSvc!" << endreq;
        return StatusCode::FAILURE;
    }    
	
	return sc;
    
}

StatusCode testRootIoSvcAlg::execute()
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    static int flag = 0;
	
    if (flag % 2) {
       bool good = m_rootIoSvc->setRunEventPair(std::pair<int,int>(20,400));
       if (good) {
          log << MSG::INFO << "Failed test, Found run/event 20,400 - which does not exist" << endreq;
       } else {
          log << MSG::INFO << "Passed test, could not find run/event 20,400" << endreq;
       }
       m_rootIoSvc->setIndex(flag);
    } else {
       bool stat = m_rootIoSvc->setRunEventPair(std::pair<int,int>(10,2));
       if (stat)
           log << MSG::INFO << "Requesting run/event (10,2) randomly" << endreq;
    }

    if (flag == 3) {
        bool retVal = m_rootIoSvc->setRootFile("","noFile.root","","");
        if (retVal == false) {
            log << MSG::INFO << "Passed ROOT file open test "
                << "failed to open non-existant file" 
                << endreq;
        }
        retVal = m_rootIoSvc->setRootFile( 
                 "",
                 "$(ROOTTESTDATADATAPATH)/vertical_surface_muons/digi.root", 
                 "$(ROOTTESTDATADATAPATH)//vertical_surface_muons/recon.root",
                 "");
        if (retVal)
             log << MSG::INFO << "Passed ROOT file open test "
                 << "Succeeded in opening new files" << endreq;
             bool good = m_rootIoSvc->setRunEventPair(std::pair<int,int>(10,1));
             if (good) 
                 log << MSG::INFO << "Requesting event 10,1 in new files" << endreq;
        else {
             log << MSG::INFO << "Failed to open new files" << endreq;
             return StatusCode::FAILURE;
        }
    }

    flag++;


    return sc;
}


StatusCode testRootIoSvcAlg::finalize()
{    
    StatusCode sc = StatusCode::SUCCESS;
    return sc;
}


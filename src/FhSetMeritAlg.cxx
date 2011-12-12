#include "RootIo/FhTool.h"
#include "ntupleWriterSvc/INTupleWriterSvc.h"
#include "TTree.h"

#include "GaudiKernel/Algorithm.h"
#include <GaudiKernel/MsgStream.h>
#include <GaudiKernel/AlgFactory.h>
#include <GaudiKernel/IService.h>
#include <GaudiKernel/IJobOptionsSvc.h>
#include <GaudiKernel/ISvcLocator.h>
#include <GaudiKernel/IDataProviderSvc.h>
#include <GaudiKernel/SmartDataPtr.h>


class FhSetMeritAlg : public Algorithm {	

public:

    FhSetMeritAlg( const std::string& name, ISvcLocator* pSvcLocator ) ;
    
    ///
    StatusCode initialize();
   
    ///
    StatusCode execute();
    
    /// prepare the common attributes and set them in available headers
    StatusCode finalize();

private:
    
} ;



//static const AlgFactory<FhSetMeritAlg> Factory ;
//const IAlgFactory& FhSetMeritAlgFactory = Factory ;
DECLARE_ALGORITHM_FACTORY(FhSetMeritAlg);

StatusCode FhSetMeritAlg::initialize() {
    MsgStream log(msgSvc(),name()) ;
    StatusCode sc ;
    
    // retrieve the headers
    IFhTool * headersTool ;
    sc = toolSvc()->retrieveTool("FhTool",headersTool) ;
    if (sc.isFailure()) {
        log<<MSG::WARNING << "Failed to retrieve headers tool" << endreq ;
        return StatusCode::FAILURE ;
    }

	return StatusCode::SUCCESS ;
}
    
FhSetMeritAlg::FhSetMeritAlg
( const std::string& name, ISvcLocator* pSvcLocator)
:  Algorithm(name, pSvcLocator) {
}

StatusCode FhSetMeritAlg::finalize() {
    StatusCode sc ;
    MsgStream log(msgSvc(),name()) ;
    // retrieve the headers
    IFhTool * m_headersTool ;
    sc = toolSvc()->retrieveTool("FhTool",m_headersTool) ;
    if (sc.isFailure()) {
        log << MSG::WARNING << "Failed to retrieve headers tool" << endreq ;
        return StatusCode::FAILURE ;
    }

    FileHeader * meritHeader = m_headersTool->meritHeader() ;

    // Retrieve tuple info an set up
    // get a pointer to RootTupleSvc
    INTupleWriterSvc *   m_rootTupleSvc;
    sc = service("RootTupleSvc", m_rootTupleSvc, true);       
    if( sc.isFailure() ) 
    {
        MsgStream log( msgSvc(), name() );
        log << MSG::WARNING << "failed to get the RootTupleSvc" << endreq;
        m_rootTupleSvc = 0;
    }

    if (m_rootTupleSvc) {
        void *treePtr;
        m_rootTupleSvc->getOutputTreePtr(treePtr,"MeritTuple");
        TTree *theTree = (TTree*)treePtr;
        if (theTree) { 
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

	return StatusCode::SUCCESS ;
}

StatusCode FhSetMeritAlg::execute() {
	return StatusCode::SUCCESS ;
}


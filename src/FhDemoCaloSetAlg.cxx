
#include "src/FhDemoCaloSetAlg.h"
#include "RootIo/FhTool.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/DigiEvent.h"

#include <GaudiKernel/MsgStream.h>
#include <GaudiKernel/AlgFactory.h>
#include <GaudiKernel/IService.h>
#include <GaudiKernel/ISvcLocator.h>
#include <GaudiKernel/IDataProviderSvc.h>
#include <GaudiKernel/SmartDataPtr.h>

#include <cstdio>

//static const AlgFactory<FhDemoCaloSetAlg> Factory ;
//const IAlgFactory& FhDemoCaloSetAlgFactory = Factory ;
DECLARE_ALGORITHM_FACTORY(FhDemoCaloSetAlg);

StatusCode FhDemoCaloSetAlg::finalize() {
    
    MsgStream log(msgSvc(),name()) ;
    StatusCode sc ;
    
    // retrieve the headers
    IFhTool * headersTool ;
    sc = toolSvc()->retrieveTool("FhTool",headersTool) ;
    if (sc.isFailure()) {
        log<<MSG::WARNING << "Failed to retrieve headers tool" << endreq ;
        return StatusCode::FAILURE ;
    }
    FileHeader * mcHeader = headersTool->mcHeader() ;
    FileHeader * digiHeader = headersTool->digiHeader() ;
    FileHeader * reconHeader = headersTool->reconHeader() ;
    
    // cal string
    if (digiHeader) {
        digiHeader->setString("CaloString","Demo String") ;
    }
    
    // end
	return StatusCode::SUCCESS ;
}
    
FhDemoCaloSetAlg::FhDemoCaloSetAlg
( const std::string& name, ISvcLocator* pSvcLocator)
:  Algorithm(name, pSvcLocator) {
}

StatusCode FhDemoCaloSetAlg::initialize() {
	return StatusCode::SUCCESS ;
}

StatusCode FhDemoCaloSetAlg::execute() {
	return StatusCode::SUCCESS ;
}
    


#include "src/FhSetAlg.h"
#include "RootIo/FhTool.h"
#include "RootIo/FhSystemEnv.h"
#include "RootIo/FhCmtConfig.h"
#include "RootIo/FhJobOptions.h"
#include "RootIo/FhChronoStatTable.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/DigiEvent.h"

#include <GaudiKernel/MsgStream.h>
#include <GaudiKernel/AlgFactory.h>
#include <GaudiKernel/IService.h>
#include <GaudiKernel/IJobOptionsSvc.h>
#include <GaudiKernel/ISvcLocator.h>
#include <GaudiKernel/IDataProviderSvc.h>
#include <GaudiKernel/SmartDataPtr.h>

#include <vector>
#include <cstdio>

#ifdef DEFECT_NO_STRINGSTREAM
#include <strstream>
#else
#include <sstream>
#endif

static const AlgFactory<FhSetAlg> Factory ;
const IAlgFactory& FhSetAlgFactory = Factory ;

StatusCode FhSetAlg::finalize() {
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
    
  
    // system environment variables
    FhSystemEnv systemEnv ;
    systemEnv.init() ;
    if (mcHeader) systemEnv.store(mcHeader) ;
    if (digiHeader) systemEnv.store(digiHeader) ;
    if (reconHeader) systemEnv.store(reconHeader) ;
    
    // cmt
    FhCmtConfig cmtConfig ;
    cmtConfig.init() ;
    if (mcHeader) cmtConfig.store(mcHeader) ;
    if (digiHeader) cmtConfig.store(digiHeader) ;
    if (reconHeader) cmtConfig.store(reconHeader) ;

    // job options
    FhJobOptions jobOptions ;
    jobOptions.init(serviceLocator()) ;
    if (mcHeader) jobOptions.store(mcHeader) ;
    if (digiHeader) jobOptions.store(digiHeader) ;
    if (reconHeader) jobOptions.store(reconHeader) ;

    // catch chrono stat summary
    FhChronoStatTable chronoStatTable ;
    chronoStatTable.init(serviceLocator()) ;
    if (mcHeader) chronoStatTable.store(mcHeader) ;
    if (digiHeader) chronoStatTable.store(digiHeader) ;
    if (reconHeader) chronoStatTable.store(reconHeader) ;

    log << MSG::INFO << "headers common data ready" << endreq ;
	return StatusCode::SUCCESS ;
}
    
FhSetAlg::FhSetAlg
( const std::string& name, ISvcLocator* pSvcLocator)
:  Algorithm(name, pSvcLocator) {
}

StatusCode FhSetAlg::initialize() {
	return StatusCode::SUCCESS ;
}

StatusCode FhSetAlg::execute() {
	return StatusCode::SUCCESS ;
}


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

StatusCode FhSetAlg::initialize() {
    MsgStream log(msgSvc(),name()) ;
    StatusCode sc ;
    
    // retrieve the headers
    IFhTool * headersTool ;
    sc = toolSvc()->retrieveTool("FhTool",headersTool) ;
    if (sc.isFailure()) {
        log<<MSG::WARNING << "Failed to retrieve headers tool" << endreq ;
        return StatusCode::FAILURE ;
    }

    // create the headers here so we can use the initialize method for the
    //  pipes
    headersTool->newMcHeader();
    headersTool->newDigiHeader();;
    headersTool->newReconHeader();
    headersTool->newGcrHeader();
    headersTool->newMeritHeader();

    FileHeader * mcHeader = headersTool->mcHeader() ;
    FileHeader * digiHeader = headersTool->digiHeader() ;
    FileHeader * reconHeader = headersTool->reconHeader() ;
    FileHeader * gcrHeader = headersTool->gcrHeader() ;
    FileHeader * meritHeader = headersTool->meritHeader() ;
    
  
    // system environment variables
    FhSystemEnv systemEnv ;
    int status = systemEnv.init() ;
    if (status < 0) {
      log << MSG::WARNING << "Failed to retrieve system env for ROOT headers" 
          << endreq;
    }else {
    if (mcHeader) systemEnv.store(mcHeader) ;
    if (digiHeader) systemEnv.store(digiHeader) ;
    if (reconHeader) systemEnv.store(reconHeader) ;
    if (gcrHeader) systemEnv.store(gcrHeader) ;
    if (meritHeader) systemEnv.store(meritHeader) ;
    }
    
    // cmt
    FhCmtConfig cmtConfig ;
    status = cmtConfig.init() ;
    if (status < 0) {
        log << MSG::WARNING << "Failed to retrieve CMT config for ROOT headers"
            << endreq;
    } else {
    if (mcHeader) cmtConfig.store(mcHeader) ;
    if (digiHeader) cmtConfig.store(digiHeader) ;
    if (reconHeader) cmtConfig.store(reconHeader) ;
    if (gcrHeader) cmtConfig.store(gcrHeader) ;
    if (meritHeader) cmtConfig.store(meritHeader) ;
    }

    // job options
    FhJobOptions jobOptions ;
    jobOptions.init(serviceLocator()) ;
    if (mcHeader) jobOptions.store(mcHeader) ;
    if (digiHeader) jobOptions.store(digiHeader) ;
    if (reconHeader) jobOptions.store(reconHeader) ;
    if (gcrHeader) jobOptions.store(gcrHeader) ;
    if (meritHeader) jobOptions.store(meritHeader) ;

    log << MSG::INFO << "headers common data ready" << endreq ;
	return StatusCode::SUCCESS ;
}
    
FhSetAlg::FhSetAlg
( const std::string& name, ISvcLocator* pSvcLocator)
:  Algorithm(name, pSvcLocator) {
}

StatusCode FhSetAlg::finalize() {
    StatusCode sc ;
    MsgStream log(msgSvc(),name()) ;
    // retrieve the headers
    IFhTool * headersTool ;
    sc = toolSvc()->retrieveTool("FhTool",headersTool) ;
    if (sc.isFailure()) {
        log << MSG::WARNING << "Failed to retrieve headers tool" << endreq ;
        return StatusCode::FAILURE ;
    }
    FileHeader * mcHeader = headersTool->mcHeader() ;
    FileHeader * digiHeader = headersTool->digiHeader() ;
    FileHeader * reconHeader = headersTool->reconHeader() ;
    FileHeader * gcrHeader = headersTool->gcrHeader() ;
    FileHeader * meritHeader = headersTool->meritHeader() ;
    // catch chrono stat SUMMARY
    FhChronoStatTable chronoStatTable ;
    chronoStatTable.init(serviceLocator()) ;
    if (mcHeader) chronoStatTable.store(mcHeader) ;
    if (digiHeader) chronoStatTable.store(digiHeader) ;
    if (reconHeader) chronoStatTable.store(reconHeader) ;
    if (gcrHeader) chronoStatTable.store(gcrHeader) ;
    if (meritHeader) chronoStatTable.store(meritHeader) ;
	return StatusCode::SUCCESS ;
}

StatusCode FhSetAlg::execute() {
	return StatusCode::SUCCESS ;
}

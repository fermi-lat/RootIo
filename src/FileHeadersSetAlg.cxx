
#include "src/FileHeadersSetAlg.h"
#include "src/FileHeadersTool.h"

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
#include <sstream>
#include <cstdio>

static const AlgFactory<FileHeadersSetAlg> Factory ;
const IAlgFactory& FileHeadersSetAlgFactory = Factory ;

StatusCode FileHeadersSetAlg::finalize() {

    MsgStream log(msgSvc(),name()) ;
    StatusCode sc ;
    
    // retrieve the headers
    IFileHeadersTool * headersTool ;
    sc = toolSvc()->retrieveTool("FileHeadersTool",headersTool) ;
    if (sc.isFailure()) {
        log<<MSG::WARNING << "Failed to retrieve headers tool" << endreq ;
        return StatusCode::FAILURE ;
    }
    FileHeader * mcHeader = headersTool->mcHeader() ;
    FileHeader * digiHeader = headersTool->digiHeader() ;
    FileHeader * reconHeader = headersTool->reconHeader() ;
    
    FILE * pipe ;
    char buffer[256] ;
  
//    better stored event by event  
//    // runId
//    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),EventModel::EventHeader);
//    if (mcHeader) mcHeader->setInteger("runId",eventHeader->run()) ;
//    if (digiHeader) digiHeader->setInteger("runId"eventHeader->run()) ;
//    if (reconHeader) reconHeader->setInteger("runId"eventHeader->run()) ;
    
    
    // unix environment variables
    TString unixEnv ;
    pipe = popen("printenv 2>&1","r") ;
    while ( fgets(buffer,256,pipe) != NULL )
     { unixEnv += buffer ; }
    pclose(pipe) ;
    if (mcHeader) mcHeader->setString("unixEnv",unixEnv) ;
    if (digiHeader) digiHeader->setString("unixEnv",unixEnv) ;
    if (reconHeader) reconHeader->setString("unixEnv",unixEnv) ;
    
    // David: for the cmt commands below, I must provide a package
    // name ; I would like to provide higher level package used for
    // the current appplication, but I found no way to guess it ; It
    // is why I am using RootIo. 
    
    // cmt show uses
    TString cmtUses ;
    pipe = popen("cmt -pack=RootIo show uses 2>&1","r") ;
    while ( fgets(buffer,256,pipe) != NULL )
     { cmtUses += buffer ; }
    pclose(pipe) ;
    if (mcHeader) mcHeader->setString("cmtUses",cmtUses) ;
    if (digiHeader) digiHeader->setString("cmtUses",cmtUses) ;
    if (reconHeader) reconHeader->setString("cmtUses",cmtUses) ;
        
    // cmt show macros
    TString cmtMacros ;
    pipe = popen("cmt -pack=RootIo show macros 2>&1","r") ;
    while ( fgets(buffer,256,pipe) != NULL )
     { cmtMacros += buffer ; }
    pclose(pipe) ;
    if (mcHeader) mcHeader->setString("cmtMacros",cmtMacros) ;
    if (digiHeader) digiHeader->setString("cmtMacros",cmtMacros) ;
    if (reconHeader) reconHeader->setString("cmtMacros",cmtMacros) ;
    
    // job options
    // I try to generate a format close to the options file,
    // but easily readable in C++, which means removing the useless
    // whitespaces
    IService * joSvc ;
    IJobOptionsSvc * joInt ;
    StatusCode joSc ;
    joSc = serviceLocator()->getService("JobOptionsSvc",joSvc,true) ;
    if (joSc.isSuccess()) {
        joSc = joSvc->queryInterface(IID_IJobOptionsSvc,(void**)&joInt) ;
    }
    if (joSc.isFailure()) {
        log << MSG::WARNING << "Failed to retrieve job options service" << endreq ;
    }
    else {
        std::ostringstream jobOptions ;
        std::vector< std::string > clients = joInt->getClients() ;
        std::vector< std::string >::iterator client ;
        std::string name, value ;
        for ( client = clients.begin() ; client != clients.end() ; ++client ) {
            const std::vector< const Property * > * properties = joInt->getProperties(*client) ;  
            std::vector< const Property * >::const_iterator property ;
            for ( property = properties->begin() ; property != properties->end() ; ++property ) {
                name = (*client) ;
                name.append(".") ;
                name.append((*property)->name()) ;
                jobOptions<<name<<" " ;
                std::ostringstream os ;
                (*property)->nameAndValueAsStream(os) ;
                os<<std::ends ;
                value = os.str() ;
                std::string::size_type pos = value.find(":") ;
                value.erase(0,pos+2) ;
                const char * c = value.data() ;
                bool closed = true ;
                while ((*c)!=0) {
                    if ((*c)=='"') {
                        if ((*(c+1))!='"') {
                            closed = ! closed ;
                        }
                        else {
                            jobOptions<<*c ;
                        }
                    }
                    else if ((*c)==' ') {
                        if (!closed) {
                            jobOptions<<'_' ;
                        }
                    }
                    else if ((*c)=='[') {
                        if (!closed) {
                            jobOptions<<*c ;
                        }
                        else {
                            jobOptions<<'{' ;
                        }
                    }
                    else if ((*c)==']') {
                        if (!closed) {
                            jobOptions<<*c ;
                        }
                        else {
                            jobOptions<<'}' ;
                        }
                    }
                    else {
                        jobOptions<<*c ;
                    }
                    ++c ;
                }
                jobOptions<<'\n' ;
            }
        }
        jobOptions<<std::ends ;
        std::string result = jobOptions.str() ;
        if (mcHeader) mcHeader->setString("jobOptions",result.c_str()) ;
        if (digiHeader) digiHeader->setString("jobOptions",result.c_str()) ;
        if (reconHeader) reconHeader->setString("jobOptions",result.c_str()) ;
    }

    // catch chrono stat summary
    std::ostringstream mylog ;
    IService * theChronoSvc ;
    StatusCode theChronoSc ;
    theChronoSc = serviceLocator()->getService("ChronoStatSvc",theChronoSvc,false) ;
    if (theChronoSc==StatusCode::SUCCESS) {
        std::ostream * os = msgSvc()->defaultStream() ;
        msgSvc()->setDefaultStream(&mylog) ;
        theChronoSvc->finalize() ;
        msgSvc()->setDefaultStream(os) ;
        std::string result = mylog.str() ;
        if (mcHeader) mcHeader->setString("chronoStatTable",result.c_str()) ;
        if (digiHeader) digiHeader->setString("chronoStatTable",result.c_str()) ;
        if (reconHeader) reconHeader->setString("chronoStatTable",result.c_str()) ;
    }

    log << MSG::INFO << "headers common data ready" << endreq ;
	return StatusCode::SUCCESS ;
}
    
FileHeadersSetAlg::FileHeadersSetAlg
( const std::string& name, ISvcLocator* pSvcLocator)
:  Algorithm(name, pSvcLocator) {
}

StatusCode FileHeadersSetAlg::initialize() {
	return StatusCode::SUCCESS ;
}

StatusCode FileHeadersSetAlg::execute() {
	return StatusCode::SUCCESS ;
}
FILE* FileHeadersSetAlg::popen(const char* command, const char *mode) {
#ifdef WIN32
	return _popen(command,mode);
#else
	return popen(command,mode);
#endif
}

int FileHeadersSetAlg::pclose(FILE* fp) {
#ifdef WIN32
	return _pclose(fp);
#else
	return pclose(fp);
#endif
}


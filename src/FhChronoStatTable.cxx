
#include "RootIo/FhChronoStatTable.h"
#include <GaudiKernel/IService.h>
#include <GaudiKernel/IMessageSvc.h>
#include <sstream>
#include <string>

FhChronoStatTable::FhChronoStatTable() {
}

void FhChronoStatTable::init( const FileHeader * header ) {
    if (!header) { m_raw = "" ; }
    else { m_raw = header->getString("ChronoStatTable") ; }
    rawToMap() ;
}

void FhChronoStatTable::init( ISvcLocator * svcLocator ) {
    m_raw = "" ;
    IService * theChronoSvc, * aSvc ;
    IMessageSvc * theMessageSvc ;
    StatusCode theChronoSc, theMessageSc ;
    theMessageSc = svcLocator->getService("MessageSvc",aSvc,false) ;
    if (theMessageSc==StatusCode::SUCCESS) {
        theMessageSc = aSvc->queryInterface(IID_IMessageSvc,
            (void**)&theMessageSvc) ;
    }
    theChronoSc = svcLocator->getService("ChronoStatSvc",theChronoSvc,false) ;
    if ((theChronoSc==StatusCode::SUCCESS)&&
        (theMessageSc==StatusCode::SUCCESS)) {
        std::ostringstream mylog ;
        std::ostream * os = theMessageSvc->defaultStream() ;
        theMessageSvc->setDefaultStream(&mylog) ;
        theChronoSvc->finalize() ;
        theMessageSvc->setDefaultStream(os) ;
        m_raw = mylog.str().c_str() ;
    }
    rawToMap() ;
}

void FhChronoStatTable::store( FileHeader * header ) const {
        if (header) header->setString("ChronoStatTable",m_raw) ;
}

void FhChronoStatTable::rawToMap() {
    m_map.DeleteAll() ;
    std::istringstream chronoStatTable(m_raw.Data()) ;
    std::string words[7], value, unit, name ;
    int iword = 0 ;
    while (chronoStatTable>>words[iword]) {
        if (words[iword]=="Tot=") {
            name = words[(iword-5+7)%7] ;
            chronoStatTable>>value>>unit ;
            value += " " ;
            value += unit ;
            m_map.add(name.c_str(),value.c_str()) ;
        }
        iword = (iword+1)%7 ;
    }
}

void FhChronoStatTable::getTagNames( TRegexp exp, TCollection & tagNames ) const {
    m_map.getNames(exp,tagNames) ;
}

TString FhChronoStatTable::getStringTime( const TObjString & tagName ) const {
    return m_map.getValue(tagName) ;
}

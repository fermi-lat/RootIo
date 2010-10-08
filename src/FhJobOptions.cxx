
#include "RootIo/FhJobOptions.h"
#include <GaudiKernel/IService.h>
#include <GaudiKernel/IJobOptionsSvc.h>
#include <GaudiKernel/Property.h>
#include <iostream>

#ifdef DEFECT_NO_STRINGSTREAM
#include <strstream>
#else
#include <sstream>
#endif

#include <string>

FhJobOptions::FhJobOptions() {
}

void FhJobOptions::init( const FileHeader * header ) {
    m_raw = "" ;
    if (header) {
        m_raw = header->getString("JobOptions") ;
    }
    rawToMap() ;
}

void FhJobOptions::init( ISvcLocator * svcLocator ) {
    m_raw = "" ;
    IService * joSvc =0 ;
    IJobOptionsSvc * joInt =0 ;
    StatusCode joSc ;
    joSc = svcLocator->getService("JobOptionsSvc",joSvc,true) ;
    if (joSc.isSuccess()) {
        joSc = joSvc->queryInterface(IJobOptionsSvc::interfaceID(),(void**)&joInt) ;
    }
    if (joSc.isSuccess()) {
        std::vector< std::string > clients = joInt->getClients() ;
        std::vector< std::string >::iterator client ;
        std::string name, value ;
        for ( client = clients.begin() ; client != clients.end() ; ++client ) {
            const std::vector< const Property * > * properties = joInt->getProperties(*client) ;  
            std::ostringstream jobOptions ;
            std::vector< const Property * >::const_iterator property ;
            for ( property = properties->begin() ; property != properties->end() ; ++property ) {
                // I insert the client name, in the same format as the property
                jobOptions<<"\""<<(*client)<<"\"." ;
                (*property)->fillStream(jobOptions) ;
                jobOptions<<'\n' ;
            }
            jobOptions<<std::ends ;
            m_raw += jobOptions.str().c_str() ;
        }
    }
    rawToMap() ;
}

void FhJobOptions::store( FileHeader * header ) const {
    static bool print=false;
    if (header) {
            header->setString("JobOptions",m_raw) ;
            if (!print) {
                std::cout<<"%%%%% JOB OPTIONS %%%%%\n"<<m_raw<<std::endl ;
                print = true;
            }
//        TOrdCollection names ;
//        m_raw.getNames(".*",names) ;
//        TIterator * iter = names.MakeIterator() ;
//        TObjString * key ;
//        TString name ;
//        while (( key = (TObjString *)(iter->Next()) )) {
//            name = "JobOption." ;
//            name += *key ;
//            header->setString(name.Data(),m_raw.GetValue(*key)) ;
//        }
//        delete iter ;
    }
}

void FhJobOptions::rawToMap() {
    // I try to generate a format close to the options file,
    // but easily readable in C++, which means removing the useless
    // whitespaces
    m_map.DeleteAll() ;
    std::istringstream is(m_raw.Data()) ;
    unsigned int buffer_size = m_raw.Length()+1  ;
    char * buffer = new char [buffer_size] ;
    std::string name, value ;
    const char * c1,* line2 ;
    char * c2 ;
    while (is.getline(buffer,buffer_size)) {
        std::string line1(buffer) ;
        c1 = line1.data() ;
        line2 = c2 = new char [(line1.size()+1)] ;
        bool closed = true ;
        while ((*c1)!=0) {
            if ((*c1)=='"') {
                if ((*(c1+1))!='"') {
                    closed = ! closed ;
                }
                else {
                    *c2++ = *c1 ;
                }
            }
            else if ((*c1)==' ') {
                if (!closed) {
                    *c2++ = '_' ;
                }
            }
            else if ((*c1)=='[') {
                if (!closed) {
                    *c2++ = *c1 ;
                }
                else {
                    *c2++ = '{' ;
                }
            }
            else if ((*c1)==']') {
                if (!closed) {
                    *c2++ = *c1 ;
                }
                else {
                    *c2++ = '}' ;
                }
            }
            else {
                *c2++ = *c1 ;
            }
            ++c1 ;
        }
        *c2++ = 0 ;
        line1 = line2 ;
        delete [] line2 ;
        std::string::size_type pos = line1.find(':') ;
        if (pos!=std::string::npos) {
            name = line1.substr(0,pos) ;
            value = line1.substr(pos+1,line1.length()-pos-1) ;            
            m_map.add(name.c_str(),value.c_str()) ;
        }
    }
    delete [] buffer ;
    buffer = 0 ;
}

void FhJobOptions::getNames( TRegexp exp, TCollection & names ) const {
    m_map.getNames(exp,names) ;
}

TString FhJobOptions::getValue( const TObjString & name ) const {
    return m_map.getValue(name) ;
}


#include "RootIo/FhSystemEnv.h"
#include <string>
#include <cstdio>

#ifdef DEFECT_NO_STRINGSTREAM
#include <strstream>
#else
#include <sstream>
#endif

FILE* fhOpenEnv() {
#ifdef WIN32
	return _popen("set","r");
#else
	return ::popen("printenv","r");
#endif
}

int fhClose(FILE* fp) {
#ifdef WIN32
	return _pclose(fp);
#else
	return ::pclose(fp);
#endif
}

FhSystemEnv::FhSystemEnv() {
}

void FhSystemEnv::init( const FileHeader * header ) {
    if (!header) { m_raw = "" ; }
    else { m_raw = header->getString("SystemEnv") ; }
    rawToMap() ;
}

void FhSystemEnv::init() {
    m_raw = "" ;
    FILE * pipe = fhOpenEnv() ;
    char buffer[256] ;
    while ( fgets(buffer,256,pipe) != NULL )
     { m_raw += buffer ; }
    fhClose(pipe) ;
    rawToMap() ;
}

void FhSystemEnv::store( FileHeader * header ) const {
    if (header) {
        header->setString("SystemEnv",m_raw) ;
    }
}

void FhSystemEnv::rawToMap() {
    m_map.DeleteAll() ;
    std::istringstream is(m_raw.Data()) ;
    unsigned int buffer_size = m_raw.Length()+1  ;
    char * buffer = new char [buffer_size] ;
    std::string name, value ;
    while (is.getline(buffer,buffer_size)) {
        std::string line(buffer) ;
        std::string::size_type pos = line.find('=') ;
        if (pos!=std::string::npos) {
            name = line.substr(0,pos) ;
            value = line.substr(pos+1,line.length()-pos-1) ;            
            m_map.add(name.c_str(),value.c_str()) ;
        }
    }
    delete [] buffer ;
    buffer = 0 ;
}

void FhSystemEnv::getVariableNames( TRegexp exp, TCollection & names ) const {
    m_map.getNames(exp,names) ;
}

TString FhSystemEnv::getValue( const TObjString & name ) const {
    return m_map.getValue(name) ;
}


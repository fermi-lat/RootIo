
#include "RootIo/FhCmtConfig.h"

#ifdef DEFECT_NO_STRINGSTREAM
#include <strstream>
#else
#include <sstream>
#endif

#include <string>
#include <cstdio>

FILE* fhPopen(const char* command, const char *mode) {
#ifdef WIN32
	return _popen(command,mode);
#else
	return ::popen(command,mode);
#endif
}

int fhPclose(FILE* fp) {
#ifdef WIN32
	return _pclose(fp);
#else
	return ::pclose(fp);
#endif
}

FhCmtConfig::FhCmtConfig() {
}

void FhCmtConfig::init( const FileHeader * header ) {
    if (!header) { m_rawUses = m_rawMacros = "" ; }
    else {
         m_rawUses = header->getString("CmtUses") ;
         m_rawMacros = header->getString("CmtMacros") ;
    }
    rawToMap() ;
}

void FhCmtConfig::init() {
    
    // David: for the cmt commands below, I must provide a package
    // name ; I would like to provide higher level package used for
    // the current appplication, but I found no way to guess it ; It
    // is why I am using RootIo. 
    
    m_rawUses = "" ;
    FILE * pipe = fhPopen("cmt -pack=RootIo show uses 2>&1","r") ;
    char buffer[256] ;
    while ( fgets(buffer,256,pipe) != NULL )
     { m_rawUses += buffer ; }
    fhPclose(pipe) ;
       
    m_rawMacros = "" ;
    pipe = fhPopen("cmt -pack=RootIo show macros 2>&1","r") ;
    while ( fgets(buffer,256,pipe) != NULL )
     { m_rawMacros += buffer ; }
    fhPclose(pipe) ;
              
    rawToMap() ;
}

void FhCmtConfig::store( FileHeader * header ) const {
        if (header) {
            header->setString("CmtUses",m_rawUses) ;
            header->setString("CmtMacros",m_rawMacros) ;
        }
}

void FhCmtConfig::rawToMap() {
    
    unsigned int buffer_size  ;
    char * buffer  ;

    // cmt show uses
    m_mapUses.DeleteAll() ;
    m_useHierarchy = "" ;
    std::istringstream isUses(m_rawUses.Data()) ;
    buffer_size = m_rawUses.Length()+1 ;
    buffer = new char [buffer_size] ;
    int step = 0 ;
    TRegexp expHierarchy("^# *use") ;
    TRegexp expComment("^#") ;
    while (isUses.getline(buffer,buffer_size)) {
        TString lineUses(buffer) ;
        if ((step==0)&&(lineUses.Index(expHierarchy)>0)) {
            lineUses.Remove(0,2) ;
            m_useHierarchy += lineUses ;
        } else if ((step<2)&&(lineUses.Index(expComment)>0)) {
            step = 1 ;
        } else {
            step = 2 ;
            std::istringstream cmtUse(buffer) ;
            std::string name, value, tmp ;
            cmtUse>>tmp>>name>>value ;
            while (cmtUse>>tmp) {
              value += " " ;
              value += tmp ;
            }
            m_mapUses.add(name.c_str(),value.c_str()) ;
        }       
    }
    delete [] buffer ;
    buffer = 0 ;

    // cmt show macros
    m_mapMacros.DeleteAll() ;
    std::istringstream isMacros(m_rawMacros.Data()) ;
    buffer_size = m_rawMacros.Length()+1 ;
    buffer = new char [buffer_size] ;
    while (isMacros.getline(buffer,buffer_size)) {
        std::string cmtMacro(buffer) ;
        std::string::size_type pos = cmtMacro.find('=') ;
        if (pos!=std::string::npos) {
            std::string name = cmtMacro.substr(0,pos) ;
            std::string value = cmtMacro.substr(pos+2,cmtMacro.length()-pos-3) ;            
            m_mapMacros.add(name.c_str(),value.c_str()) ;
        }
    }
    delete [] buffer ;
    buffer = 0 ;
}

void FhCmtConfig::getPackageNames( TRegexp exp, TCollection & names ) const {
    m_mapUses.getNames(exp,names) ;
}

TString FhCmtConfig::getPackageSelection( const TObjString & name ) const {
    return m_mapUses.getValue(name) ;
}

TString FhCmtConfig::getUseHierarchy() const {
    return m_useHierarchy ;
}

void FhCmtConfig::getMacroNames( TRegexp exp, TCollection & names ) const {
    m_mapMacros.getNames(exp,names) ;
}

TString FhCmtConfig::getMacroValue( const TObjString & name ) const {
    return m_mapMacros.getValue(name) ;
}

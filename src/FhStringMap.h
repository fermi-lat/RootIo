
#ifndef FhStringMap_H
#define FhStringMap_H

#include <TRegexp.h>
#include <TObjString.h>
#include <TMap.h>
class ISvcLocator ;

/*! 

 @class FhStringMap
 @brief Set of dynamically named strings

 @author David Chamont - CNRS IN2P3 LLR Ecole Polytechnique

*/

class FhStringMap {

public:

    FhStringMap() ;
    ~FhStringMap() ;

    void DeleteAll() ;
 
    // David: I built the interface with TObjString for the names,
    // rather than char * or TString, so that the user will
    // implicitly know that the getNames methods is storing
    // TObjString in the result collection.
 
    void add( const TObjString & name, const TString & value ) ;
    TString getValue( const TObjString & name ) const ;
    void getNames( TRegexp exp, TCollection & result ) const ;

private:
    TMap m_strings ;
} ; 
 
#endif







#ifndef FhChronoStatTable_H
#define FhChronoStatTable_H

#include "src/FhStringMap.h"
#include "commonRootData/FileHeader.h"
#include <GaudiKernel/ISvcLocator.h>

/*!

 @class FhChronoStatTable   
 @brief Prepare/interpret a chrono stat table for/from a file header

 @author David Chamont - CNRS IN2P3 LLR Ecole Polytechnique

*/

class FhChronoStatTable {	

public:

    FhChronoStatTable() ;
    
    // global read/write methods
    void init( const FileHeader *  ) ;
    void init( ISvcLocator * ) ;
    void store( FileHeader * ) const ;
   
    // David: I keep an interface dependant on ROOT classes,
    // rather than std ones, because I want to
    // use TRegexp.
   
    // get values
    void getTagNames( TRegexp exp, TCollection & names ) const ;
    TString getStringTime( const TObjString & name ) const ;
    
private:

    void rawToMap() ;
	TString m_raw ;
	FhStringMap m_map ;
    
} ;

#endif



#ifndef FhJobOptions_H
#define FhJobOptions_H

#include "src/FhStringMap.h"
#include "commonRootData/FileHeader.h"
#include <GaudiKernel/ISvcLocator.h>
#include <string>

/*!

 @class FhJobOptions   
 @brief Prepare/interpret job options for/from a file header

 @author David Chamont - CNRS IN2P3 LLR Ecole Polytechnique

*/

class FhJobOptions {	

public:

    FhJobOptions() ;
    
    // global read/write methods
    void init( const FileHeader *  ) ;
    void init( ISvcLocator * ) ;
    void store( FileHeader * ) const ;
   
    // David: I keep an interface dependant on ROOT classes,
    // rather than std ones, because I want to
    // use TRegexp.
   
    // get values
    void getNames( TRegexp exp, TCollection & names ) const ;
    TString getValue( const TObjString & name ) const ;
    
private:

    void rawToMap() ;
	TString m_raw ;
	FhStringMap m_map ;
    
} ;

#endif


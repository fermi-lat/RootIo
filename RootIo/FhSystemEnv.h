
#ifndef FhSystemEnv_H
#define FhSystemEnv_H

#include "src/FhStringMap.h"
#include "commonRootData/FileHeader.h"
#include <GaudiKernel/ISvcLocator.h>

/*!

 @class FhSystemEnv   
 @brief Prepare/interpret system variables for/from a file header

 @author David Chamont - CNRS IN2P3 LLR Ecole Polytechnique

*/

class FhSystemEnv {	

public:

    FhSystemEnv() ;
    
    // global read/write methods
    void init( const FileHeader *  ) ;
    int init() ;
    void store( FileHeader * ) const ;
   
    // David: I keep an interface dependant on ROOT classes,
    // rather than std ones, because I want to
    // use TRegexp.
   
    // get values
    void getVariableNames( TRegexp exp, TCollection & names ) const ;
    TString getValue( const TObjString & name ) const ;
    
private:

    void rawToMap() ;
	TString m_raw ;
	FhStringMap m_map ;
    
} ;

#endif


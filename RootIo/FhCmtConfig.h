
#ifndef FhCmtConfig_H
#define FhCmtConfig_H

#include "src/FhStringMap.h"
#include "commonRootData/FileHeader.h"
#include <GaudiKernel/ISvcLocator.h>

/*!

 @class FhCmtConfig   
 @brief Prepare/interpret cmt configuration for/from a file header

 @author David Chamont - CNRS IN2P3 LLR Ecole Polytechnique

*/

class FhCmtConfig {	

public:

    FhCmtConfig() ;
    
    // global read/write methods
    void init( const FileHeader *  ) ;
    void init() ;
    void store( FileHeader * ) const ;
   
    // David: I keep an interface dependant on ROOT classes,
    // rather than std ones, because I want to
    // use TRegexp.
   
    // get values
    void getPackageNames( TRegexp exp, TCollection & names ) const ;
    TString getPackageSelection( const TObjString & name ) const ;
    TString getUseHierarchy() const ;
    void getMacroNames( TRegexp exp, TCollection & names ) const ;
    TString getMacroValue( const TObjString & name ) const ;
    
private:

    void rawToMap() ;
	TString m_rawUses ;
	TString m_useHierarchy ;
	FhStringMap m_mapUses ;
	TString m_rawMacros ;
	FhStringMap m_mapMacros ;
    
} ;

#endif



#include "src/FhStringMap.h"
#include <iostream>

FhStringMap::FhStringMap() {
}

FhStringMap::~FhStringMap() {
	DeleteAll() ;
}

void FhStringMap::DeleteAll() {
	m_strings.DeleteAll() ;
}

void FhStringMap::add( const TObjString & name, const TString & value ) {
    m_strings.Add(name.Clone(),new TObjString(value)) ;
}

TString FhStringMap::getValue( const TObjString & name ) const {
	TObject * value = m_strings.GetValue(name.GetName()) ;
	if (!value)
	  return "Empty String" ;
	else
	  return ((TObjString *)value)->String() ;
}
 
void FhStringMap::getNames(  TRegexp exp, TCollection & result ) const {
    TIterator * iter = m_strings.MakeIterator() ;
    TObjString * key ;
    Ssiz_t len ;
    while (( key = (TObjString *)(iter->Next()) )) {
        exp.Index(key->String(),&len,0) ;
        if (len>0) {
            result.Add(key) ;
        }
    }
    delete iter ;
}
 

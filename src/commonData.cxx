#include "commonData.h"


std::map<const Event::McParticle*, TRef> commonData::m_mcPartMap;
std::map<const Event::McIntegratingHit*, TRef> commonData::m_mcIntHitMap;
std::map<const Event::McPositionHit*, TRef> commonData::m_mcPosHitMap;
std::map<const Event::TkrDigi*, TRef> commonData::m_tkrDigiMap;
std::map<const Event::CalDigi*, TRef> commonData::m_calDigiMap;

McEvent* commonData::m_mcEvt;
DigiEvent* commonData::m_digiEvt;


void commonData::clear() { 
    m_mcPartMap.clear(); 
    m_mcIntHitMap.clear();
    m_mcPosHitMap.clear();
    m_calDigiMap.clear(); 
}

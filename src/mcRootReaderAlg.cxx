#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "GlastEvent/TopLevel/Event.h"
#include "GlastEvent/TopLevel/MCEvent.h"
#include "GlastEvent/TopLevel/EventModel.h"
#include "GlastEvent/MonteCarlo/McParticle.h"
#include "GlastEvent/MonteCarlo/McIntegratingHit.h"
#include "GlastEvent/MonteCarlo/McPositionHit.h"


#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TObjArray.h"
#include "TCollection.h"  // Declares TIter

#include "mcRootData/McEvent.h"

#include "facilities/Util.h"

#include <map>

/** @class mcRootReaderAlg
 * @brief Reads Monte Carlo data from a persistent ROOT file and stores the
 * the data in the TDS.
 *
 * @author Heather Kelly
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/mcRootReaderAlg.cxx,v 1.2 2002/04/29 14:17:31 heather Exp $
 */

class mcRootReaderAlg : public Algorithm
{	
public:
    
    mcRootReaderAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    /// Handles setup by opening ROOT file in read mode and creating a new TTree
    StatusCode initialize();
   
    /// Orchastrates reading from ROOT file and storing the data on the TDS for each event
    StatusCode execute();
    
    /// Closes the ROOT file and cleans up
    StatusCode finalize();
            
        
private:

    /// Retrieves McParticles from the ROOT file and stores them on the TDS
    StatusCode readMcParticles();

    /// Retrieves McPositionHits from the ROOT file and stores them on the TDS
    StatusCode readMcPositionHits();

    /// Retrieves McIntegratingHits from the ROOT file and stores them on the TDS
    StatusCode readMcIntegratingHits();

    /// Converts from ROOT's VolumeIdentifier to idents::VolumeIdentifier 
    void convertVolumeId(VolumeIdentifier rootVolId, idents::VolumeIdentifier &tdsVolId);

    /// Closes the ROOT file
    void close();
   
    /// ROOT file pointer
    TFile *m_mcFile;
    /// ROOT tree pointer
    TTree *m_mcTree;
    /// Top-level Monte Carlo ROOT object
    McEvent *m_mcEvt;
    /// name of the output ROOT file
    std::string m_fileName;
    /// name of the Monte Carlo TTree stored in the ROOT file
    std::string m_treeName;

    /// Keep track of MC particles as we retrieve them from the ROOT file
    /// The id is the unique id assigned by ROOT for every TObject.
    std::map<UInt_t, mc::McParticle*> m_particleMap;
    
};

static const AlgFactory<mcRootReaderAlg>  Factory;
const IAlgFactory& mcRootReaderAlgFactory = Factory;


mcRootReaderAlg::mcRootReaderAlg(const std::string& name, ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator)
{
    // Input pararmeters that may be set via the jobOptions file
    declareProperty("mcRootFile",m_fileName="mc.root");
    declareProperty("mcTreeName", m_treeName="Mc");

    m_particleMap.clear();
}

StatusCode mcRootReaderAlg::initialize()
{
    // Purpose and Method:  Called once before the run begins.  This method
    //    opens a new ROOT file and prepares for reading.

    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    // Use the Job options service to set the Algorithm's parameters
    // This will retrieve parameters set in the job options file
    setProperties();

    facilities::Util::expandEnvVar(&m_fileName);
    
    // Save the current directory for the ntuple writer service
    TDirectory *saveDir = gDirectory;   
    m_mcFile = new TFile(m_fileName.c_str(), "READ");
    if (!m_mcFile->IsOpen()) sc = StatusCode::FAILURE;
    m_mcFile->cd();
    m_mcTree = (TTree*)m_mcFile->Get(m_treeName.c_str());
    m_mcEvt = 0;//new McEvent();
    m_mcTree->SetBranchAddress("McEvent", &m_mcEvt);
    
    saveDir->cd();
    return sc;
    
}

StatusCode mcRootReaderAlg::execute()
{
    // Purpose and Method:  Called once per event.  This method calls
    //   the appropriate methods to read data from the ROOT file and store
    //   data on the TDS.

    MsgStream log(msgSvc(), name());

    StatusCode sc = StatusCode::SUCCESS;

    static UInt_t evtId = 0;
    m_mcTree->GetEvent(evtId);

    m_particleMap.clear();

    MCEvent* mcEventObj = 
        SmartDataPtr<MCEvent>(eventSvc(), EventModel::MC::Event);
    if (!mcEventObj) {
        log << "Failed to retrieve McEvent" << endreq;
        return sc;
    }

    sc = readMcParticles();
    if (sc.isFailure()) return sc;

    sc = readMcPositionHits();
    if (sc.isFailure()) return sc;

    sc = readMcIntegratingHits();
    if (sc.isFailure()) return sc;
   
    m_mcEvt->Clear();
    evtId++;
    
    return sc;
}


StatusCode mcRootReaderAlg::readMcParticles() {
    // Purpose and Method:  Retrieve the McParticle collection from the ROOT file 
    //    and fill the TDS McParticle collection
    // Input:  ROOT McParticle collection
    // TDS Output:  EventModel::MC::McParticleCol

    MsgStream log(msgSvc(), name());

    StatusCode sc = StatusCode::SUCCESS;

    TObjArray *particles = m_mcEvt->getMcParticleCol();
    if (!particles) return sc;
    TIter partIter(particles);

    // create the TDS location for the McParticle Collection
    mc::McParticleCol* pTdsCol = new mc::McParticleCol;
    sc = eventSvc()->registerObject(EventModel::MC::McParticleCol, pTdsCol);
    if (sc.isFailure()) {
        log << "Failed to register McParticle Collection" << endreq;
        return sc;
    }

    McParticle *pRoot;
    // Create the map of ROOT unique ids and TDS McParticle objects
    while (pRoot = (McParticle*)partIter.Next()) {
        mc::McParticle *pTds = new mc::McParticle();
        m_particleMap[pRoot->GetUniqueID()] = pTds;
    }

    // reset the iterator to the beginning of the list of McParticles
    partIter.Reset();

    // Now that the map is available, we initialize all of the TDS McParticles
    while (pRoot = (McParticle*)partIter.Next()) {

        mc::McParticle *pTds = m_particleMap[pRoot->GetUniqueID()];

        int idTds = pRoot->getParticleId();

        unsigned int statusBitsTds = pRoot->getStatusFlags();
        
        TLorentzVector initialMomRoot = pRoot->getInitialFourMomentum();
        
        HepLorentzVector initialMomTds(initialMomRoot.X(), initialMomRoot.Y(),
            initialMomRoot.Z(), initialMomRoot.T());

        TLorentzVector finalMomRoot = pRoot->getFinalFourMomentum();
        HepLorentzVector finalMomTds(finalMomRoot.X(), finalMomRoot.Y(), 
            finalMomRoot.Z(), finalMomRoot.T());

        TVector3 finalPosRoot = pRoot->getFinalPosition();
        HepPoint3D finalPosTds(finalPosRoot.X(), finalPosRoot.Y(), finalPosRoot.Z());

        const McParticle *momRoot = pRoot->getMother();

        mc::McParticle *momTds = 0;

        if (momRoot != 0) {
            // The case where this is the primary particle
            if (momRoot->GetUniqueID() == pRoot->GetUniqueID()) {
                momTds = pTds;
            } else if (m_particleMap.find(momRoot->GetUniqueID()) != m_particleMap.end()) {
                momTds = m_particleMap[momRoot->GetUniqueID()];
            } else {
                log << MSG::INFO << "Failed to find the McParticle in the map!!" << endreq;
            }

        }

        // Setup the TDS version fo the McParticle
        pTds->init(momTds, idTds, statusBitsTds, initialMomTds, finalMomTds, finalPosTds);

        // Add the TDS McParticle to the TDS collection of McParticles
        pTdsCol->push_back(pTds);

    }

    return sc;
}

StatusCode mcRootReaderAlg::readMcPositionHits() {
    // Purpose and Method:  Retrieve the McPositionHit collection from the ROOT
    //  file and fill the TDS McPositionHit collection.
    // Input:  ROOT McPositionHit collectoin
    // TDS Output:  EventModel::MC::McPositionHitCol

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    //SmartDataPtr<McPositionHitVector> posHits(eventSvc(), EventModel::MC::McPositionHitCol);    
    
    const TObjArray *posHits = m_mcEvt->getMcPositionHitCol();
    if (!posHits) return sc;
    TIter hitIter(posHits);

    // create the TDS location for the McParticle Collection
    McPositionHitVector* pTdsCol = new McPositionHitVector;
    sc = eventSvc()->registerObject(EventModel::MC::McPositionHitCol, pTdsCol);
    if (sc.isFailure()) {
        log << "Failed to register McPositionHit Collection" << endreq;
        return sc;
    }

    const McPositionHit *posHitRoot;
    while (posHitRoot = (McPositionHit*)hitIter.Next()) {
    
        mc::McPositionHit *posHitTds = new mc::McPositionHit();

        VolumeIdentifier volIdRoot = posHitRoot->getVolumeId();
        idents::VolumeIdentifier volIdTds;
        convertVolumeId(volIdRoot, volIdTds);

        TVector3 entryRoot = posHitRoot->getEntryPosition();
        HepPoint3D entryTds(entryRoot.X(), entryRoot.Y(), entryRoot.Z());

        TVector3 exitRoot = posHitRoot->getExitPosition();
        HepPoint3D exitTds(exitRoot.X(), exitRoot.Y(), exitRoot.Z());
        
        double edepTds= posHitRoot->getDepositedEnergy();
        
        double epartTds = posHitRoot->getParticleEnergy();
        posHitTds->setParticleEnergy(epartTds);
        
        double tofTds = posHitRoot->getTimeOfFlight();
        posHitTds->setTimeOfFlight(tofTds);
        
        const McParticle *mcRoot = posHitRoot->getMcParticle();
        mc::McParticle *mcTds = 0;
        if (mcRoot != 0) {
            mcTds = m_particleMap[mcRoot->GetUniqueID()];
            posHitTds->setMcParticle(mcTds);
        }
      
        const McParticle *originRoot = posHitRoot->getOriginMcParticle();
        mc::McParticle *originTds = 0;
        if (originRoot != 0) {
            originTds = m_particleMap[originRoot->GetUniqueID()];
            posHitTds->setOriginMcParticle(originTds);
        }

        // setup the TDS McPositionHit
        posHitTds->init(edepTds, volIdTds, entryTds, exitTds);
        // add the McPositionHit to the TDS collection of McPositionHits
        pTdsCol->push_back(posHitTds);
    }
    
    return sc;
}

StatusCode mcRootReaderAlg::readMcIntegratingHits() {
    // Purpose and Method:  Retrieve the McIntegratingHit collection from the
    //     ROOT file and fill the TDS McIntegratingHit collection.
    // Input:  ROOT McIntegratingHit collection
    // TDS Output:  EventModel::MC::McIntegratingHitCol

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    const TObjArray *intHits = m_mcEvt->getMcIntegratingHitCol();
    if (!intHits) return sc;
    TIter hitIter(intHits);

    // create the TDS location for the McParticle Collection
    McIntegratingHitVector* pTdsCol = new McIntegratingHitVector;
    sc = eventSvc()->registerObject(EventModel::MC::McIntegratingHitCol, pTdsCol);
    if (sc.isFailure()) {
        log << "Failed to register McIntegratingHit" << endreq;
        return sc;
    }

    McIntegratingHit *intHitRoot;
    while (intHitRoot = (McIntegratingHit*)hitIter.Next()) {

        mc::McIntegratingHit *intHitTds = new mc::McIntegratingHit();

        const VolumeIdentifier idRoot = intHitRoot->getVolumeId();
        idents::VolumeIdentifier idTds;
        convertVolumeId(idRoot, idTds);

        intHitTds->setVolumeID(idTds);

        const McIntegratingHit::energyDepositMap mcPartMapRoot = intHitRoot->getItemizedEnergy();
        McIntegratingHit::energyDepositMap::const_iterator rootMapIt;
        log << MSG::DEBUG << "EnergyMap size: " << mcPartMapRoot.size() << endreq;
        // Can't seem to read this back in due to TRef problem in Root 3.02.03
        /*
        for (rootMapIt = mcPartMapRoot.begin(); rootMapIt != mcPartMapRoot.end(); rootMapIt++){
            McParticle* mcPartRoot = rootMapIt->first;
            mc::McParticle *mcPartTds = m_particleMap[mcPartRoot->GetUniqueID()];
            double e = rootMapIt->second;
            TVector3 posRoot = mcPartRoot->getFinalPosition();
            HepPoint3D posTds(posRoot.X(), posRoot.Y(), posRoot.Z());
            intHitTds->addEnergyItem(e, mcPartTds, posTds);
        }
        */

        // Add the TDS McIntegratingHit to the TDS McIntegratingHit collection
        pTdsCol->push_back(intHitTds);
    }

    return sc;
}

void mcRootReaderAlg::convertVolumeId(VolumeIdentifier rootVolId, 
                                      idents::VolumeIdentifier& tdsVolId) 
{
    // Purpose and Method:  We must store the volume ids as two 32 bit UInt_t
    //     in the ROOT class.  The idents::VolumeIdentifier class stores the
    //     data in one 64 bit word.  We must convert from the two 32 bit words
    //     into the 64 bit word.  We perform the conversion by iterating over
    //     all of the ids in the ROOT VolumeIdentifier and appending them to
    //     the TDS idents::VolumeIdentifier.
    // Input:  ROOT VolumeIdentifier
    // Ouput:  idents::VolumeIdentifier

    unsigned int index;
    for (index = 0; index < rootVolId.size(); index++) {
        tdsVolId.append(rootVolId.operator [](index));
    }

}


void mcRootReaderAlg::close() 
{
    // Purpose and Method:  Writes the ROOT file at the end of the run.
    //    The TObject::kOverWrite parameter is used in the Write method
    //    since ROOT will periodically write to the ROOT file when the bufSize
    //    is filled.  Writing would create 2 copies of the same tree to be
    //    stored in the ROOT file, if we did not specify kOverwrite.

    TDirectory *saveDir = gDirectory;
    m_mcFile->cd();
    m_mcFile->Close();
    saveDir->cd();
}

StatusCode mcRootReaderAlg::finalize()
{
    close();
    
    StatusCode sc = StatusCode::SUCCESS;
    return sc;
}
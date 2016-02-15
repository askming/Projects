/*========================================================================================
 *
 *  qreg.h - definitions/declarations for QREG class
 *
 *========================================================================================*/

#ifndef _INCLUDED_QREG_
#define _INCLUDED_QREG_


enum {GIBBS, RW};
enum {NOSCALE, SCALE}


class QREG
{
    // data
    decl m_vY;
    decl m_mX;
    decl m_vXB;
    
    // dimensions
    decl m_iN;
    decl m_iK;

    // samplers
    decl m_iMethod;
    decl m_iScale;

    // parameters
    decl m_dProb;
    decl m_dTheta;
    decl m_dTau2;

    decl m_vBeta;
    decl m_dSigma;
    decl m_vV;

    // prior
    decl m_vBeta0;
    decl m_mInvB0;
    decl m_dN0;
    decl m_dS0;

    // MCMC
    decl m_iBurn;
    decl m_iDraw;
    decl m_mPostBeta;
    decl m_mPostSigma;
    decl m_dAccept;
    decl m_dScale;


    // member functions
    QREG( const dP );
    LoadData( const sFile );
    InfoData();
    SetPrior( const vBeta0, const vSigma0 );
    SetMCMC( const iMethod, const iScale );
    RunMCMC( const vOpt );
    PostSum();

    traceMCMC( const iter, const nstep );
    genBeta();
    genSigma();
    genV();
    
    genBeta_MH();
    genSigma_MH();
    logPostBeta( const vBeta );
}


#endif

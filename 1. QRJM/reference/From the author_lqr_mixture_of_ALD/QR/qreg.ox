/*========================================================================================
 *
 *  qreg.ox - QREG class
 *
 *========================================================================================*/
#include <oxstd.h>
#include <oxprob.h>
#include <oxfloat.h>
#include <oxdraw.h>
#import  <maximize>
#include "qreg.h"
#include "subroutine.ox"


const decl EPS = 1.0e-8;


QREG::QREG( const dP )
{
    // set parameters
	m_dProb  = dP;
    m_dTheta = (1.0 - 2.0*m_dProb)/(m_dProb - m_dProb*m_dProb);
    m_dTau2  = 2.0/(m_dProb - m_dProb*m_dProb);

    // set format
    format( 1024 );
    format( "%#10.3f" );

    // print
    println( DLINE, DLINE );
    println( "Bayesian Analysis of Quantile Regression" );
    println( "%80s", time() ~ "/" ~ date() );
}


QREG::LoadData( const sFile )
{
    decl data;

    data = loadmat( sFile );
    m_vY = data[][0];
    m_mX = 1 ~ data[][1:];

    // dimensions
    m_iN = rows( m_vY );
    m_iK = columns( m_mX );
}


QREG::InfoData()
{
    decl i, j, ydat, xdat, nc = 7;

    title( "DATA SUMMARY" );

    // y
    ydat = statistics( m_vY );

    // x
    xdat = zeros( m_iK, nc );
    for (i = 0; i < m_iK; i++) {
        xdat[i][] = statistics( m_mX[][i] );
    }

    // print
    println( "number of observations =", "%10d",   m_iN );
    println( "number of regressors   =", "%10d",   m_iK );
    println( "probability            =", "%10.3f", m_dProb );
    print( "\n", "%10s", "" );
    for (i = 0; i < nc; i++) {
        print( "%10s", CSTAT[i] );
    }
    println( "\n", LINE, LINE );
    print( "%-10s", "y" );
    for (i = 0; i < nc; i++) {
        print( ydat[i] );
    }
    print( "\n" );
    for (i = 0; i < m_iK; i++) {
        print( "%-10s", sprint( "x", i ) );
        for (j = 0; j < nc; j++) {
            print( xdat[i][j] );
        }
        print( "\n" );
    }
    println( LINE, LINE );
    println( NL );
}


QREG::SetPrior( const vBeta0, const vSigma0 )
{
    title( "PRIOR" );

    // beta
    m_vBeta0 = ones( m_iK, 1 ) * vBeta0[0];
    m_mInvB0 = unit( m_iK ) / vBeta0[1];

    // sigma
    m_dN0 = vSigma0[0];
    m_dS0 = vSigma0[1];

    // print
    println( "beta  ~  N(", vBeta0[0], ",", vBeta0[1], " )" );
    println( "sigma ~ IG(", m_dN0, "/2,", m_dS0, "/2 )" );
    println( NL );
}


QREG::SetMCMC( const iMethod, const iScale )
{
    decl i;

    // set scale and the algorithm
    m_iMethod = iMethod;
    m_iScale  = iScale;

    // initial values
    m_vBeta  = ones( m_iK, 1 );
    m_vXB    = m_mX * m_vBeta;
    m_dSigma = 1.0;
    if (m_iMethod == GIBBS) {
       m_vV = 0.1 * ones( m_iN, 1 );
    }
}


QREG::RunMCMC( const vOpt )
{
    decl it, nn, nstep;

	title( "MCMC" );

    // set iterations
    m_iBurn = vOpt[0];
    m_iDraw = vOpt[1];
    nn      = m_iBurn + m_iDraw;
    nstep   = int( 0.1 * nn );

    // set scale parameter for RW
    if (m_iMethod == RW) m_dScale = vOpt[2];

    // set parameters
    m_mPostBeta  = zeros( m_iK, nn );
    m_mPostSigma = zeros( 1,    nn );
    m_dAccept = 0;

    // print
    println( "burn-in period =", "%10d", m_iBurn );
    println( "no. of draws   =", "%10d", m_iDraw );
    println( "scale          =", "%10s", m_iScale == SCALE ? "YES" : "NO", "\n" );

    // main loop
    for (it = 0; it < nn; it++) {
        genV();
        genBeta();
        if (m_iScale == SCALE) {
            genSigma();
        }
        traceMCMC( it, nstep );
        
        m_mPostBeta[][it]  = m_vBeta;
        m_mPostSigma[][it] = m_dSigma;
    }
    println( LINE, LINE );
    println( "acceptance rate = ", "%10.2f", 100*m_dAccept/nn, " %" );  
    println( NL );
}


QREG::traceMCMC( const iter, const nstep )
{
    if (iter==0) {
        println( 
            "%8s", "iter", 
            "%10s", "beta[0]", 
            "%10s", "beta[1]", 
            "%10s", "beta[2]", 
            "%10s", "sigma" );
        println( LINE, LINE );
    }
    if ( !fmod(iter, nstep) || iter==-m_iBurn ) {
        println( "%8d", iter, "%10.3f", m_vBeta[0], "%10.3f", m_vBeta[1], "%10.3f", m_vBeta[2],
        "%10.3f", m_dSigma );
    }
}


QREG::PostSum()
{
    decl pmntBeta, pmntSigma;
    decl i, j, bm, asCol;

    title( "POSTERIOR SUMMARY" );

    bm = m_iDraw * 0.1;
    asCol = {"mean", "std", "95% CI", "IF"};

    // posterior summary
    pmntBeta  = postsummary( m_mPostBeta[][m_iBurn:],  bm );
    pmntSigma = postsummary( m_mPostSigma[][m_iBurn:], bm );

    // print
    print( "%10s", "" );
	for (i = 0; i < 4; i++) {
		if (i == 2) {
			print( "%17s", asCol[i], "%6s", "" );
		}
		else {
			print( "%10s", asCol[i] );
		}
	}
	println( "\n", LINE, LINE );

	// beta
	for (i = 0; i < m_iK; i++) {
		print( "%-10s", sprint( "beta[", i, "]" ) );
		for (j = 0; j < 2; j++) {
			print( pmntBeta[i][j] );
		}
		print( "%5s", "(", "%8.3f", pmntBeta[i][2], "," );
		print( "%8.3f", pmntBeta[i][3], ")" );
		println( pmntBeta[i][4] );
	}
	
	// sigma
	print( "%-10s", "sigma" );
	for (j = 0; j < 2; j++) {
		print( pmntSigma[j] );
	}
	print( "%5s", "(", "%8.3f", pmntSigma[2], "," );
	print( "%8.3f", pmntSigma[3], ")" );
	println( pmntSigma[4] );
	println( LINE, LINE );
    println( NL );
}


QREG::genBeta()
{
    decl vY, mX, bm, B;

    // set X
    mX = m_mX ./ (m_dSigma * m_dTau2 * m_vV);

    // mean and covariance
    B  = invertsym( mX'm_mX + m_mInvB0 );
    bm = B * ( mX'(m_vY - m_dTheta*m_vV)  + m_mInvB0*m_vBeta0 );

    // update beta
    m_vBeta = bm + choleski( B ) * rann( m_iK, 1 );
    m_vXB   = m_mX * m_vBeta;
}


QREG::genSigma()
{
	decl n1, e, e2, s1;

	e  = (m_vY - m_vXB - m_dTheta*m_vV) ./ sqrt(m_vV);
	e2 = sumsqrc( e ) / m_dTau2;

	// set parameters
	n1 = 0.5*(m_iN + m_dN0) + m_iN;
	s1 = 0.5*(e2   + m_dS0) + sumc( m_vV );
	
	// update sigma
	m_dSigma = 1.0 / rangamma( 1, 1, n1, s1 );
}


QREG::genV()
{
    decl i, nu, aa, bb, st2, v, u;

    // set parameter
    st2 = m_dSigma * m_dTau2;
    nu  = 0.5;
    bb  = sqrt( m_dTheta^2/st2 + 2.0/m_dSigma );
    u   = m_vY - m_vXB;

    // update v
    for (i = 0; i < m_iN; i++) {
        aa = sqrt( u[i]^2 / st2 );
        v = rangig( 1, 1, nu, aa, bb );
        if ( isnan( v ) || v < EPS ) {
            continue;
        }
        m_vV[i] = v;
    }
}


QREG::genBeta_MH()
{
    decl beta_o, beta_n, logpost_o, logpost_n;
    decl a, u;

    // old
    beta_o    = m_vBeta;
    logpost_o = logPostBeta( beta_o );

    // candidate
    beta_n    = beta_o + m_dScale*rann( m_iK, 1 );
    logpost_n = logPostBeta( beta_n );

    a = exp( logpost_n - logpost_o );
    u = ranu( 1, 1 );
    if ( u < a ) {
        m_vBeta = beta_n;
        ++m_dAccept;
    }
}


QREG::logPostBeta( const vBeta )
{
    decl u, rho, lnf, lnp;

    u   = m_vY - m_mX*vBeta;
    rho = 0.5*(fabs(u) + (2.0*m_dProb - 1.0)*u) / m_dSigma;
    lnf = -sumc( rho );
    lnp = -0.5*(vBeta - m_vBeta0)'*m_mInvB0*(vBeta - m_vBeta0);

    return lnf + lnp;
}


QREG::genSigma_MH()
{
    decl u, sse, n1, s1;
    
    u   = m_vY - m_vXB;	
    sse = sumc( (fabs(u) + (2.0*m_dProb - 1.0)*u) ); 

    // set parameters
    n1 = m_iN + 0.5*m_dN0;
    s1 = 0.5*(sse  + m_dS0);
    
    // update sigma
    m_dSigma = 1.0 / rangamma( 1, 1, n1, s1 );
}


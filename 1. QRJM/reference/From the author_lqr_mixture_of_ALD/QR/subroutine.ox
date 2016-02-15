/*========================================================================================
 *
 *  subroutine.ox - Ox file for subrouines
 *
 *========================================================================================*/
#ifndef _SUBROUTINE_INCLUDED_
#define _SUBROUTINE_INCLUDED_

#include <oxprob.h>
#include <oxfloat.h>


// constant variables
const decl SEED  = 1234;
const decl DLINE = "========================================";
const decl LINE  = "----------------------------------------";
const decl DOT   = "........................................";
const decl NL    = "\n";
const decl CSTAT = {"mean", "std", "min", "5%", "50%", "95%", "max"};


static title( const sT )
{
    println( DLINE, DLINE );
    println( sT );
    println( DLINE, DLINE );
}



// basic statistics
static statistics( const vX )
{
	decl nc   = 7;
	decl q    = <0.025, 0.5, 0.975>;
    decl xdat = zeros( 1, nc );

    xdat[0]   = meanc( vX );
    xdat[1]   = sqrt( varc( vX ) );
    xdat[2]   = min( vX );
    xdat[3:5] = quantilec( vX, q );
    xdat[6]   = max( vX );

    return xdat;
}


static fTvar( const mX, const dBm )
{
	decl msp = periodogram( mX', dBm, 2, 1 ) / columns(mX);

	return M_2PI*( msp[0][] )';
}


// function for posterior summary
static postsummary( const mX, const dBm )
{
    decl npar = rows( mX );
    decl xdat = zeros( npar, 5 );

    xdat[][0] = meanr( mX );
    xdat[][1] = sqrt( varr( mX ) );
    xdat[][2] = quantiler( mX, 0.025 );
    xdat[][3] = quantiler( mX, 0.975 );
    xdat[][4] = fTvar( mX, dBm ) ./ varr( mX );

    return xdat;
}


// truncated normal s.t. x > L
static Rantnl( const vMu, const vS2, const vLeft )
{
    decl s  = sqrt( vS2 );
    decl u  = ranu( sizer( vMu ), 1 );
    decl pr = 1.0;
    decl pl = probn( (vLeft  - vMu) ./ s );
    decl x  = vMu + s .* quann( u .* (pr - pl) + pl );

    return x;
}



// truncated normal s.t. x < R
static Rantnr( const vMu, const vS2, const vRight )
{
    decl s  = sqrt( vS2 );
    decl u  = ranu( sizer( vMu ), 1 );
    decl pr = probn( (vRight - vMu) ./ s );
    decl pl = 0.0;
    decl x  = vMu + s .* quann( u .* (pr - pl) + pl );

    return x;
}



// truncated normal s.t. L < x < R
static Rantnlr( const vMu, const vS2, const vLeft, const vRight )
{
    decl s  = sqrt( vS2 );
    decl u  = ranu( sizer( vMu ), 1 );
    decl pl = probn( (vLeft  - vMu) ./ s );
    decl pr = probn( (vRight - vMu) ./ s );
    decl x  = vMu + s .* quann( u .* (pr - pl) + pl );

    return x;
}



// multivariate t distribution
static Ranmvt( const vMu, const mV, const dDF )
{
    decl n = rows( vMu );
    decl s = rangamma( 1, 1, 0.5*dDF, 0.5*dDF );
    decl x = vMu + choleski( mV/s ) * rann( n, 1 );

    return x;
}


// Wishart distribution
static Ranwishart( const mS, const dDF )
{
    decl n = sizer( mS );
    decl c = choleski( mS );
    decl x = zeros( n, n );
    decl h, i;

    for (i = 0; i < n; i++) {
        x[i][i] = sqrt( ranchi(1, 1, dDF-i) );
        if (i < n-1) x[i][i+1:] = rann(1, n-i-1);
    }
    h = x*c';

    return h'h;
}


// log density of multivarite t distribution
static logdensMvt( const vX, const vMu, const mV, const dDF )
{
    decl p   = rows( vX );
    decl det = determinant( mV );
    decl sse = (vX - vMu)'invert( mV )*(vX - vMu);
    decl lnf = - 0.5*log( det ) - 0.5*(dDF + p)*log( 1.0 + sse/dDF );

    return double( lnf );
}


// log density of t distribution
static logdenst( const x, const mu, const s2, const df )
{
    decl tmp = 1 + (x - mu)^2 / (df*s2);
    decl f   = -0.5*(df + 1)*log( tmp );

    return double(f);
}


// discrete distribution
static Ranprob( const vProb )
{
    decl u, s;

    // check vProb
    if ( any( vProb .< 0.0 ) )
        eprint( "ERROR! vProb must be non-negative." );

    u = ranu(1,1);
    s = cumulate( vProb/sumc( vProb ) );

    return( sumc( s .<= u ) );
}

#endif
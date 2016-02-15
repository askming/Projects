/*========================================================================================
 *
 *  xqreg.ox
 *
 *========================================================================================*/
#include <oxstd.h>
#include <oxprob.h>
#include <oxfloat.h>
#include <oxdraw.h>
#import  <maximize>
#include "qreg.ox"
#include "subroutine.ox"

main()
{
	decl beta0  = <0, 100>;        // hyper-parameters for beta
	decl sigma0 = <3, 0.1>;        // hyper-parameters for sigma^2
	decl opt    = <1000, 5000>;    // burn-in period, # of draws for posterior inference
	decl obj, time;

    time = timer();
    obj =  new QREG( 0.5 );
    obj -> LoadData( "test.csv" );
    obj -> InfoData();
    obj -> SetPrior( beta0, sigma0 );
    obj -> SetMCMC( GIBBS, SCALE ); // choice of methods, whether sigma^2 is included or not
    obj -> RunMCMC( opt );
    obj -> PostSum();
    delete obj;

    // time
    println( DLINE, DLINE );
    println( "Execution time:   ", "%62s", timespan( time ) );
    println( "Program finished: ", "%62s", date() );
    println( DLINE, DLINE );
}


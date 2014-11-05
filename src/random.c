//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Filename:  random.c
// 
// Purpose:   Random utility procedures for BayeSys3.
// 
// History:   Random.c  17 Nov 1994 - 13 Sep 2003
//
// Acknowledgments:
//   "Numerical Recipes", Press et al, for ideas
//   "Handbook of Mathematical Functions", Abramowitz and Stegun, for formulas
//    Peter J Acklam website, for inverse normal code
//-----------------------------------------------------------------------------
/*
    Copyright (c) 1994-2003 Maximum Entropy Data Consultants Ltd,
                            114c Milton Road, Cambridge CB4 1XE, England

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include "license.txt"
*/
#include <stddef.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include "random.h"

static const double SHIFT32   = 1.0 / 4294967296.0;          // 2^-32
static const double HalfPi    = 1.57079632679489661922;      // pi/2
static const double TwoxPi    = 6.28318530717958647688;      // 2*pi
static const double SqrPi2    = 1.25331413731550025120;      // sqrt(pi/2)
static const double Sqr2Pi    = 2.50662827463100050240;      // sqrt(2*pi)
static const double LogSqr2Pi = 0.91893853320467274177;      // log(sqrt(2*pi))

// Internal prototypes
static void   Positive2   (Rand_t, double, double, double, double*, double*);
static void   Positive2A  (double, double, double, double,
                           double*, double*, double*);
static void   Positive2B  (double, double, double, double,
                           double*, double*, double*);
static void   Positive2C  (double, double, double, double,
                           double*, double*, double*);
static void   Positive2D  (double, double, double*, double*, double*);
static void   Posneg2     (Rand_t, double, double, double, double, double,
                           double*, double*);
static void   Posneg2A    (double, double, double, double, double, double,
                           double*, double*, double*);
static void   Posneg2B    (double, double, double, double, double, double,
                           double*, double*, double*);
static void   Posneg2C    (double, double, double, double, double, double,
                           double*, double*, double*);
static void   Posneg2D    (double, double, double, double, double, double,
                           double*, double*, double*);
static double eigenerf2   (double, double, double, double, double);
static double erfint2     (double, double);
static double wedge       (double, double, double);

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  RanInit
//
// Purpose:   Initialise generator of 4 unsigneds for other Random procedures
//
// Format:    Rand[0] + Rand[1] * 2^32 = # calls to Ranint since initialisation
//            Rand[2] = generator offset, re-randomised every 2^32 calls
//            Rand[3] = extra random integer available after each call
//
// History:   John Skilling   15 Jan 2002, 31 Oct 2002
//-----------------------------------------------------------------------------
// 
int RanInit(          //   O  Seed, either from input or time
Rand_t  Rand,         //   O  Random generator state          [4]
int     seed)         // I    Seed: +ve = value, -ve = time seed
{
    unsigned  j, k;

    k = 1;
    for( j = 0; k; ++j )
        k += k;
    if( j != 32 )                      // Check 32-bit arithmetic
        return E_RAN_ARITH;
    if( seed < 0 )
    {
        seed = (int)time(NULL);
        if( seed < 0 )
            seed = ~seed;              // still OK after A.D.2030
    }
    Rand[0] = Rand[1] = 0;                     // 64-bit counter
    Rand[2] = 1013904223 + 1664525 * seed;     // sticky offset
    Rand[3] = 1013904223 + 1664525 * Rand[2];  // extra random integer
    return  seed;
}
#if 0
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Kernel of alternative initialiser #1, crude congruential.
//-----------------------------------------------------------------------------
void  RanInit1(
Rand_t  Rand,         //   O  Random generator state         [1]
int     seed)         // I    Non-negative seed
{
    Rand[0] = seed;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Kernel of alternative initialiser #2, double congruential with shuffling.
//-----------------------------------------------------------------------------
void  RanInit2(
Rand_t  Rand,         //   O  Random generator state         [35]
int     seed)         // I    Non-negative seed
{
    int  i, j;
    i = (seed > 2147483397) ? seed - 2147483397 : seed + 1; // [1...2147483398]
    Rand[33] = i;
    for( j = 39; j >= 0; j-- )
    {
        i = 40014 * i - 2147483563 * (i / 53668);
        if( i < 0 )
            i += 2147483563;                                // [1...2147483562]
        if( j < 32 )
            Rand[j] = i;
    }
    Rand[34] = Rand[0];
    Rand[32] = i;
} 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Kernel of alternative initialiser #3, Knuth's subtractive.
//-----------------------------------------------------------------------------
void  RanInit3(
Rand_t  Rand,         //   O  Random generator state         [112]
int     seed)         // I    Non-negative seed
{
    unsigned  i, j, k, m;
    m = (unsigned)(161803398 - seed);
    Rand[54] = m;
    Rand[109] = 0;
    k = 1;
    for( i = 1; i < 55; ++i )
    {
        j = (21 * i - 1) % 55;
        Rand[j] = k;
        k = m - k;
        m = Rand[j];
        Rand[54+i] = i;
    }
    for( k = 0; k < 4; ++k )
        for( i = 0; i < 55; ++i )
            Rand[i] -= Rand[(i + 31) % 55];
    Rand[110] = 0;
    Rand[111] = 31;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Kernel of alternative initialiser #4, 64-bit hashing.
//-----------------------------------------------------------------------------
void  RanInit4(
Rand_t  Rand,         //   O  Random generator state         [3]
int     seed)         // I    Non-negative seed
{
    Rand[0] = Rand[1] = 0;
    Rand[2] = seed;
}
#endif

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Ranint
//
// Purpose:   Random integer sample, in [-2^31,2^31).
//
// Method:    64-bit hashing is on #calls + offset to produce 32-bit candidate
//            with algorithm from ran4 in "Numerical Recipes".
//            After each 2^32 calls, the 64-bit offset is incremented by
//            congruential generators on its 32-bit halves.
//
// Notes: (1) Period is just less than 2^95, extremely large.
//        (2) All random calls are directed through this procedure.
//        (3) After a call, Rand[3] is available as another random integer.
//
// History:   John Skilling   28 Jan 2002, 31 Oct 2002, 17 Dec 2002
//-----------------------------------------------------------------------------
// 
int Ranint(           //   O  32-bit integer
Rand_t  Rand)         // I O  Random generator state         [4]
{
    unsigned  m, n;   // 64-bit register
    unsigned  u, v, w;
    int       i;
// 64-bit counter, for hashing
    if( !(Rand[0] ++) )
    {
        Rand[1] ++;
// occasional update offset in Rand[2]
        i = Rand[2] >> 1;                    // [0...2147483647]
        i -= i / 24683721;                   // [0...2147483561]
        i++;                                 // [1...2147483562]
        i = 40014 * i - 2147483563 * (i / 53668);
        if( i < 0 )                          // (40014 * i) % 2147483563
            i += 2147483563;                 // [1...2147483562]
        i--;                                 // [0...2147483561]
        i += i / 24683720;                   // [0...2147483647] (with holes)
        Rand[2] = i << 1;                    // even
        if( 1013904223 + 1664525 * i < 0 )   // chance
            Rand[2]++;                       // even or odd
    }
// Two double steps of 64-bit hash
    n = Rand[0] + Rand[2];
    m = Rand[1] + Rand[2];
    w = n ^ 0xbaa96887;
    v = w >> 16;
    w &= 0xffff;
    u = (v - w) * (v + w);
    m ^= (((u >> 16) | (u << 16)) ^ 0xb4f0c4a7) + w * v;
    w = m ^ 0x1e17d32c;
    v = w >> 16;
    w &= 0xffff;
    u = (v - w) * (v + w);
    n ^= (((u >> 16) | (u << 16)) ^ 0x178b0f3c) + w * v;
    w = n ^ 0x03bcdc3c;
    v = w >> 16;
    w &= 0xffff;
    u = (v - w) * (v + w);
    m ^= (((u >> 16) | (u << 16)) ^ 0x96aa3a59) + w * v;
    w = m ^ 0x0f33d1b2;
    v = w >> 16;
    w &= 0xffff;
    u = (v - w) * (v + w);
    n ^= (((u >> 16) | (u << 16)) ^ 0xaa5835b9) + w * v;
    Rand[3] = m;
    return  n;
}
#if 0
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Kernel of alternative generator #1, crude congruential.
//-----------------------------------------------------------------------------
int Ranint1(          //   O  32-bit integer
Rand_t  Rand)         // I O  Random generator state         [1]
{
    Rand[0] = Rand[0] * 1664525 + 1013904223;
    return Rand[0];
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Kernel of alternative generator #2, double congruential with shuffling.
//-----------------------------------------------------------------------------
int Ranint2(      //   O  Odd integer with 85 holes: ONLY (2^31 - 85) VALUES
Rand_t  Rand)     // I O  Random generator state             [35]
{
    int i, j, k;
// First congruential sequence
    i = Rand[32];
    i = 40014 * i - 2147483563 * (i / 53668);
    if( i < 0 )
        i += 2147483563;
    Rand[32] = i;                      // [1...2147483562]
// Second congruential sequence
    k = Rand[33];
    k = 40692 * k - 2147483399 * (k / 52774);
    if( k < 0 )
        k += 2147483399;
    Rand[33] = k;                      // [1...2147483398]
// Subtract from shuffled table
    j = Rand[34] / 67108862;           // [0...31]
    k = Rand[j] - k;                   // [-2147483397...2147483561]
    if( k < 0 )
       k += 2147483563;                // (2^31 - 85) values
    Rand[34] = k;                      // [0...2147483562]
    Rand[j] = i;
// Expand result k
    k += k / 24970740;            // [0...2^31) with holes at 24970740*[1...85]
    return  (k << 1) | 1;
} 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Kernel of alternative generator #3, Knuth's subtractive.
//-----------------------------------------------------------------------------
int Ranint3(          //   O  32-bit integer
Rand_t  Rand)         // I O  Random generator state        [112]
{
    unsigned  j, k;
    j = Rand[110];    Rand[110] = Rand[55+j];
    k = Rand[111];    Rand[111] = Rand[55+k];
    return (Rand[j] -= Rand[k]);
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Kernel of alternative generator #4, 64-bit hashing.
//-----------------------------------------------------------------------------
int   Ranint4(        //   O  32-bit integer 
Rand_t  Rand)         // I O  Random generator state        [3]
{
    unsigned  u, v, w, m, n;
// 64-bit counter, for hashing
    if( !(Rand[0] ++) )
        Rand[1] ++;
// 64-bit hash
    n = Rand[0] + Rand[2];
    m = Rand[1] + Rand[2];
    w = n ^ 0xbaa96887;
    v = w >> 16;
    w &= 0xffff;
    u = (v - w) * (v + w);
    m ^= (((u >> 16) | (u << 16)) ^ 0xb4f0c4a7) + w * v;
    w = m ^ 0x1e17d32c;
    v = w >> 16;
    w &= 0xffff;
    u = (v - w) * (v + w);
    n ^= (((u >> 16) | (u << 16)) ^ 0x178b0f3c) + w * v;
    w = n ^ 0x03bcdc3c;
    v = w >> 16;
    w &= 0xffff;
    u = (v - w) * (v + w);
    m ^= (((u >> 16) | (u << 16)) ^ 0x96aa3a59) + w * v;
    w = m ^ 0x0f33d1b2;
    v = w >> 16;
    w &= 0xffff;
    u = (v - w) * (v + w);
    n ^= (((u >> 16) | (u << 16)) ^ 0xaa5835b9) + w * v;
    return  n;
}
#endif

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Ranfloat
//
// Purpose:   Random single-precision floating-point sample.
//            The 2^23 allowed values are odd multiples of 2^-24,
//            symmetrically placed in strict interior of (0,1).
//
// Notes: (1) Tuned to 23-bit mantissa in "float" format.
//        (2) Uses Ranint.
//
// History:   John Skilling   20 Oct 2002
//-----------------------------------------------------------------------------
// 
float Ranfloat(        //   O  Value
Rand_t   Rand)         // I O  Random generator state
{
    unsigned u = ((unsigned)Ranint(Rand)
                   & 0xfffffe00) ^ 0x00000100;  // switch lowest (2^-24) bit ON
    return (float)u * (float)SHIFT32;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Randouble
//
// Purpose:   Random double-precision floating-point sample.
//            The 2^52 allowed values are odd multiples of 2^-53,
//            symmetrically placed in strict interior of (0,1).
//
// Notes: (1) Tuned to 52-bit mantissa in "double" format.
//        (2) Uses one call to Ranint to get 64 random bits, with extra
//            random integer available in Rand[3].
//        (3) All floating-point random calls are directed through this code,
//            except Rangauss which uses the extra random integer in Rand[3].
//
// History:   John Skilling   6 May 1995, 3 Dec 1995, 24 Aug 1996
//                           20 Oct 2002, 17 Dec 2002
//-----------------------------------------------------------------------------
// 
double Randouble(       //   O  Value within (0,1)
Rand_t  Rand)           // I O  Random generator state
{
    unsigned  hi, lo;
    hi = (unsigned)Ranint(Rand);                // top 32 bits
    lo = (Rand[3]                               // low bits
                  & 0xfffff000) ^ 0x00000800;   // switch lowest (2^-53) bit ON
    return  ((double)hi + (double)lo * SHIFT32) * SHIFT32;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Rangauss
// 
// Purpose:   Sample from Gaussian N(0,1)
// 
// Notes: (1) Cannot overflow;  |value| < 10 standard deviations
//        (2) 2nd sample y*sqrt(...) could be available, best with y = ...+ 0.5
// 
// History:   JS          3 Oct 2002  Half of Box-Muller
//                       17 Dec 2002  Use Rand[3] for extra 40% speedup
//-----------------------------------------------------------------------------
// 
double Rangauss(       //   O  Value
Rand_t   Rand)         // I O  Random generator state
{
    static const double RR = 4611686018427387904.0;   // 2^62
    double x, y, r;
    do
    {
        x = Ranint(Rand) + 0.5;              // (-2^31,2^31), not 0
        y = (int)Rand[3];                    // [-2^31,2^31)
        r = x * x + y * y;
    } while( r >= RR );                      // disc, radius 2^31
    return  x * sqrt(2.0 * log(RR / r) / r); // Waste 2nd Box-Muller sample
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Rancauchy
// 
// Purpose:   Draw sample from
//                prob(x)  =  1 / (1 + x * x)
// 
// Notes:     Cannot overflow;  |value| <= 2^53/pi
// 
// History:   JS         3 Oct 2002
//-----------------------------------------------------------------------------
// 
double  Rancauchy(     //   O  Value
Rand_t   Rand)         // I O  Random generator state
{
    double  t = tan(HalfPi * (Randouble(Rand) - 0.5));
    return (1.0 - t * t) / (2.0 * t);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Rangamma
// 
// Purpose:   Draw sample from
//                             -1+c  -x
//                prob(x)  =  x     e  / Gamma(c)   for x in [0, inf)
// 
// Notes:     x=0 can be returned, especially if c is small, but x=inf cannot.
// 
// History:   JS         24 Jan 1994, 19 Oct 1995, 24 Aug 1996, 3 Oct 2002
//-----------------------------------------------------------------------------
// 
double  Rangamma(      //   O  Value
Rand_t   Rand,         // I O  Random generator state
double   c)            // I    Exponent
{
    double a, q, r, g, e;

    g = c - 1.0;
    if (g > 0.0)
    {
        e = sqrt(g + c) / g;
        do
        {
            do
            {
                r = Rancauchy(Rand);
                q = e * r + 1.0;
            } while (q <= 0.0);
        } while( log(Randouble(Rand) / (1.0 + r * r))
                > g * (log(q) - q + 1.0) );             // Accept?
        q *= g;
    }
    else
    {
        r = 1.0 / (1.0 + c * exp(-1.0));                // prob(bounding p < 1)
        do
        {
            if( Randouble(Rand) >= r )                  // Bounding x > 1
            {
                q = 1.0 - log(Randouble(Rand));         // exp(-q) in (1,inf)
                a = g * log(q);                         // log(Accept ratio)
            }
            else                                        // Bounding x <= 1
            {
                q = exp(log(Randouble(Rand)) / c);      // q^(-1+c) in (0,1)
                a = -q;                                 // log(Accept ratio)
            }
        } while( log(Randouble(Rand)) > a );
    }
    return q;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Ranpoiss
// 
// Purpose:   Draw sample from Poisson
//
//                           -c  j
//                prob(j) = e   c / j!               ( 0 <= c <~ 4000000000 )
// 
// Notes:     Distribution is truncated at 2^32-1,
//            so input parameter c should not approach 2^32.
//
// History:   JS         15 May 1998, 24 Mar 2001, 3 Oct 2002
//-----------------------------------------------------------------------------
// 
unsigned   Ranpoiss(   //   O  Value
Rand_t   Rand,         // I O  Random generator state
double   c)            // I    Mean
{
    unsigned j, k;
    double   a, b, r, t, y;

    if( c <= 0.0 )
        return 0;
    if( c < 70.0 )                             // Direct is faster for small c
    {
        j = k = (unsigned)c;
        t = exp(j * log(c) - c - logfactorial(j));
        for( ; ; )
        {
            r = Randouble(Rand);
            if( t >= r )
                return j;
            a = b = t;
            do
            {
                if( a > b )
                {
                    if( (t += a *= j-- / c) >= r )
                        return j;
                }
                else
                {
                    if( (t += b *= c / ++k) >= r )
                        return k;
                }
            } while( b > 0.0 );            // Trap failure from rounding error
        }
    }
    else
    {
        b = sqrt(2.0 * c);
        do
        {
            do
            {
                y = Rancauchy(Rand);
                r = c + b * y;
            } while( r < 0.0 );
            j = (unsigned)r;
            t = j * log(c) - logfactorial(j) - c;
        } while( log(Randouble(Rand) / (1.731 * (1.0 + y * y) * b)) > t );
        return j;
    }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Ranbinom
// 
// Purpose:   Draw sample from binomial
//
//                             n!     j      n-j
//                prob(j) = -------- p  (1-p)              0 <= j <= n
//                          j!(n-j)!
// 
// History:   JS         24 Mar 2001, 3 Oct 2002
//-----------------------------------------------------------------------------
// 
unsigned  Ranbinom(    //   O  Value
Rand_t   Rand,         // I O  Random generator state
unsigned n,            // I    Range
double   p)            // I    Mean/Range
{
    unsigned j;
    double   a, b, g, q, r, t, u, w, y;

    if( n == 0 )
        return 0;
    if( p <= 0.0 )
        return 0;
    if( p >= 1.0 )
        return n;
    q = 1.0 - p;

    if( n < 100 )
    {
        u = w = j = (unsigned)(p * n + p);
        t = exp(j * log(p) + (n-j) * log(q)
                 + logfactorial(n) - logfactorial(j) - logfactorial(n-j));
        for( ; ; )
        {
            r = Randouble(Rand);
            if( t >= r )
                return j;
            a = b = t;
            do
            {
                if( a > b )
                {
                    a *= u / p;
                    u -= 1.0;
                    a *= q / (n - u);
                    t += a;
                    if( t >= r )
                        return (unsigned)u;
                }
                else
                {
                    b *= (n - w) / q;
                    w += 1.0;
                    b *= p / w;
                    t += b;
                    if( t >= r )
                        return (unsigned)w;
                }
            } while( a + b > 0.0 );      // Trap failure from rounding error
        }
    }
    else
    {
        b = sqrt(2.0 * n * p * q + 0.25);
        g = logfactorial(n);
        do
        {
            do
            {
                y = Rancauchy(Rand);
                r = n * p + b * y + 0.5;
            } while( r < 0.0 || r >= n + 1.0 );
            j = (unsigned)r;
            t = g - logfactorial(j) - logfactorial(n - j)
               + j * log(p) + (n - j) * log(q);
        } while( log(Randouble(Rand) / (0.8190 * (1.0 + y * y) * b)) > t );
        return j;
    }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Ranbeta
//
// Purpose:   Draw sample from beta distribution
//                             r       n-r
//                prob(x)  =  x (1 - x)   / Z ,     Z = (n+1)! / r! (n-r)!
//
// Method:    Rejection beneath Cauchy of phenomenologically adequate width.
// 
// History:   John Skilling   28 Nov 2000
//-----------------------------------------------------------------------------
// 
double Ranbeta(       //   O  Value
Rand_t   Rand,        // I O  Random generator state
int      r,           // I    Number of successes >= 0
int      n)           // I    Total number        >= r
{
    double p, q, x, y, c;
    x = (double)n;
    if( r <= 0 )
        x = 1.0 - exp(log(Randouble(Rand)) / (x + 1.0));
    else if( r >= n )
        x = exp(log(Randouble(Rand)) / (x + 1.0));
    else
    {    
        p = (double)r / x;
        q = 1.0 - p;
        c = x / (3.0 * p * q);             // less than curvature at top
        do
        {
            do
            {
                y = Rancauchy(Rand) / sqrt(c);
                x = p + y;                 // x from Cauchy upper bound 
            }
            while( x <= 0.0 || x >= 1.0 ); // x within range
            y = r * log(x / p) + (n-r) * log((1.0 - x) / q)   // log(true)
               + log(1.0 + c * y * y);                        // -log(Cauchy)
        } while( log(Randouble(Rand)) > y );     // accept?
    }
    return x;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Rangeom
// 
// Purpose:   Draw integer from truncated geometric distribution 
//                                              j
//                    prob(j)   proportional   a    for j in [0,N-1]
//
//                    where  a = alpha/(alpha+1),   alpha > 0.0
//
//            Special case N = 0 means N = 2^32 for which
//
//                    <j> = alpha ,   var(j) = alpha(alpha+1)
// 
// History:   John Skilling   3 Apr 2001  Linearly truncated version
//                           17 Dec 2002  Hard truncation
//-----------------------------------------------------------------------------
// 
unsigned Rangeom(     //   O  Value
Rand_t   Rand,        // I O  Random generator state
double   alpha,       // I    Parameter
unsigned N)           // I    Supremum (0 interpreted as 2^32)
{
    double   a, r;

    if( N == 1 || alpha <= 0.0 )
        return 0;
    if( alpha * DBL_EPSILON >= 0.5 )
        return Rangrid(Rand, N);

    a = alpha / (alpha + 1.0);
    do
    {
        r = Randouble(Rand);
        if( N )
            r *= 1.0 - pow(a, (double)N);
        r = log(1.0 - r) / log(a);
    } while( r >= N );            // safety
    return (int)r;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Rangrid
//
// Purpose:   Random integer from [0,bound-1].
//
// Method:    Development of Ranint.
//
// History:   John Skilling   6 May 1995 - 15 Jan 2002
//-----------------------------------------------------------------------------
// 
unsigned   Rangrid(    //   O  Value
Rand_t   Rand,         // I O  Random generator state
unsigned bound)        // I    Supremum (0 interpreted as 2^32)
{
    unsigned   i, M, N;
    if( bound == 0 )
        return (unsigned)Ranint(Rand);
    N = (unsigned)(-1) / bound;
    M = N * bound;
    do  i = (unsigned)Ranint(Rand);
    while( i >= M );
    return i / N;
 }

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Ranperm
//
// Purpose:   Random permutation of {0,1,2,...,n-1}.
//
// History:   John Skilling   23 May 1998, 24 Jan 1999, 9 Nov 2000
//-----------------------------------------------------------------------------
// 
void Ranperm( 
Rand_t   Rand,        // I O  Random generator state
int      n,           // I    Dimension
int*     perm)        //   O  Output permutation
{
    int    i, m;
    for( m = 0; m < n; ++m )
    {
        i = Rangrid(Rand, (unsigned)(m + 1));
        if( i < m )
            perm[m] = perm[i];
        perm[i] = m;
    }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Ran2gauss
// 
// Purpose:   Draw sample from normalised form of
//
//                 exp( - g.x  -  x.A.x / 2 )             ( A non-negative )
//
// History:   JS         31 Dec 2001
//-----------------------------------------------------------------------------
// 
void Ran2gauss(
Rand_t   Rand,        // I O  Random generator state
double   g1,          // I    x gradient at origin
double   g2,          // I    y gradient at origin
double   A11,         // I    xx curvature >= 0
double   A12,         // I    xy curvature, A12*A12 <= A11*A22
double   A22,         // I    yy curvature >= 0
double*  x,           //   O  x sample position
double*  y)           //   O  y sample position
{
    double d = A11 * A22 - A12 * A12;
    double xbar = (g2 * A12 - g1 * A22) / d;
    double ybar = (g1 * A12 - g2 * A11) / d;
    *x = Rangauss(Rand) * sqrt(A22 / d);
    *y = Rangauss(Rand) / sqrt(A22) - *x * A12 / A22;
    *x += xbar;
    *y += ybar;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Ran1pos
// 
// Purpose:   Draw positive sample (x +ve) from normalised form of
//
//                   exp( - g x - A x^2 / 2 )
//
// History:   JS         30 Nov 1997
//-----------------------------------------------------------------------------
// 
double Ran1pos(    //   O  sample value
Rand_t  Rand,      // I O  random generator state
double  g,         // I    coeff of x
double  A)         // I    coeff of x^2,  >= 0
{
    double a, x;

    a = sqrt(A);
    if( 0.3 * a > g )                    // Reject from completed normal curve
    {
        do  x = Rangauss(Rand) / a - g / A;
        while( x < 0.0 );
    }
    else                                 // Reject from bounding exponential
    {
        do  x = -log(Randouble(Rand)) / g;
        while( -log(Randouble(Rand)) < A * x * x / 2.0 );
    }
    return x;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Ran2pos
// 
// Purpose:   Draw sample in positive quadrant from normalised form of
//
//              exp( - g1 x - g2 y - (A11 x^2 + 2 A12 x y + A22 y^2) / 2 )
//
// History:   JS         31 Dec 2001, 3 Oct 2002
//-----------------------------------------------------------------------------
// 
void Ran2pos(
Rand_t  Rand,     // I O  random generator state
double  g1,       // I    coeff of |x|,  > 0
double  g2,       // I    coeff of |y|,  > 0
double  A11,      // I    coeff of x^2,  >= 0
double  A12,      // I    coeff of xy,   |A12| < sqrt(A11*A22)
double  A22,      // I    coeff of y^2,  >= 0
double* x,        //   O  sample
double* y)        //   O  sample
{
    double  a;

    if( A12 * A12 <= DBL_EPSILON * A11 * A22 )
    {
        *x = Ran1pos(Rand, g1, A11);
        *y = Ran1pos(Rand, g2, A22);
        return;
    }

    g1 /= sqrt(A11);
    g2 /= sqrt(A22);
    a = A12 / sqrt(A11 * A22);
    if( a > 1.0 )
        a = 1.0;
    if( a < -1.0 )
        a = -1.0;

    if( g2 <= g1 )
        Positive2(Rand, g1, g2, -a, x, y);
    else
        Positive2(Rand, g2, g1, -a, y, x);

    *x /= sqrt(A11);
    *y /= sqrt(A22);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Positive2
// 
// Purpose:   Draw sample in positive quadrant from normalised form of
//
//                exp( - u x - v y - (x^2 - 2 a x y + y^2) / 2 )
//
//            with  u >= v  as well as  |a| < 1
//
// Method:
//   Apart from the simple special case  u > 0.3 or so, v > 0.3 or so:
//
//   Sample x from | exp -x^2/2  - u x + (a x - v)^2/2 ,  ax > v
//                 | exp -x^2/2  - u x                 ,  ax < v
//
//   then
//   Sample y from | exp -(y - a x + v)^2 / 2          ,  ax > v
//                 | exp - y^2 / 2                     ,  ax < v
//   with
//   Acceptance    | 1                                 ,  ax > v
//                 | exp (a x - v)y                    ,  ax < v
//   
//   Provided  u >= v, these algorithms always sample with O(1) efficiency.
//
//   General sampler for ANY logconcave function of x (such as the above) is
//            exp( height + ycap + q * (x - xcap) ),  x > xcap  (q < 0)
//            exp( height + ycap + p * (x - xcap) ),  x < xcap  (p > 0)
//   where (xcap,ycap) is apex of required curve over all x, height = O(1) +ve,
//   and p,q define left and right exponentials tangent to the required curve
//   A x^2 + B x + C (with A < 0), according to
//    p = Lslope(x,y, A,B,C) = 2*A*x + B + 2 * sqrt(-A *(y - A*x*x - B*x - C))
//    q = Rslope(x,y, A,B,C) = 2*A*x + B - 2 * sqrt(-A *(y - A*x*x - B*x - C))
//   This method always samples with O(1) efficiency.
// 
//   Subsidiary codes crafted to evaluate all sqrt arguments definitely +ve,
//   left slope p definitely +ve, and right slope q definitely -ve.
//        
// History:   JS         31 Dec 2001, 3 Oct 2002
//-----------------------------------------------------------------------------
// 
static void Positive2(
Rand_t  Rand,     // I O  random generator state
double  u,        // I    coeff of x
double  v,        // I    coeff of y
double  a,        // I    coeff of xy,   0 < |a| <= 1
double* x0,       //   O  sample
double* y0)       //   O  sample
{
    double height = 0.5000;
    double x, y, accept, xcap, p, q, r, s, dx;

    if( v > 0.3000 )
// Special case: sample around (0,0) from  exp - u x - v y
    {
        do
        {
            x = -log(Randouble(Rand)) / u;
            y = -log(Randouble(Rand)) / v;
            accept = a * x * y - 0.5 * (x * x + y * y);
        } while( accept < log(Randouble(Rand)) );
    }
    else
// General sampler ....
    {
// Bi-exponential for x
        for( ; ; )
        {
            if( a > 0.0 )
            {
                if( v < 0.0 )
                    Positive2A(u, v, a, height, &xcap, &p, &q);
                else
                    Positive2B(u, v, a, height, &xcap, &p, &q);
            }
            else
            {
                if( v < 0.0 )
                    Positive2C(u, v, a, height, &xcap, &p, &q);
                else
                    Positive2D(u,       height, &xcap, &p, &q);
            }

// Sample x
            for( ; ; )
            {
                s = (xcap == 0.0 || Randouble(Rand) < p / (p - q)) ? q : p;
                dx = log(Randouble(Rand)) / s;
                x = dx + xcap;
                if( x <= 0.0 )
                    continue;
                accept = -dx * (u + xcap + 0.5 * dx);
                r = a * x - v;
                if( r > 0.0 )
                    accept += a * dx * (0.5 * a * (x + xcap) - v);
                accept -= height + s * dx;
                if( log(Randouble(Rand)) < accept )
                    break;
            }

// Sample y|x
            if( - v + a * x > 0.0 )
            {
                y = - v + a * x + Rangauss(Rand);
                if( y >= 0.0 )
                    break;
            }
            else
            {
                y = Rangauss(Rand);
                if( y >= 0.0 && log(Randouble(Rand)) < (- v + a * x) * y )
                    break;
            }
        }
    }
// Exit
    *x0 = x;
    *y0 = y;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Positive2A
// 
// Purpose:   Set bi-exponential sampler for Positive2 in case a > 0, v <= 0
//      x = inf       |      
//                    | exp  - x^2 / 2 - u x + (a x - v)^2 / 2
//      x = 0         |
//                 
// History:   JS         31 Dec 2001, 3 Oct 2002
//-----------------------------------------------------------------------------
// 
static void Positive2A(
double  u,        // I    >= v
double  v,        // I    <= 0
double  a,        // I    > 0
double  height,   // I    height of sampler cusp above target maximum
double* xc,       //   O  target maximum
double* p0,       //   O  sampler  exp( yc + p0 * (x-xc) ),  x < xc  (p0 > 0)
double* q0)       //   O  sampler  exp( yc + q0 * (x-xc) ),  x > xc  (q0 < 0)
{
    double  y0, m0, xcap, ycap;
    double  d, p, q;
    d = 1.0 - a * a;
    if( d < DBL_EPSILON )
        d = DBL_EPSILON;

// Set boundaries  (0,y0)...
// and slopes        m0
    y0 = 0.5 * v * v;
    m0 = -(u + a * v);

// Find apex from slope=0 
    if( m0 <= 0.0 )                      // apex at 0
    {
        xcap = 0.0;
        ycap = height + y0;

        p = 0.0;   // arbitrary
        q = m0 - sqrt(2.0 * d * height); // Rslope(xcap, ycap, -d/2, m0, v*v/2)
    }
    else                                 // apex in (0,inf)
    {
        xcap = m0 / d;
        ycap = height + 0.5 * (m0 * xcap + v * v);

        if( ycap - y0 <= m0 * xcap )
            p = sqrt(2.0 * d * height);  // Lslope(xcap, ycap, -d/2, m0, v*v/2)
        else
            p = (ycap - y0) / xcap;
        q = -sqrt(2.0 * d * height);     // Rslope(xcap, ycap, -d/2, m0, v*v/2)
    }

    *xc = xcap;
    *p0 = p;
    *q0 = q;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Positive2B
// 
// Purpose:   Set bi-exponential sampler for Positive2 in case a > 0, v >= 0
//                 
//                    | exp  - x^2 / 2 - u x + (a x - v)^2 / 2
//     x2 = v/a       |
//                    | exp  - x^2 / 2 - u x
//      x = 0         |
//                 
// History:   JS         31 Dec 2001, 3 Oct 2002
//-----------------------------------------------------------------------------
// 
static void Positive2B(
double  u,        // I    >= v
double  v,        // I    >= 0
double  a,        // I    > 0
double  height,   // I    height of sampler cusp above target maximum
double* xc,       //   O  target maximum
double* p0,       //   O  sampler  exp( yc + p0 * (x-xc) ),  x < xc  (p0 > 0)
double* q0)       //   O  sampler  exp( yc + q0 * (x-xc) ),  x > xc  (q0 < 0)
{
    double  x2, y2, m2;
    double  d, q, r;
    d = 1.0 - a * a;
    if( d < DBL_EPSILON )
        d = DBL_EPSILON;

// Set boundaries  (0,y0)...(x2,y2)...
// and slopes     +ve,-ve      m2
    x2 = v / a;
    y2 = x2 * (- u - 0.5 * x2);
    m2 = - u - x2;

// Apex necessarily at 0 
    r = (height - y2) + m2 * x2;
    if( r <= 0.0 )                    // Rslope(0, height, -1/2., -u, 0)
        q = -u - sqrt(2.0 * height);
    else                              // Rslope(0, height, -d/2, -u-a*v, v*v/2)
        q = -(u + a * v) - sqrt(2.0 * d * (d * height + a * a * r));

    *xc = 0.0;
    *p0 = 0.0;   // arbitrary
    *q0 = q;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Positive2C
// 
// Purpose:   Set bi-exponential sampler for Positive2 in case a < 0, v <= 0
//                 
//                    | exp  - x^2 / 2 - u x
//     x1 = -v/-a     |
//                    | exp  - x^2 / 2 - u x + (a x - v)^2 / 2
//      x = 0         |
//                 
// History:   JS         31 Dec 2001, 3 Oct 2002
//-----------------------------------------------------------------------------
// 
static void Positive2C(
double  u,        // I    >= v
double  v,        // I    <= 0
double  a,        // I    < 0
double  height,   // I    height of sampler cusp above target maximum
double* xc,       //   O  target maximum
double* p0,       //   O  sampler  exp( yc + p0 * (x-xc) ),  x < xc  (p0 > 0)
double* q0)       //   O  sampler  exp( yc + q0 * (x-xc) ),  x > xc  (q0 < 0)
{
    double  x1, y0, y1, m0, m1, xcap, ycap;
    double  d, p, q, r;
    d = 1.0 - a * a;
    if( d < DBL_EPSILON )
        d = DBL_EPSILON;

// Set boundaries  (0,y0)...(x1,y1)...
// and slopes         m0      m1
    x1 = v / a;
    y1 = x1 * (- u - 0.5 * x1);
    m1 = - u - x1;                   // -ve
    y0 = 0.5 * v * v;
    m0 = -(u + a * v);

// Find apex from slope=0 
    if( m0 <= 0.0 )                  // apex at 0
    {
        xcap = 0.0;
        ycap = height + 0.5 * v * v;

        p = 0.0;   // arbitrary
        r = (ycap - y1) + m1 * x1;
        if( r <= 0.0 )               // Rslope(0, ycap, -d/2, -u-a*v, v*v/2)
            q = m0 - sqrt(2.0 * d * height);
        else                         // Rslope(0, ycap, -1/2., -u, 0)
            q = (v - u) + 2.0 * height / (v - sqrt(2.0 * ycap));
    }
    else                             // apex in (0,x1)
    {
        xcap = m0 / d;
        ycap = height + 0.5 * (d * xcap * xcap + v * v);

        r = (ycap - y0) - m0 * xcap;
        if( r <= 0.0 )               // Lslope(xcap, ycap, -d/2, -u-a*v, v*v/2)
            p = sqrt(2.0 * d * height);
        else
            p = (ycap - y0) / xcap;
        r = (ycap - y1) - m1 * (xcap - x1);
        if( r <= 0.0 )               // Rslope(xcap, ycap, -d/2, -u-a*v, v*v/2)
            q = -sqrt(2.0 * d * height);
        else                         // Rslope(xcap, ycap, -1/2., -u, 0)
            q = (a * (u - v) + (1.0 + a) * v) * a / d
                 - sqrt(2.0 * r + (xcap - x1) * (xcap - x1));
    }

    *xc = xcap;
    *p0 = p;
    *q0 = q;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Positive2D
// 
// Purpose:   Set bi-exponential sampler for Positive2 in case a < 0, v >= 0
//      x = inf       |
//                    | exp  - x^2 / 2 - u x
//      x = 0         |
//                 
// History:   JS         31 Dec 2001, 3 Oct 2002
//-----------------------------------------------------------------------------
// 
static void Positive2D(
double  u,        // I    >= v >= 0
double  height,   // I    height of sampler cusp above target maximum
double* xc,       //   O  target maximum
double* p0,       //   O  sampler  exp( yc + p0 * (x-xc) ),  x < xc  (p0 > 0)
double* q0)       //   O  sampler  exp( yc + q0 * (x-xc) ),  x > xc  (q0 < 0)
{

// Apex necessarily at 0 
    *xc = 0.0;
    *p0 = 0.0;   // arbitrary
    *q0 = -u - sqrt(2.0 * height);           // Rslope(0, height, -1/2., -u, 0)
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Ran1posneg
// 
// Purpose:   Draw sample from normalised form of
//
//                   exp( - u|x| - f x - A x^2 / 2 )
//
// Method:
//   General sampler for ANY logconcave function of x (such as the above) is
//            exp( height + ycap + q * (x - xcap) ),  x > xcap  (q < 0)
//            exp( height + ycap + p * (x - xcap) ),  x < xcap  (p > 0)
//   where (xcap,ycap) is apex of required curve over all x, height = O(1) +ve,
//   and p,q define left and right exponentials tangent to the required curve
//   a x^2 + b x + c (with A < 0), according to
//    p = Lslope(x,y, a,b,c) = 2*a*x + b + 2 * sqrt(-a *(y - a*x*x - b*x - c))
//    q = Rslope(x,y, a,b,c) = 2*a*x + b - 2 * sqrt(-a *(y - a*x*x - b*x - c))
//   This method always samples with O(1) efficiency.
//
// History:   JS         31 Dec 2001, 3 Oct 2000
//-----------------------------------------------------------------------------
// 
double Ran1posneg( //   O  sample value
Rand_t  Rand,      // I O  random generator state
double  f,         // I    coeff of x
double  u,         // I    coeff of |x|,  > 0
double  A)         // I    coeff of x^2,  >= 0
{
static const double height = 0.5000;
    double  x, xcap, ycap;
    double  accept, p, q, r, m9, m0;

// Special case of no curvature
    if( A <= 0.0 )
    {
        if( Randouble(Rand) < 0.5 * (1.0 + f / u) )
            return log(Randouble(Rand)) / (u - f);
        else
            return -log(Randouble(Rand)) / (u + f);
    }
// Scale general case to unit curvature
    f /= sqrt(A);
    u /= sqrt(A);
    m9 = u - f;
    m0 = -u - f;
// Set dominating bi-exponential sampler
    if( m9 < 0.0 )          // Apex in x -ve
    {
        xcap = m9;
        ycap = height + 0.5 * xcap * xcap;
        p = sqrt(2.0 * height);             // Lslope(xcap, ycap, -1/2., m9, 0)
        if( ycap <= m9 * xcap )
            q = -sqrt(2.0 * height);        // Rslope(xcap, ycap, -1/2., m9, 0)
        else
        {
            r = ycap - m0 * xcap;
            if( r <= 0.0 )
                q = ycap / xcap;
            else                            // Rslope(xcap, ycap, -1/2., m0, 0)
                q = -2.0 * u - sqrt(2.0 * r + m9 * m9);
        }
    }
    else if( m0 <= 0.0 )   // Apex at x = 0
    {
        xcap = 0.0;
        ycap = height;
        p = m9 + sqrt(2.0 * height);        // Lslope(xcap, ycap, -1/2., m9, 0)
        q = m0 - sqrt(2.0 * height);        // Rslope(xcap, ycap, -1/2., m0, 0)
    }
    else                   // Apex in x +ve
    {
        xcap = m0;
        ycap = height + 0.5 * xcap * xcap;
        if( ycap <= m0 * xcap )
            p = sqrt(2.0 * height);         // Lslope(xcap, ycap, -1/2., m0, 0)
        else
        {
            r = ycap - m9 * xcap;
            if( r <= 0.0 )
                p = ycap / xcap;
            else                            // Lslope(xcap, ycap, -1/2., m9, 0)
                p = 2.0 * u + sqrt(2.0 * r + m0 * m0);
        }
        q = -sqrt(2.0 * height);            // Rslope(xcap, ycap, -1/2., m0, 0)
    }
// Sample x
    do
    {
// accept = -log(sampler)
        if( Randouble(Rand) < p / (p - q) )
        {
            x = xcap + log(Randouble(Rand)) / q;
            accept = -(ycap + q * (x - xcap));
        }
        else
        {
            x = xcap + log(Randouble(Rand)) / p;
            accept = -(ycap + p * (x - xcap));
        }
// accept += log(true)  
        if( x > 0.0 )
            accept += x * (- u - f - 0.5 * x);
        else
            accept += x * (u - f - 0.5 * x);
    } while( log(Randouble(Rand)) > accept );
// Exit
    return x / sqrt(A);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Ran2posneg
// 
// Purpose:   Draw sample from normalised form of
//
//    exp( - u|x| - v|y| - f x - g y - (A11 x^2 + 2 A12 x y + A22 y^2) / 2 )
//
// History:   JS         31 Dec 2001, 3 Oct 2002
//-----------------------------------------------------------------------------
// 
void Ran2posneg(
Rand_t  Rand,     // I O  random generator state
double  f,        // I    coeff of x
double  g,        // I    coeff of y
double  u,        // I    coeff of |x|,  > 0
double  v,        // I    coeff of |y|,  > 0
double  A11,      // I    coeff of x^2,  >= 0
double  A12,      // I    coeff of xy,   |A12| < sqrt(A11*A22)
double  A22,      // I    coeff of y^2,  >= 0
double* x,        //   O  sample
double* y)        //   O  sample
{
    double  a;
    if( A12 * A12 <= DBL_EPSILON * A11 * A22 )
    {
        *x = Ran1posneg(Rand, f, u, A11);
        *y = Ran1posneg(Rand, g, v, A22);
        return;
    }
    f /= sqrt(A11);
    u /= sqrt(A11);
    g /= sqrt(A22);
    v /= sqrt(A22);
    a = A12 / sqrt(A11 * A22);
    if( a > 1.0 )
        a = 1.0;
    if( a < -1.0 )
        a = -1.0;

    if( f <= 0.0 )
    {
        if( g <= 0.0 )
        {
            if( -g - v >= -f - u )
                Posneg2(Rand, -f, -g, u, v, -a, x, y);
            else
                Posneg2(Rand, -g, -f, v, u, -a, y, x);
        }
        else
        {
            if( g - v >= -f - u )
                Posneg2(Rand, -f, g, u, v, a, x, y);
            else
                Posneg2(Rand, g, -f, v, u, a, y, x);
            *y = - *y;
        }
    }
    else
    {
        if( g <= 0.0 )
        {
            if( -g - v >= f - u )
                Posneg2(Rand, f, -g, u, v, a, x, y);
            else
                Posneg2(Rand, -g, f, v, u, a, y, x);
        }
        else
        {
            if( g - v >= f - u )
                Posneg2(Rand, f, g, u, v, -a, x, y);
            else
                Posneg2(Rand, g, f, v, u, -a, y, x);
            *y = - *y;
        }
        *x = - *x;
    }

    *x /= sqrt(A11);
    *y /= sqrt(A22);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Posneg2
// 
// Purpose:   Draw sample from normalised form of
//
//        exp( - u|x| - v|y| + f x + g y - (x^2 - 2 a x y + y^2) / 2 )
//
//            with  f >= 0, g >= 0, f - u <= g - v
//            as well as u > 0, v > 0, 0 < |a| <= 1
//
// Method:
//   Apart from the simple special case  u - f > 0.3 or so, v - g > 0.3 or so:
//
//                 | exp -x^2/2 + f x - u|x| + (g - v + a x)^2/2 ,  g-v+ax > 0
//   Sample x from | exp -x^2/2 + f x - u|x|                     ,  |g+ax| < v
//                 | exp -x^2/2 + f x - u|x| + (g + v + a x)^2/2 ,  g+v+ax < 0
//   then
//                 | exp -(y - g + v - a x)^2 / 2                ,  g-v+ax > 0
//   Sample y from | exp - y^2 / 2                               ,  |g+ax| < v
//                 | exp -(y - g - v - a x)^2 / 2                ,  g+v+ax < 0
//   with
//                 | exp - v (|y| - y)                           ,  g-v+ax > 0
//   Acceptance    | exp - v |y| + (g + a x) y                   ,  |g+ax| < v
//                 | exp - v (|y| + y)                           ,  g+v+ax < 0
//
//   Provided  f-u <= g-v, these algorithms always sample with O(1) efficiency.
//
//   General sampler for ANY logconcave function of x (such as the above) is
//            exp( height + ycap + q * (x - xcap) ),  x > xcap  (q < 0)
//            exp( height + ycap + p * (x - xcap) ),  x < xcap  (p > 0)
//   where (xcap,ycap) is apex of required curve over all x, height = O(1) +ve,
//   and p,q define left and right exponentials tangent to the required curve
//   A x^2 + B x + C (with A < 0), according to
//    p = Lslope(x,y, A,B,C) = 2*A*x + B + 2 * sqrt(-A *(y - A*x*x - B*x - C))
//    q = Rslope(x,y, A,B,C) = 2*A*x + B - 2 * sqrt(-A *(y - A*x*x - B*x - C))
//   This method always samples with O(1) efficiency.
// 
//   Subsidiary codes crafted to evaluate all sqrt arguments definitely +ve,
//   left slope p definitely +ve, and right slope q definitely -ve.
//        
// History:   JS         31 Dec 2001, 3 Oct 2002
//-----------------------------------------------------------------------------
// 
static void Posneg2(
Rand_t  Rand,     // I O  random generator state
double  f,        // I    coeff of x,    >= 0
double  g,        // I    coeff of y,    >= 0
double  u,        // I    coeff of |x|,  > 0
double  v,        // I    coeff of |y|,  > 0
double  a,        // I    coeff of xy,   -1 < a < 1
double* x0,       //   O  sample
double* y0)       //   O  sample
{
    double height = 0.5000;
    double x, y, accept, xcap, p, q, r, s, dx;

    if( v - g > 0.3000 )
// Special case: sample around (0,0) from exp  f x + g y - u|x| - v|y|
    {
        do
        {
            x = (Randouble(Rand) > 0.5 * (1.0 - f / u))
               ? log(Randouble(Rand)) / (f - u)
               : log(Randouble(Rand)) / (f + u);
            y = (Randouble(Rand) > 0.5 * (1.0 - g / v))
               ? log(Randouble(Rand)) / (g - v)
               : log(Randouble(Rand)) / (g + v);
            accept = a * x * y - 0.5 * (x * x + y * y);
        } while( accept < log(Randouble(Rand)) );
        *x0 = x;
        *y0 = y;
        return;
    }

// General bi-exponential for x
    do
    {
        if( a > 0.0 )
        {
            if( g > v )
                Posneg2A(f+u, f-u, g+v, g-v, a, height, &xcap, &p, &q);
            else
                Posneg2B(f+u, f-u, g+v, g-v, a, height, &xcap, &p, &q);
        }
        else
        {
            if( g > v )
                Posneg2C(f+u, f-u, g+v, g-v, a, height, &xcap, &p, &q);
            else
                Posneg2D(f+u, f-u, g+v, g-v, a, height, &xcap, &p, &q);
        }

// Sample x
        do
        {
            s = (Randouble(Rand) < p / (p - q)) ? q : p;
            dx = log(Randouble(Rand)) / s;
            x = dx + xcap;
            if( x > 0.0 )
                accept = dx * (f - u - 0.5 * (xcap + x));
            else
                accept = dx * (f + u - 0.5 * (xcap + x));
            r = g - v + a * x;
            if( r > 0.0 )
                accept += a * dx * (r - 0.5 * a * dx);
            r = g + v + a * x;
            if( r < 0.0 )
                accept += a * dx * (r - 0.5 * a * dx);
            accept -= height + s * dx;
        } while( log(Randouble(Rand)) > accept );

// Sample y|x
        if( g - v + a * x > 0.0 )
        {
            y = g - v + a * x + Rangauss(Rand);
            accept = (y >= 0.0) ? 0.0 : 2.0 * v * y; 
        }
        else if( g + v + a * x < 0.0 )
        {
            y = g + v + a * x + Rangauss(Rand);
            accept = (y > 0.0) ? -2.0 * v * y : 0.0;
        }
        else
        {
            y = Rangauss(Rand);
            accept = (g + a * x) * y - v * fabs(y);
        }
    } while( accept < 0.0 && accept < log(Randouble(Rand)) );
// Exit
    *x0 = x;
    *y0 = y;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Posneg2A
// 
// Purpose:   Set bi-exponential sampler for Posneg2 in case a > 0, g >= v
//                 
//                    | exp  - x^2 / 2 + f x - u x + (g - v + a x)^2 / 2
//      x = 0         |
//                    | exp  - x^2 / 2 + f x + u x + (g - v + a x)^2 / 2
//     x2 = -(g-v)/a  |
//                    | exp  - x^2 / 2 + f x + u x
//     x1 = -(g+v)/a  |
//                    | exp  - x^2 / 2 + f x + u x + (g + v + a x)^2 / 2
//         
// History:   JS         31 Dec 2001, 3 Oct 2002
//-----------------------------------------------------------------------------
// 
static void Posneg2A(
double  fpu,      // I    f + u  > 0
double  fmu,      // I    f - u  <= gmv
double  gpv,      // I    g + v  > 0
double  gmv,      // I    g - v  >= 0
double  a,        // I    > 0
double  height,   // I    height of sampler cusp above target maximum
double* xc,       //   O  target maximum
double* p0,       //   O  sampler  exp( yc + p0 * (x-xc) ),  x < xc  (p0 > 0)
double* q0)       //   O  sampler  exp( yc + q0 * (x-xc) ),  x > xc  (q0 < 0)
{
    double  x1, x2, y0, y1, y2, m0, m1, m2, m9, xcap, ycap;
    double  d, p, q, r, s, t, u;
    d = 1.0 - a * a;
    if( d < DBL_EPSILON )
        d = DBL_EPSILON;

// Set boundaries  ...(x1,y1)...(x2,y2)...(0,y0)...
// and slopes            m1        m2     m9,m0
    x1 = -gpv / a;
    y1 = x1 * (fpu - 0.5 * x1);
    m1 = fpu - x1;

    x2 = -gmv / a;
    y2 = x2 * (fpu - 0.5 * x2);
    m2 = fpu - x2;

    y0 = 0.5 * gmv * gmv;
    m9 = fpu + a * gmv;                 // +ve
    m0 = fmu + a * gmv;

// Find apex from slope=0 
    if( m0 <= 0.0 )                     // apex at 0
    {
        xcap = 0.0;
        ycap = height + y0;

        if(ycap <= y2 - m2*x2)    // Lslope(xcap,ycap,-d/2,fpu+a*gmv,gmv*gmv/2)
            p = m9 + sqrt(2.0 * d * height);
        else
        {
            r = (ycap - y1) + m1 * x1;
            if( r <= 0.0 )        // Lslope(xcap, ycap, -1/2., fpu, 0)
                p = fpu + sqrt(2.0 * ycap);
            else                  // Lslope(xcap,ycap,-d/2,fpu+a*gpv,gpv*gpv/2)
                p = (fpu + a * gpv) + sqrt(2.0 * d * (d * ycap + a * a * r));
        }                         // Rslope(xcap,ycap,-d/2,fmu+a*gmv,gmv*gmv/2)
        q = m0 - sqrt(2.0 * d * height);
    }
    else                                // apex in (0,inf)
    {
        xcap = m0 / d;
        ycap = height + 0.5 * (d * xcap * xcap + gmv * gmv);

        if(ycap - y0 <= m0*xcap)  // Lslope(xcap,ycap,-d/2,fmu+a*gmv,gmv*gmv/2)
            p = sqrt(2.0 * d * height);
        else
        {
            r = (ycap - y0) - m9 * xcap;
            if( r <= 0.0 )
                p = (ycap - y0) / xcap;
            else
            {
                s = (ycap - y2) - m2 * (xcap - x2);
                if( s <= 0.0 )    // Lslope(xcap,ycap,-d/2,fpu+a*gmv,gmv*gmv/2)
                    p = (fpu - fmu) + sqrt(d * (2.0 * r + d * xcap * xcap));
                else
                {
                    t = (ycap - y1) - m1 * (xcap - x1);
                    if( t <= 0.0 )
                    {             // Lslope(xcap,ycap,-1/2.,fpu,0)
                        u = 2.0 * s - 2.0 * xcap * x2 + x2 * x2;
                        p = fpu + u / (xcap + sqrt(u + xcap * xcap));
                    }
                    else
                    {             // Lslope(xcap,ycap,-d/2,fpu+a*gpv,gpv*gpv/2)
                        u = fmu + a * gmv + d * gpv / a;
                        p = ((fpu - fmu) + a * (gpv - gmv))
                                            + sqrt(2.0 * d * t + u * u);
                    }
                }
            }
        }                         // Rslope(xcap,ycap,-d/2,fmu+a*gmv,gmv*gmv/2)
        q = -sqrt(2.0 * d * height);
    }

    *xc = xcap;
    *p0 = p;
    *q0 = q;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Posneg2B
// 
// Purpose:   Set bi-exponential sampler for Posneg2 in case a > 0, g <= v
//                 
//                    | exp  - x^2 / 2 + f x - u x + (g - v + a x)^2 / 2
//     x2 = (v-g)/a   |
//                    | exp  - x^2 / 2 + f x - u x
//      x = 0         |
//                    | exp  - x^2 / 2 + f x + u x
//     x1 = -(v+g)/a  |
//                    | exp  - x^2 / 2 + f x + u x + (g + v + a x)^2 / 2
//                 
// History:   JS         31 Dec 2001, 3 Oct 2002
//-----------------------------------------------------------------------------
// 
static void Posneg2B(
double  fpu,      // I    f + u  > 0
double  fmu,      // I    f - u  <= gmv
double  gpv,      // I    g + v  > 0
double  gmv,      // I    g - v  <= 0
double  a,        // I    > 0
double  height,   // I    height of sampler cusp above target maximum
double* xc,       //   O  target maximum
double* p0,       //   O  sampler  exp( yc + p0 * (x-xc) ),  x < xc  (p0 > 0)
double* q0)       //   O  sampler  exp( yc + q0 * (x-xc) ),  x > xc  (q0 < 0)
{
    double  x1, x2, y1, y2, m1, m2;
    double  d, p, q, t;
    d = 1.0 - a * a;
    if( d < DBL_EPSILON )
        d = DBL_EPSILON;

// Set boundaries  ...(x1,y1)...(0,y0)...(x2,y2)...
// and slopes            m1    +ve,-ve      m2
    x1 = -gpv / a;
    y1 = x1 * (fpu - 0.5 * x1);
    m1 = fpu - x1;

    x2 = -gmv / a;
    y2 = x2 * (fmu - 0.5 * x2);
    m2 = fmu - x2;

// Apex necessarily at 0 
    t = m1 * x1 - (y1 - height);
    if( t <= 0.0 )                 // Lslope(0,height,-1/2.,fpu,0)
        p = fpu + sqrt(2.0 * height);
    else                           // Lslope(0,height,-d/2,fpu+a*gpv,gpv*gpv/2)
        p = (fpu + a * gpv) + sqrt(d * (2.0 * t + d * x1 * x1));

    t = m2 * x2 - (y2 - height);
    if( t <= 0.0 )                 // Rslope(0,height,-1/2.,fmu,0)
        q = fmu - sqrt(2.0 * height);
    else                           // Rslope(0,height,-d/2,fmu+a*gmv,gmv*gmv/2)
        q = (fmu + a * gmv) - sqrt(d * (2.0 * t + d * x2 * x2));

    *xc = 0.0;
    *p0 = p;
    *q0 = q;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Posneg2C
// 
// Purpose:   Set bi-exponential sampler for Posneg2 in case a < 0, g >= v
//                 
//                    | exp  - x^2 / 2 + f x - u x + (g + v + a x)^2 / 2
//     x2 = (g+v)/-a  |
//                    | exp  - x^2 / 2 + f x - u x
//     x1 = (g-v)/-a  |
//                    | exp  - x^2 / 2 + f x - u x + (g - v + a x)^2 / 2
//      x = 0         |
//                    | exp  - x^2 / 2 + f x + u x + (g - v + a x)^2 / 2
//                 
// History:   JS         31 Dec 2001, 3 Oct 2002
//-----------------------------------------------------------------------------
// 
static void Posneg2C(
double  fpu,      // I    f + u  > 0
double  fmu,      // I    f - u  <= gmv
double  gpv,      // I    g + v  > 0
double  gmv,      // I    g - v  >= 0
double  a,        // I    < 0
double  height,   // I    height of sampler cusp above target maximum
double* xc,       //   O  target maximum
double* p0,       //   O  sampler  exp( yc + p0 * (x-xc) ),  x < xc  (p0 > 0)
double* q0)       //   O  sampler  exp( yc + q0 * (x-xc) ),  x > xc  (q0 < 0)
{
    double  x1, x2, y0, y1, y2, m0, m1, m2, m9, xcap, ycap;
    double  d, p, q, r, s, t, u, v, w;
    d = 1.0 - a * a;
    if( d < DBL_EPSILON )
        d = DBL_EPSILON;

// Set boundaries  ...(0,y0)...(x1,y1)...(x2,y2)...
// and slopes         m9,m0      m1        m2
    x1 = -gmv / a;
    y1 = x1 * (fmu - 0.5 * x1);
    m1 = fmu - x1;                      // -ve

    x2 = -gpv / a;
    y2 = x2 * (fmu - 0.5 * x2);
    m2 = fmu - x2;

    y0 = 0.5 * gmv * gmv;
    m9 = fpu + a * gmv;
    m0 = fmu + a * gmv;

// Find apex from slope=0 
    if( m9 < 0.0 )                      // apex in (-inf,0)
    {
        xcap = m9 / d;
        ycap = height + 0.5 * (d * xcap * xcap + gmv * gmv);

        p = sqrt(2. * d * height);// Lslope(xcap,ycap,-d/2,fpu+a*gmv,gmv*gmv/2)
        r = (ycap - y0) - m9 * xcap;
        if( r <= 0.0 )            // Rslope(xcap,ycap,-d/2,fpu+a*gmv,gmv*gmv/2)
            q = -sqrt(2.0 * d * height);
        else
        {
            s = ycap - y0 - m0 * xcap;
            if( s <= 0.0 )
                q = (ycap - y0) / xcap;
            else
            {
                t = (ycap - y1) - m1 * (xcap - x1);
                if( t <= 0.0 )    // Rslope(xcap,ycap,-d/2,fmu+a*gmv,gmv*gmv/2)
                    q = (fmu - fpu) - sqrt(2.0 * d * s + m9 * m9);
                else
                {
                    u = (ycap - y2) - m2 * (xcap - x2);
                    if( u <= 0.0 )
                    {             // Rslope(xcap,ycap,-1/2.,fmu,0)
                        v = 2.0 * t + (x1 - xcap) * (x1 - xcap);
                        w = (x1 * (1.0 - a) - 2.0 * xcap) * x1 * (1.0 + a);
                        q = (fmu - gmv)
                           - (2.0 * t + w) / (gmv - xcap + sqrt(v));
                    }
                    else
                    {             // Rslope(xcap,ycap,-d/2,fmu+a*gpv,gpv*gpv/2)
                        v = d * (2.0 * u + d * (x2 - xcap) * (x2 - xcap));
                        w = (fmu - fpu) - a * (gmv - gpv);
                        q = w - sqrt(v);
                    }
                }
            }
        }
    }
    else if( m0 <= 0.0 )                 // apex at 0
    {
        xcap = 0.0;
        ycap = height + 0.5 * gmv * gmv;
                                  // Lslope(xcap,ycap,-d/2,fpu+a*gmv,gmv*gmv/2)
        p = m9 + sqrt(2.0 * d * height);

        r = ycap - y1 + m1 * x1;
        if( r <= 0.0 )            // Rslope(xcap,ycap,-d/2,fmu+a*gmv,gmv*gmv/2)
            q = m0 - sqrt(2.0 * d * height);
        else
        {
            s = ycap - y2 + m2 * x2;
            if( s <= 0.0 )        // Rslope(xcap,ycap,-1/2.,fmu,0)
                q = (fmu - gmv) - 2.0 * height / (gmv + sqrt(2.0 * ycap));
            else                  // Rslope(xcap,ycap,-d/2,fmu+a*gpv,gpv*gpv/2)
                q = (a * (gpv - gmv) + m0) - sqrt(d * (2.0 * s + d * x2 * x2));
        }
    }
    else                                // apex in (0,x1)
    {
        xcap = m0 / d;
        ycap = height + 0.5 * (d * xcap * xcap + gmv * gmv);

        r = ycap - y0 - m0 * xcap;
        if( r <= 0.0 )            // Lslope(xcap,ycap,-d/2,fmu+a*gmv,gmv*gmv/2)
            p = sqrt(2.0 * d * height);
        else
        {
            s = ycap - y0 - m9 * xcap;
            if( s <= 0.0 )
                p = (ycap - y0) / xcap;
            else                  // Lslope(xcap,ycap,-d/2,fpu+a*gmv,gmv*gmv/2)
                p = (fpu - fmu) + sqrt(d * (2.0 * s + d * xcap * xcap));
        }

        r = ycap - y1 - m1 * (xcap - x1);
        if( r <= 0.0 )            // Rslope(xcap,ycap,-d/2,fmu+a*gmv,gmv*gmv/2)
            q = -sqrt(2.0 * d * height);
        else
        {
            s = ycap - y2 - m2 * (xcap - x2);
            if( s <= 0.0 )
            {                     // Rslope(xcap,ycap,-1/2.,fmu,0)
                v = 2.0 * r + (fmu - x1) * (fmu - x1) * (1.0 + a * a) / d;
                w = (a * (gmv - fmu) - (1.0 + a) * gmv) * a / d;
                q = -v / (w + sqrt(2.0 * r + (xcap - x1) * (xcap - x1)));
            }
            else                  // Rslope(xcap,ycap,-d/2,fmu+a*gpv,gpv*gpv/2)
                q = a * (gpv - gmv)
                   - sqrt(d * (2.0 * s + d * (xcap - x2) * (xcap - x2)));
        }
    }

    *xc = xcap;
    *p0 = p;
    *q0 = q;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Posneg2D
// 
// Purpose:   Set bi-exponential sampler for Posneg2 in case a < 0, g <= v
//                 
//                    | exp  - x^2 / 2 + f x - u x + (g + v + a x)^2 / 2
//     x2 = (v+g)/-a  |
//                    | exp  - x^2 / 2 + f x - u x
//      x = 0         |
//                    | exp  - x^2 / 2 + f x + u x
//     x1 = -(v-g)/-a |
//                    | exp  - x^2 / 2 + f x + u x + (g - v + a x)^2 / 2
//                 
// History:   JS         31 Dec 2001, 3 Oct 2002
//-----------------------------------------------------------------------------
// 
static void Posneg2D(
double  fpu,      // I    f + u  > 0
double  fmu,      // I    f - u  <= gmv
double  gpv,      // I    g + v  > 0
double  gmv,      // I    g - v  <= 0
double  a,        // I    < 0
double  height,   // I    height of sampler cusp above target maximum
double* xc,       //   O  target maximum
double* p0,       //   O  sampler  exp( yc + p0 * (x-xc) ),  x < xc  (p0 > 0)
double* q0)       //   O  sampler  exp( yc + q0 * (x-xc) ),  x > xc  (q0 < 0)
{
    double  x1, x2, y1, y2, m1, m2;
    double  d, p, q, t;
    d = 1.0 - a * a;
    if( d < DBL_EPSILON )
        d = DBL_EPSILON;

// Set boundaries  ...(x1,y1)...(0,y0)...(x2,y2)...
// and slopes            m1    +ve,-ve      m2
    x1 = -gmv / a;
    y1 = x1 * (fpu - 0.5 * x1);
    m1 = fpu - x1;

    x2 = -gpv / a;
    y2 = x2 * (fmu - 0.5 * x2);
    m2 = fmu - x2;

// Apex necessarily at 0 
    t = m1 * x1 - (y1 - height);
    if( t <= 0.0 )                 // Lslope(0, height, -1/2., fpu, 0)
        p = fpu + sqrt(2.0 * height);
    else                           // Lslope(0,height,-d/2,fpu+a*gmv,gmv*gmv/2)
        p = (fpu + a * gmv) + sqrt(d * (2.0 * t + d * x1 * x1));

    t = m2 * x2 - (y2 - height);
    if( t <= 0.0 )                 // Rslope(0, height, -1/2., fmu, 0)
        q = fmu - sqrt(2.0 * height);
    else                           // Rslope(0,height,-d/2,fmu+a*gpv,gpv*gpv/2)
        q = (fmu + a * gpv) - sqrt(d * (2.0 * t + d * x2 * x2));

    *xc = 0.0;
    *p0 = p;
    *q0 = q;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  RanPos1
// 
// Purpose:   Draw sample in x >= 0 from normalised form of  p(x) exp(-L(x))
//
//                p(x)  =  s * delta(x)  +  1
//
//                L(x)  =  g * x  +  A * x * x / 2
//
//            where g = gradL(0)
//                  A = gradgradL
//
// History:   JS         29 Oct 1997, 30 Nov 1997, 3 Oct 2002
//-----------------------------------------------------------------------------
// 
double  RanPos1(      //   O  Sample value on exit
Rand_t   Rand,        // I O  Random generator state
double   s,           // I    Spike amplitude
double   g,           // I    Gaussian gradient at origin
double   A)           // I    Gaussian curvature >= 0
{
    double acc, mu, prob0, z;

    acc = sqrt(A);
    if( 0.3 * acc > g )                   // Reject from completed normal curve
    {
        mu = -g / A;
        prob0 = exp(g * mu / 2.0) * sqrt(A / TwoxPi);
        prob0 = 1.0 / (1.0 + s * prob0);
        do
        {
            if( Randouble(Rand) > prob0 )
                z = 0.0;
            else
                z = mu + Rangauss(Rand) / acc;
        } while( z < 0.0 );
    }
    else                                  // Reject from bounding exponential
    {
        prob0 = 1.0 / (1.0 + s * g);
        do
        {
            if( Randouble(Rand) > prob0 )
            {
                z = 0.0;
                break;
            }
            else
            {
                z = -log(Randouble(Rand)) / g;
            }
        } while( -log(Randouble(Rand)) < A * z * z / 2.0 );
    }
    return z;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  RanPos01
// 
// Purpose:   Draw sample in 0 < x < 1 from normalised form of
//
//                 exp( - g * x  -  A * x * x / 2 )
//
// History:   JS          9 Jun 1998
//-----------------------------------------------------------------------------
// 
double  RanPos01(     //   O  Sample value on exit
Rand_t   Rand,        // I O  Random generator state
double   g,           // I    Gaussian gradient at origin
double   A)           // I    Gaussian curvature >= 0
{
    double a, b, E, x, x0;
    int    flag = 0;

    if( g < -A / 2.0 )
    {
        g = -A - g;
        flag = 1;
    }
    if( A <= sqrt(DBL_EPSILON) )
    {
        if( g * g < DBL_EPSILON )
            x = Randouble(Rand);
        else
        {
            E = 1.0 - exp(-g);
            do x = - log(1.0 - E * Randouble(Rand)) / g;
            while( x >= 1.0 );
        }
    }
    else
    {
        A = sqrt(A);
        a = g / A;
        if( g > 0.0 )
        {
            x0 = 1.0 / (sqrt(1.0 + 0.25 * a * a) + 0.5 * a);
            if( x0 > A )
                x0 = A;
            b = x0 + a;
            E = 1.0 - exp(-A * b);
            if( E * E < DBL_EPSILON )
                do
                {
                    do x = A * Randouble(Rand);
                    while( log(Randouble(Rand)) > -0.5 * (x - x0) * (x - x0) );
                } while( x >= A );
            else
                do
                {
                    do x = - log(1.0 - E * Randouble(Rand)) / b;
                    while( log(Randouble(Rand)) > -0.5 * (x - x0) * (x - x0) );
                } while( x >= A );
        }
        else
        {
            b = a + A;
            x0 = (b < 1.0) ? b : 1.0;
            E = 1.0 - exp(-b * x0);
            if( E * E < DBL_EPSILON )
                do
                {
                    do x = b * Randouble(Rand);
                    while( log(Randouble(Rand)) > -0.5 * (x - x0) * (x - x0) );
                    x = (Ranint(Rand) < 0) ? - x - a : x - a;
                } while( x * (x - A) > 0.0 );
            else
                do
                {
                    do x = -log(1.0 - E * Randouble(Rand)) / x0;
                    while( log(Randouble(Rand)) > -0.5 * (x - x0) * (x - x0) );
                    x = (Ranint(Rand) < 0) ? - x - a : x - a;
                } while( x * (x - A) > 0.0 );
        }
        x /= A;
    }
    return flag ? 1.0 - x : x;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Ran1posX
// 
// Purpose:   Draw sample in x >= 0 from normalised form of
//
//                 x * exp( - g * x  -  A * x * x / 2 )
//
// History:   JS         14 Jan 2003
//-----------------------------------------------------------------------------
// 
double Ran1posX(      //   O  Sample value on exit
Rand_t   Rand,        // I O  Random generator state
double   g,           // I    Gaussian gradient at origin
double   A)           // I    Gaussian curvature >= 0
{
    double a, x, y, xcap, accept;

    if( A < 0.0 )
        return 0.0;
    if( A == 0.0 && g <= 0.0 )
        return 0.0;
    xcap = (g >= 0.0) ? 2.0 / (sqrt(g * g + 4.0 * A) + g)
                      : (sqrt(g * g + 4.0 * A) - g) / (2.0 * A);      // Apex
    a = 1.7000 / sqrt(A + 1.0 / (xcap * xcap));      // Adequate Cauchy width
    do                                         // Reject from bounding Cauchy
    {
        do
        {
            y = Rancauchy(Rand);
            x = xcap + y * a;
        } while( x <= 0.0 );
        accept = log((1.0 + y * y) * x / xcap)
                + (xcap - x) * (g + 0.5 * A * (xcap + x));
    } while( accept < log(Randouble(Rand)) ); 
    return  x;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:
//                           infinity   - g*x - A*x*x/2
//            logerf := log  INTEGRAL  e                dx
//                             x=0
//
//                    = log( sqrt(pi/2A) exp(g*g/2A) (1 - erf(g/sqrt(2A)) )
//
// History:   JS  25 Nov 1997     From "erf" code
//                31 Mar 1998, 9 Jun 1998
//                23 Jun 1998     Polish
//-----------------------------------------------------------------------------
// 
double logerf(   //   O Output value
double  g,       // I   Gradient at origin
double  A)       // I   Curvature at origin
{
    double  a, a0, a1, b0, b1, f, q, r, t, y;

    if( A <= 0.0 )
        return -log(g);

    a = sqrt(g * g / A);
    if( a < 2.0 )               // before which continued fraction too slow
    {
// Power series, proportional error = DBL_EPSILON * exp(a*a/2)
        t = a * a / 2.0;
        q = r = a;
        for( f = 1.5; q > r * DBL_EPSILON; f += 1.0 )
            r += q *= t / f;
        if( g >= 0.0 )
            y = log(SqrPi2 * exp(t) - r);
        else
            y = log(SqrPi2 * exp(t) + r);
    }
    else if( a < 8.6 )   // after which asymptotic series is fully accurate
    {
// Continued fraction, proportional error = DBL_EPSILON
        b0 = q = 0.0;
        a0 = b1 = f = 1.0;
        a1 = t = a * a / 2.0;
        for( f = 1.0; a0 / a1 - b0 / b1 > DBL_EPSILON; f = 1.0 / a1 )
        {
            q += 0.5;
            a0 = (a1 + a0 * q) * f;
            b0 = (b1 + b0 * q) * f;
            q += 0.5;
            r = q * f;
            a1 = t * a0 + r * a1;
            b1 = t * b0 + r * b1;
        }
        y = b1 * a / (2.0 * a1);
        if( g >= 0.0 )
            y = log(y);
        else
            y = log(2.0 * SqrPi2 * exp(t) - y);
    }
    else
    {
// Asymptotic series, proportional error = MAX(e(-a*a/2), DBL_EPSILON)
        if( g >= 0.0 )
        {
            t = 2.0 / (a * a);
            q = r = 1.0;
            for( f = 0.5 * t; fabs(q) > DBL_EPSILON; f += t )
            {
                r += q *= -f;
                if( f > 1.0 )   break;           // unnecessary if range OK
            }
            y = log((r - q * f / 2.0) / a);
        }
        else
            y = log(2.0 * SqrPi2) + a * a / 2.0;
    }
    return  y - log(A) / 2.0;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:
//                           inf, inf        - g.x - x.A.x/2
//            logerf2 := log INTEGRAL dx dy e                  (A non-negative)
//                           x=0, y=0
//
//            where  g.x = g1*x + g2*y,   x.A.x = A11*x*x + 2*A12*x*y + A22*y*y
//
// History:   JS   12 Jul 1998
//-----------------------------------------------------------------------------
// 
double logerf2(      //   O  Output value
double   g1,         // I    x gradient at origin
double   g2,         // I    y gradient at origin
double   A11,        // I    xx curvature >= 0
double   A12,        // I    xy curvature, A12*A12 <= A11*A22
double   A22)        // I    yy curvature >= 0
{
    double  ctheta;       // x coord of ivec =  y coord of jvec
    double  stheta;       // y coord of ivec = -x coord of jvec
    double  val_i;        // Eigenvalue for ivec
    double  val_j;        // Eigenvalue for jvec
    double  origi, origj; // Origin relative to eigencentre
    double  xveci, xvecj; // x-axis relative to eigencoords
    double  yveci, yvecj; // y-axis relative to eigencoords
    double  u, v;         // Eigencoords
    double  a;            // Radius to wedge corner in eigencoords
    double  s, t;         // local
    double  TINY = pow(DBL_EPSILON, 2.0/3.0);

// Separable if  |A12| * <x(separable)> * <y(separable)>  <=  sqrt(EPSILON)
    u = (g1 > 0.0) ? 1.0 / (A11 + g1 * g1) : (A11 + g1 * g1) / (A11 * A11);
    v = (g2 > 0.0) ? 1.0 / (A22 + g2 * g2) : (A22 + g2 * g2) / (A22 * A22);
    if( A12 * A12 * u * v <= DBL_EPSILON )
        return  logerf(g1, A11) + logerf(g2, A22);

// Generic case ....
// Scale Gaussian to standard form exp(-u*u/2-v*v/2)
    if( A11 >= A22 )
    {
        t = (A11 - A22) / 2.0;
        t += sqrt(t * t + A12 * A12);
        s = sqrt(t * t + A12 * A12);
        val_i = t + A22;
        ctheta = t / s;
        stheta = A12 / s;
    }
    else
    {
        t = (A22 - A11) / 2.0;
        t += sqrt(t * t + A12 * A12);
        s = sqrt(t * t + A12 * A12);
        val_i = t + A11;
        ctheta = A12 / s;
        stheta = t / s;
    }
    val_j = (A11 * A22 - A12 * A12) / val_i;
    if( val_j < TINY * sqrt(A11 * A22) )
        val_j = TINY * sqrt(A11 * A22);      // Protect singular eval
    val_i = sqrt(val_i);
    val_j = sqrt(val_j);
    origi = ( g1 * ctheta + g2 * stheta) / val_i;
    origj = (-g1 * stheta + g2 * ctheta) / val_j;
    xveci =  val_i * ctheta;
    xvecj = -val_j * stheta;
    yveci =  val_i * stheta;
    yvecj =  val_j * ctheta;
    
    a = sqrt(origi * origi + origj * origj);
    if( a > 0.0 )
    {
        origi /= a;
        origj /= a;
    }
    s = eigenerf2(a,
               origi * xveci + origj * xvecj, -origj * xveci + origi * xvecj,
               origi * yveci + origj * yvecj, -origj * yveci + origi * yvecj);
    s = s + a * a / 2.0 - log(val_i * val_j);
    return s;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:
//                                        -r*r/2
//            eigenerf2 :=  log INTEGRAL e       dArea
//                               wedge
//                                              yyy..
//                                           yyy.....
//                                        yyy........
//                                     yyy...........
//                         |        yyy..............
//                       --O-->   a..................
//                         |        xxxxx............
//                                       xxxxx.......
//                                            xxxxx..
//
//            0 < Angle(xxxx, yyyy) < Pi
//
// History:   JS  12 Jul 1998
//-----------------------------------------------------------------------------
// 
static double eigenerf2(// Output value
double  a,          // I   Radius to corner of wedge
double  xcos,       // I   x-edge direction is (xcos,xsin) relative to radius a
double  xsin,       // I   x-edge direction is (xcos,xsin) relative to radius a
double  ycos,       // I   y-edge direction is (ycos,ysin) relative to radius a
double  ysin)       // I   y-edge direction is (ycos,ysin) relative to radius a
{
#undef  PLUS
#define PLUS(x,y)  ((x>y) ? x+log(1.0+exp(y-x)) : y+log(1.0+exp(x-y)))
#undef  MINUS
#define MINUS(x,y) ((x>y) ? x+log(1.0-exp(y-x)) : log(DBL_MIN))
#undef  erfint
#define erfint(x)  (logerf((x), 1.0) - (x) * (x) / 2.0)

    double t, u, v;

// Normalisation to unit vectors
    t = sqrt(xcos * xcos + xsin * xsin);
    xcos /= t;
    xsin /= t;
    t = sqrt(ycos * ycos + ysin * ysin);
    ycos /= t;
    ysin /= t;
// Force upward orientation
    if( xsin + ysin < 0.0 )
    {
        t = xcos;    xcos = ycos;    ycos = t;
        t = xsin;    xsin = -ysin;   ysin = -t;
    }

// Select case
    if( ycos > 0.5 )
    {                                                             // case P
            u = wedge(a, xsin, ysin);
    }
    else if( xcos < -0.5 )
    {                                                             // case R
            u = wedge(a, -xsin, -ysin);
            t = LogSqr2Pi + erfint2(a * ysin, a * xsin);
            u = PLUS(u, t);
    }
    else if( ycos + 1.732 * ysin < 0.0 )
    {
        if( ysin < -0.866 )
        {                                                         // case W
            v = wedge(a, xcos, -ycos);
            t = erfint(a * xsin) + erfint(a * xcos);
            v = PLUS(v, t);
            t = erfint(-a * ysin) + erfint(a * ycos);
            v = PLUS(v, t);
            u = 2.0 * LogSqr2Pi;
            u = MINUS(u, v);
        }
        else if( 1.732 * xsin < xcos )
        {                                                         // case V
            u = LogSqr2Pi + erfint(a * ysin);
            t = wedge(a, -xsin, ysin);
            u = MINUS(u, t);
        }
        else
        {                                                         // case T
            u = LogSqr2Pi + erfint(a * ysin);
            v = erfint(a * xsin) + erfint(a * xcos);
            if( -ysin > -xcos )
            {
                t = wedge(a, -xcos, -ysin);
                u = PLUS(u, t);
            }
            else
            {
                t = wedge(a, -ysin, -xcos);
                v = PLUS(v, t);
            }
            u = MINUS(u, v);
        }
    }
    else
    {
        if( xsin < -0.866 )
        {                                                         // case U
            u = erfint(a * ysin) + erfint(a * ycos);
            t = erfint(-a * xsin) + erfint(a * xcos);
            u = PLUS(u, t);
            t = wedge(a, -xcos, ycos);
            u = MINUS(u, t);
        }
        else if( 1.732 * xsin > xcos )
        {                                                         // case Q
            u = wedge(a, -xcos, -ycos);
            t = erfint(a * ysin) + erfint2(a * ycos, a * xcos);
            u = PLUS(u, t);
            if( ysin > xsin )
            {
                t = erfint(a * xcos) + erfint2(a * xsin, a * ysin);
                u = MINUS(u, t);
            }
            else
            {
                t = erfint(a * xcos) + erfint2(a * ysin, a * xsin);
                u = PLUS(u, t);
            }
        }
        else
        {                                                         // case S
            u = erfint(a * ysin) + erfint(a * ycos);
            if( -ycos >= xsin )
            {
                t = wedge(a, xsin, -ycos);
                u = PLUS(u, t);
            }
            else
            {
                t = wedge(a, -ycos, xsin);
                u = MINUS(u, t);
            }
        }
    }
    return u;
#undef erfint
#undef MINUS
#undef PLUS
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:
//                              x=b     -x*x/2
//            erfint2 :=  log INTEGRAL e      dx,     a < b
//                              x=a
//
// History:   JS  12 Jul 1998
//-----------------------------------------------------------------------------
// 
static double erfint2(  //   O Output value
       double  a,       // I   Lower limit
       double  b)       // I   Upper limit
{
    double  aa, bb, e1, e2, ee, s, t;

    aa = fabs(a);
    bb = fabs(b);
    if( aa < bb )
    {
        s = aa;    aa = bb;    bb = s;
    }
    e1 = DBL_EPSILON / (aa + 1.0); 
    t = b - a;
    ee = (aa + fabs(t)) * t;
    e2 = ee * ee * fabs(t);
    if( e2 < e1 )
        return log((1.0 - (a + t / 3.0) * t / 2.0) * t) - a * a / 2.0;

    s = logerf(aa, 1.0) - aa * aa / 2.0;
    t = logerf(bb, 1.0) - bb * bb / 2.0;
    if( a * b >= 0.0 )
        return t + log(1.0 - exp(s-t));
    else
        return log(Sqr2Pi - exp(s) - exp(t));
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:
//
//            wedge  := log(w(a,t) - w(a,s))    (a >= 0, -1 < s < t < 1)  where
// 
//                      infinity
//            w(a,x) := INTEGRAL (arcsin(x) - arcsin(a*x/r)) exp(-r*r/2) r dr
//                        r=a
//
// Notes:     Series abs cgt  for all |s| < 1, |t| < 1,   but cgce slow
//            and more storage would be needed if |s| ~ 1 or |t| ~ 1 .
//            This version allows values up to sqrt(3)/2, as from eigenerf2.
//
// History:   JS  12 Jul 1998, 4 Sep 1999, 12 Sep 1999
//-----------------------------------------------------------------------------
// 
static double wedge(    //   O Output value
       double  a,       // I   Input radius
       double  s,       // I   Lower angle (as sine)
       double  t)       // I   Upper angle (as sine)
{
#define STORE 120           // Want s^STORE < DBL_EPSILON for full accuracy
    double w[STORE];        // arcsin_terms(s,t)  =  SUM w[k] (a/r)^(k+k+1)
    double p;    // p[k] = INT[a,inf] (1 - (a/r)^(k+k+1)) exp(a*a/2-r*r/2) r dr
    double q;    // max p[k] when recurrence relation for p started internally
    double b, m, r, u, v, sum;
    int    j, k, n;

// Series {x, (1/2)x^3/3, (1*3/2*4)x^5/5, (1*3*5/2*4*6)x^7/7),..} difference
    n = 0;
    m = 1.0;
    if( s * t > 0.0 )
    {                      // Implicit difference between s and t
        r = u = v = 1.0;
        while( (w[n++] = r * u * (t - s) / m) > DBL_EPSILON )
        {
            if( n == STORE )
                break;                  // Memory overflow trap
            r *= 1.0 - 0.5 / n;
            v *= t;
            u = u * s + v;
            v *= t;
            u = u * s + v;
            m += 2.0;
        }
    }
    else
    {                      // Explicit difference between s and t
        r = 1.0;
        u = s;
        v = t;
        while( (w[n++] = r * (v - u) / m) > DBL_EPSILON )
        {
            if( n == STORE )
                break;                  // Memory overflow trap
            r *= 1.0 - 0.5 / n;
            u *= s * s;
            v *= t * t;
            m += 2.0;
        }
    }

// Initialise recurrence relation for exponential parts
    b = a * a / 2.0;
    sum = 0.0;

    if( a < 7.0 )      // Balance point for error control
    {
// Start at beginning unless errors would amplify too far
        p = 1.0 - a * exp(logerf(a, 1.0));
        for( k = 0; k < n; ++k )
        {
            sum += w[k] * p;
            p = 1.0 - b * p / (k + 0.5);
        }
    }
    else
    {
// Start in middle at maximum term if errors would otherwise amplify
        if( b < STORE )      // (= dimension limit)
        {
            k = (int)b;
            u = k * k;
            p = 0.5 + (0.031 - 0.4 / u) / u;  // Phenomenological asymptote
            j = k + 1;
            u = j * j;
            q = 0.5 + (0.031 - 0.4 / u) / u;  // Phenomenological asymptote
            q = (j - 0.5) * (1.0 - q) / j;
            u = b - k;
            q = u * q + (1.0 - u) * p;        // Interpolated maximum term
        }
        else
        {
            j = STORE;
            q = (j + 0.5) / (j + b);          // Phenomenological approximant
        }
// Count down from maximum term
        p = q;
        for( k = j - 1; k >= n; --k )
            p = (k - 0.5) * (1.0 - p) / b;
        for( ; k >= 0; --k )
        {
            sum += w[k] * p;
            p = (k - 0.5) * (1.0 - p) / b;
        }
// Count up from maximum term
        p = q;
        for( k = j; k < n; ++k )
        {
            p = 1.0 - b * p / (k - 0.5);
            sum += w[k] * p;
        }
    }
    return  log(sum) - b;
#undef STORE
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  logfactorial
// 
// Purpose:   log( k! ) = lookup table up to 99, then logGamma(k+1)
//
// History:   JS         15 May 1998, 24 Mar 2001, 3 Oct 2002, 18 Aug 2003
//-----------------------------------------------------------------------------
// 
double logfactorial(      //   O  Value = log(k!)
unsigned   k)             // I    Argument
{
static const double  logfact[100] = {  0.0000000000000000,
            0.0000000000000000,  0.6931471805599453,  1.7917594692280550,
            3.1780538303479456,  4.7874917427820460,  6.5792512120101010,
            8.5251613610654143, 10.6046029027452502, 12.8018274800814696,
           15.1044125730755153, 17.5023078458738858, 19.9872144956618861,
           22.5521638531234229, 25.1912211827386815, 27.8992713838408916,
           30.6718601060806728, 33.5050734501368889, 36.3954452080330536,
           39.3398841871994940, 42.3356164607534850, 45.3801388984769080,
           48.4711813518352239, 51.6066755677643736, 54.7847293981123192,
           58.0036052229805199, 61.2617017610020020, 64.5575386270063311,
           67.8897431371815350, 71.2570389671680090, 74.6582363488301644,
           78.0922235533153106, 81.5579594561150372, 85.0544670175815174,
           88.5808275421976788, 92.1361756036870925, 95.7196945421432025,
           99.3306124547874269,102.9681986145138126,106.6317602606434591,
          110.3206397147573954,114.0342117814617032,117.7718813997450715,
          121.5330815154386340,125.3172711493568951,129.1239336391272149,
          132.9525750356163099,136.8027226373263685,140.6739236482342594,
          144.5657439463448860,148.4777669517730321,152.4095925844973578,
          156.3608363030787852,160.3311282166309070,164.3201122631951814,
          168.3274454484276523,172.3527971391628016,176.3958484069973517,
          180.4562914175437711,184.5338288614494905,188.6281734236715912,
          192.7390472878449024,196.8661816728899940,201.0093163992815267,
          205.1681994826411985,209.3425867525368356,213.5322414945632612,
          217.7369341139542273,221.9564418191303340,226.1905483237275933,
          230.4390435657769523,234.7017234428182677,238.9783895618343230,
          243.2688490029827142,247.5729140961868839,251.8904022097231944,
          256.2211355500095255,260.5649409718632093,264.9216497985528010,
          269.2910976510198225,273.6731242856937041,278.0675734403661429,
          282.4742926876303960,286.8931332954269940,291.3239500942703076,
          295.7666013507606240,300.2209486470141318,304.6868567656687155,
          309.1641935801469219,313.6528299498790618,318.1526396202093268,
          322.6634991267261769,327.1852877037752172,331.7178871969284731,
          336.2611819791984770,340.8150588707990179,345.3794070622668541,
          349.9541180407702369,354.5390855194408088,359.1342053695753988};
    double s, y;

    if( k < 100 )
        return logfact[k];
    y = k + 1.0;
    s = y + 2.269488974204959960;
    s = y + 1.517473649153287398 / s;
    s = y + 1.011523068126841711 / s;
    s = y + 0.525606469002695417 / s;
    s = y + 0.252380952380952380 / s;
    s = y + 0.033333333333333333 / s;
    s =     0.083333333333333333 / s;
    s = s + 0.91893853320467 - y + (y - 0.5) * log(y);
    return s;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  logGamma
// 
// Purpose:   log( Gamma(x) )
//
// History:   JS        18 Aug 2003  From gammln 8 Dec 1994
//                                   Continued fraction (Abramowitz and Stegun)
//-----------------------------------------------------------------------------
// 
double logGamma(            //   O  Value = log(Gamma(x))
double   x)                 // I    Argument
{
    int    k;
    double s, t, y;

    y = x;
    for( k = 16; y < 16.0; --k )
    {
        y += 1.0;
    }
// Continued fraction for s = log(Gamma(y))
    s = y + 2.269488974204959960;
    s = y + 1.517473649153287398 / s;
    s = y + 1.011523068126841711 / s;
    s = y + 0.525606469002695417 / s;
    s = y + 0.252380952380952380 / s;
    s = y + 0.033333333333333333 / s;
    s =     0.083333333333333333 / s;
    s = s + 0.91893853320467 - y + (y - 0.5) * log(y);
// Return to s = log(Gamma(original x))
    t = 1.0;
    for( ++k; k <= 16; ++k )
    {
        y -= 1.0;
        t *= y;
    }
    s -= log(t);
// Required output
    return  s;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  InvNorm
// 
// Purpose:   Inverse of normal distribution, returning x(y) where
//                 y(x) = Normal(x)
//                           x
//                      = INTEGRAL exp(-u*u/2) du / sqrt(2*pi)
//                         -inf
//                      = ( 1 + erf(x/sqrt(2)) ) / 2
//
//            Rational approximations claimed to be correct to 1.15e-9
//
// History:   JS        13 Sep 2003
//            From Peter J Acklam www//home.online.no/~pjacklam/notes/invnorm
//-----------------------------------------------------------------------------
// 
double InvNorm(   //   O Inverse normal (= number of standard deviations)
double  y)        // I   Argument       (= cumulative probability)
{
    double  r, s;
    if( y < 0.02425 )
    {
        s = sqrt(-2.0 * log(y));
        return  + (2.938163982698783  + 
                  (4.374664141464968  -
                  (2.549732539343734  +
                  (2.400758277161838  +
                  (0.3223964580411365 +
                   0.007784894002430293 * s) * s) * s) * s) * s)
                / (1.0                +
                  (3.754408661907416  +
                  (2.445134137142996  +
                  (0.3224671290700398 +
                   0.007784695709041462 * s) * s) * s) * s);
    }
    else if( y > 0.97575 )
    {
        s = sqrt(-2.0 * log(1.0 - y));
        return  - (2.938163982698783  + 
                  (4.374664141464968  -
                  (2.549732539343734  +
                  (2.400758277161838  +
                  (0.3223964580411365 +
                   0.007784894002430293 * s) * s) * s) * s) * s)
                / (1.0                +
                  (3.754408661907416  +
                  (2.445134137142996  +
                  (0.3224671290700398 +
                   0.007784695709041462 * s) * s) * s) * s);
    }
    else
    {
        s = y - 0.5;
        r = s * s;
        return  s * (2.506628277459239 -
                    (30.66479806614716 -
                    (138.3577518672690 -
                    (275.9285104469687 -
                    (220.9460984245205 -
                     39.69683028665376 * r) * r) * r) * r) * r)
                  / (1.0 -
                    (13.28068155288572 -
                    (66.80131188771972 -
                    (155.6989798598866 -
                    (161.5858368580409 -
                     54.47609879822406 * r) * r) * r) * r) * r);
    }
}


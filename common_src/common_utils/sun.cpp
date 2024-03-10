#include "common_utils/sun.h"
#include <math.h>

    #define PI          3.14159265358979323846
    #define RADEG       (180.0/PI)
    #define DEGRAD      (PI/180.0)
    #define sind(x)     sin((x)*DEGRAD)
    #define cosd(x)     cos((x)*DEGRAD)
    #define tand(x)     tan((x)*DEGRAD)
    #define asind(x)    (RADEG*asin(x))
    #define acosd(x)    (RADEG*acos(x))
    #define atand(x)    (RADEG*atan(x))
    #define atan2d(y,x) (RADEG*atan2((y),(x)))


    double rev( double x )
    {
        return  x - floor(x/360.0)*360.0;
    }


    double cbrt( double x )
    {
        if ( x > 0.0 )
            return exp( log(x) / 3.0 );
        else if ( x < 0.0 )
            return -cbrt(-x);
        else /* x == 0.0 */
            return 0.0;
    }
float3 Sun::sun_direction(EnvironmentParameters &params)
{
    //taken from here http://stjarnhimlen.se/comp/tutorial.html

    double d = 367*params.year - (7*(params.year + ((params.month+9)/12)))/4 + (275*params.month)/9 + params.day - 730530;
    
    double w = 282.9404 + 4.70935E-5   * d;    //(longitude of perihelion deg)
    double a = 1.000000;                               //(mean distance, a.u.)
    double e = 0.016709 - 1.151E-9             * d;    //(eccentricity)
    double M = 356.0470 + 0.9856002585 * d;    //(mean anomaly deg)
    double oblecl = 23.4393 - 3.563E-7 * d; //obliquity of the ecliptic deg
    double L = w + M;
    double E = M + (180/PI) * e * sind(M) * (1 + e * cosd(M));
    double x = cosd(E) - e;
    double y = sind(E) * sqrt(1 - e*e);
    double r = sqrt(x*x + y*y);
    double v = atan2d( y, x );
    double lon = v + w;
    
    x = r*cosd(lon);
    y = r*sind(lon);
    double z = 0;

    double xequat = x;
    double yequat = y * cosd(oblecl) - z * sind(oblecl);
    double zequat = y * sind(oblecl) + z * cosd(oblecl);

    r    =  sqrt( xequat*xequat + yequat*yequat + zequat*zequat );
    double RA   =  atan2d( yequat, xequat );
    double Decl =  atan2d( zequat, sqrt( xequat*xequat + yequat*yequat) );
    double RA_hours = RA/15.0;

    double UT = params.hours + params.minutes/60.0 + params.seconds/3600.0;
    double GMST0 = rev(L + 180)/15;
    double SIDTIME = GMST0 + UT + params.longitude_deg/15;
    if (SIDTIME < 0)
        SIDTIME += 24;
    else if (SIDTIME > 24)
        SIDTIME -= 24;
    
    double HA = (SIDTIME + RA_hours)*15; // in degrees

    x = cosd(HA)*cosd(Decl);
    y = sind(HA)*cosd(Decl);
    z = sind(Decl);

    double lat = params.latitude_deg;
    double xhor = x * sind(lat) - z * cosd(lat);
    double yhor = y;
    double zhor = x * cosd(lat) + z * sind(lat);
    
    return normalize(float3(xhor,yhor,zhor));
}

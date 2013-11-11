/* 
 * Partial copyright by Doug Mink, Harvard-Smithsonian Center for Astrophysics
 * and by Patrick Wallace
 */

#include <math.h>
#include <stdio.h>	/* for fprintf() and sprintf() */
#include <ctype.h>
#include <string.h>
#include "common.h"

#define raddeg(X) (Q3C_RADEG*X)
#define degrad(X) (Q3C_DEGRA*X)

static void v2s3(q3c_coord_t [], q3c_coord_t *, q3c_coord_t *, q3c_coord_t *);
static void s2v3 (q3c_coord_t,q3c_coord_t , q3c_coord_t ,q3c_coord_t [] );
static char *eqstrn (q3c_coord_t, q3c_coord_t);


/* Convert geocentric equatorial rectangular coordinates to
   right ascension, declination, and distance */

static void v2s3(	q3c_coord_t pos[3],	/* x,y,z geocentric equatorial position of object */
			q3c_coord_t *rra,	/* Right ascension in radians (returned) */
			q3c_coord_t *rdec,	/* Declination in radians (returned) */
			q3c_coord_t *r	/* Distance to object in same units as pos (returned) */
		)
{
    q3c_coord_t x,y,z,rxy,rxy2,z2;

    x = pos[0];
    y = pos[1];
    z = pos[2];

    *rra = q3c_atan2 (y, x);

    /* Keep RA within 0 to 2pi range */
    if (*rra < 0.0)
	*rra = *rra + (2.0 * Q3C_PI);
    if (*rra > 2.0 * Q3C_PI)
	*rra = *rra - (2.0 * Q3C_PI);

    rxy2 = x*x + y*y;
    rxy = q3c_sqrt (rxy2);
    *rdec = q3c_atan2 (z, rxy);

    z2 = z * z;
    *r = q3c_sqrt (rxy2 + z2);

    return;
}

static void s2v3 (q3c_coord_t rra,	/* Right ascension in radians */
	q3c_coord_t rdec,	/* Declination in radians */
	q3c_coord_t r,	/* Distance to object in same units as pos */
	q3c_coord_t pos[3]	/* x,y,z geocentric equatorial position of object (returned) */
	)
{
    pos[0] = r * cos (rra) * cos (rdec);
    pos[1] = r * sin (rra) * cos (rdec);
    pos[2] = r * sin (rdec);

    return;
}

/* Return string with right ascension in hours and declination in degrees */

static char *eqstrn (
	q3c_coord_t	dra,		/* Right ascension in degrees */
	q3c_coord_t	ddec		/* Declination in degrees */
	)

{
	char	*eqcoor;	/* ASCII character string of position (returned) */
	char	decp;
	int	rah,irm,decd,decm;
	q3c_coord_t	xpos,ypos,xp,yp,ras,decs;

    /*  Right ascension to hours, minutes, and seconds */
    xpos = dra / 15.0;
    rah = (int) xpos;
    xp = (q3c_coord_t) 60.0 * (xpos - (q3c_coord_t) rah);
    irm = (int) xp;
    ras = (q3c_coord_t) 60.0 * (xp - (q3c_coord_t) irm);

    /* Declination to degrees, minutes, seconds */
    if (ddec < 0) {
	ypos = -ddec;
	decp = '-';
	}
    else {
	decp = '+';
	ypos = ddec;
	}
    decd = (int) ypos;
    yp = (q3c_coord_t) 60.0 * (ypos - (q3c_coord_t) decd);
    decm = (int) yp;
    decs = (q3c_coord_t) 60.0 * (yp - (q3c_coord_t) decm);

    eqcoor = malloc (32);
    (void)sprintf (eqcoor,"%02d:%02d:%06.3f %c%02d:%02d:%05.2f",
		   rah,irm,ras,decp,decd,decm,decs);
    if (eqcoor[6] == ' ')
	eqcoor[6] = '0';
    if (eqcoor[20] == ' ')
	eqcoor[20] = '0';

    return (eqcoor);
}



int	idg=0;

/*  l2,b2 system of galactic coordinates
 *  p = 192.25       ra of galactic north pole (mean b1950.0)
 *  q =  62.6        inclination of galactic to mean b1950.0 equator
 *  r =  33          longitude of ascending node
 *  p,q,r are degrees

 *  Equatorial to galactic rotation matrix
    (The Eulerian angles are p, q, 90-r)
	+cp.cq.sr-sp.cr	+sp.cq.sr+cp.cr	-sq.sr
	-cp.cq.cr-sp.sr	-sp.cq.cr+cp.sr	+sq.cr
	cp.sq		+sp.sq		+cq
 */


/*  l2,b2 system of galactic coordinates
    p = 192.25       ra of galactic north pole (mean b1950.0)
    q =  62.6        inclination of galactic to mean b1950.0 equator
    r =  33          longitude of ascending node
    p,q,r are degrees */

/*  Equatorial to galactic rotation matrix
    The eulerian angles are p, q, 90-r
	+cp.cq.sr-sp.cr     +sp.cq.sr+cp.cr     -sq.sr
	-cp.cq.cr-sp.sr     -sp.cq.cr+cp.sr     +sq.cr
	+cp.sq              +sp.sq              +cq		*/

static q3c_coord_t jgal[3][3] =
	{{-0.054875539726,-0.873437108010,-0.483834985808},
	{0.494109453312,-0.444829589425, 0.746982251810},
	{-0.867666135858,-0.198076386122, 0.455983795705}};

/* Transform J2000 equatorial coordinates to IAU 1958 galactic coordinates */

void fk52gal (q3c_coord_t *dtheta, q3c_coord_t *dphi)
{
	/* J2000 right ascension in degrees
		   Galactic longitude (l2) in degrees (returned) */
	/* J2000 declination in degrees
		Galactic latitude (b2) in degrees (returned) */
	/* Rotation matrices by P.T.Wallace, Starlink eqgal and galeq, March 1986 */
	/*  Input equatorial coordinates are J2000 FK5.
    Use gal2fk4() if converting from B1950 FK4 coordinates.
    Reference: Blaauw et al, MNRAS,121,123 (1960) */

    q3c_coord_t pos[3],pos1[3],r,dl,db,rl,rb,rra,rdec,dra,ddec;
    char *eqcoor;
    int i;

    /*  Spherical to cartesian */
    dra = *dtheta;
    ddec = *dphi;
    rra = degrad (dra);
    rdec = degrad (ddec);
    r = 1.0;
    s2v3 (rra,rdec,r,pos);

    /*  Rotate to galactic */
    for (i = 0; i < 3; i++) {
		pos1[i] = pos[0]*jgal[i][0] + pos[1]*jgal[i][1] + pos[2]*jgal[i][2];
	}

    /*  Cartesian to spherical */
    v2s3 (pos1,&rl,&rb,&r);

    dl = raddeg (rl);
    db = raddeg (rb);
    *dtheta = dl;
    *dphi = db;

    /*  Print result if in diagnostic mode */
    if (idg) {
		eqcoor = eqstrn (dra,ddec);
		fprintf (stderr,"FK52GAL: J2000 RA,Dec= %s\n",eqcoor);
		fprintf (stderr,"FK52GAL: long = %.5f lat = %.5f\n",dl,db);
		free (eqcoor);
	}

    return;
}


/*--- Transform IAU 1958 galactic coordinates to J2000 equatorial coordinates */

void gal2fk5 (q3c_coord_t *dtheta, q3c_coord_t *dphi)
{
	/* Galactic longitude (l2) in degrees
		   J2000.0 ra in degrees (returned) */
	/* Galactic latitude (b2) in degrees
		   J2000.0 dec in degrees (returned) */
	/*  Output equatorial coordinates are J2000.
	Use gal2fk4() to convert to B1950 coordinates.
    Reference: Blaauw et al, MNRAS,121,123 (1960) */

    q3c_coord_t pos[3],pos1[3],r,dl,db,rl,rb,rra,rdec,dra,ddec;
    int i;
    char *eqcoor;

    /*  Spherical to Cartesian */
    dl = *dtheta;
    db = *dphi;
    rl = degrad (dl);
    rb = degrad (db);
    r = 1.0;
    s2v3 (rl,rb,r,pos);

    /*  Rotate to equatorial coordinates */
    for (i = 0; i < 3; i++) {
	    pos1[i] = pos[0]*jgal[0][i] + pos[1]*jgal[1][i] + pos[2]*jgal[2][i];
	    }

    /*  Cartesian to Spherical */
    v2s3 (pos1,&rra,&rdec,&r);
    dra = raddeg (rra);
    ddec = raddeg (rdec);
    *dtheta = dra;
    *dphi = ddec;

    /*  Print result if in diagnostic mode */
    if (idg) {
		fprintf (stderr,"GAL2FK5: long = %.5f lat = %.5f\n",dl,db);
		eqcoor = eqstrn (dra,ddec);
		fprintf (stderr,"GAL2FK5: J2000 RA,Dec= %s\n",eqcoor);
		free (eqcoor);
	}

    return;
}



#include <math.h>
#include <stdio.h>
#include <stddef.h>

/* Subroutine */ int ranlux_0_(int n__, float *rvec, int *lenv, int *
							   isdext, int *lout, int *inout, int *k1, int *k2, 
							   int *lux, int *ins)
{
    /* Initialized data */
	
    static int notyet = 1;
    static int luxlev = 3;
    static int in24 = 0;
    static int kount = 0;
    static int mkount = 0;
    static int i24 = 24;
    static int j24 = 10;
    static float carry = 0.f;
    static int ndskip[5] = { 0,24,73,199,365 };
	
    /* System generated locals */
    int i__1, i__2;
	
    /* Local variables */
    static int i__, k, lp, isd, isk;
    static float uni;
    static int ilx, ivec, izip, next[24], izip2, jseed;
    static float seeds[24];
    static int inner, nskip;
    static float twom12, twom24;
    static int inseed, iseeds[24], iouter;
	
	
	/*         Subtract-and-borrow random number generator proposed by */
	/*         Marsaglia and Zaman, implemented by F. James with the name */
	/*         RCARRY in 1991, and later improved by Martin Luescher */
	/*         in 1993 to produce "Luxury Pseudorandom Numbers". */
	/*     Fortran 77 coded by F. James, 1993 */
	
	/*       references: */
	/*  M. Luscher, Computer Physics Communications  79 (1994) 100 */
	/*  F. James, Computer Physics Communications 79 (1994) 111 */
	
	/*   LUXURY LEVELS. */
	/*   ------ ------      The available luxury levels are: */
	
	/*  level 0  (p=24): equivalent to the original RCARRY of Marsaglia */
	/*           and Zaman, very long period, but fails many tests. */
	/*  level 1  (p=48): considerable improvement in quality over level 0, */
	/*           now passes the gap test, but still fails spectral test. */
	/*  level 2  (p=97): passes all known tests, but theoretically still */
	/*           defective. */
	/*  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible */
	/*           correlations have very small chance of being observed. */
	/*  level 4  (p=389): highest possible luxury, all 24 bits chaotic. */
	
	/* !!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
	/* !!!  Calling sequences for RANLUX:                                  ++ */
	/* !!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++ */
	/* !!!                   32-bit random floating point numbers between  ++ */
	/* !!!                   zero (not included) and one (also not incl.). ++ */
	/* !!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++ */
	/* !!!               one 32-bit int INT and sets Luxury Level LUX  ++ */
	/* !!!               which is int between zero and MAXLEV, or if   ++ */
	/* !!!               LUX .GT. 24, it sets p=LUX directly.  K1 and K2   ++ */
	/* !!!               should be set to zero unless restarting at a break++ */
	/* !!!               point given by output of RLUXAT (see RLUXAT).     ++ */
	/* !!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four ints++ */
	/* !!!               which can be used to restart the RANLUX generator ++ */
	/* !!!               at the current point by calling RLUXGO.  K1 and K2++ */
	/* !!!               specify how many numbers were generated since the ++ */
	/* !!!               initialization with LUX and INT.  The restarting  ++ */
	/* !!!               skips over  K1+K2*E9   numbers, so it can be long.++ */
	/* !!!   A more efficient but less convenient way of restarting is by: ++ */
	/* !!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++ */
	/* !!!                   ISVEC of 25 32-bit ints (see RLUXUT)      ++ */
	/* !!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++ */
	/* !!!                 32-bit int seeds, to be used for restarting ++ */
	/* !!!      ISVEC must be dimensioned 25 in the calling program        ++ */
	/* !!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
    /* Parameter adjustments */
    if (rvec) {
		--rvec;
	}
    if (isdext) {
		--isdext;
	}
	
    /* Function Body */
    switch(n__) {
		case 1: goto L_rluxin;
		case 2: goto L_rluxut;
		case 3: goto L_rluxat;
		case 4: goto L_rluxgo;
	}
	
	/*                               default */
	/*  Luxury Level   0     1     2   *3*    4 */
	/* orresponds to p=24    48    97   223   389 */
	/*     time factor 1     2     3     6    10   on slow workstation */
	/*                 1    1.5    2     3     5   on fast mainframe */
	
	/*  NOTYET is .TRUE. if no initialization has been performed yet. */
	/*              Default Initialization by Multiplicative Congruential */
    if (notyet) {
		notyet = 0;
		jseed = 314159265;
		inseed = jseed;
		printf("Ranlux default initialization %d \n",jseed);
		luxlev = 3;
		nskip = ndskip[luxlev];
		lp = nskip + 24;
		in24 = 0;
		kount = 0;
		mkount = 0;
		printf("Ranlux default luxury level = %d, p = %d \n", luxlev, lp);
		twom24 = 1.f;
		for (i__ = 1; i__ <= 24; ++i__) {
			twom24 *= .5f;
			k = jseed / 53668;
			jseed = (jseed - k * 53668) * 40014 - k * 12211;
			if (jseed < 0) {
				jseed += 2147483563;
			}
			iseeds[i__ - 1] = jseed % 16777216;
			/* L25: */
		}
		twom12 = twom24 * 4096.f;
		for (i__ = 1; i__ <= 24; ++i__) {
			seeds[i__ - 1] = (float) iseeds[i__ - 1] * twom24;
			next[i__ - 1] = i__ - 1;
			/* L50: */
		}
		next[0] = 24;
		i24 = 24;
		j24 = 10;
		carry = 0.f;
		if (seeds[23] == 0.f) {
			carry = twom24;
		}
    }
	
	/*          The Generator proper: "Subtract-with-borrow", */
	/*          as proposed by Marsaglia and Zaman, */
	/*          Florida State University, March, 1989 */
	
    i__1 = *lenv;
    for (ivec = 1; ivec <= i__1; ++ivec) {
		uni = seeds[j24 - 1] - seeds[i24 - 1] - carry;
		if (uni < 0.f) {
			uni += 1.f;
			carry = twom24;
		} else {
			carry = 0.f;
		}
		seeds[i24 - 1] = uni;
		i24 = next[i24 - 1];
		j24 = next[j24 - 1];
		rvec[ivec] = uni;
		/*  small numbers (with less than 12 "significant" bits) are "padded". */
		if (uni < twom12) {
			rvec[ivec] += twom24 * seeds[j24 - 1];
			/*        and zero is forbidden in case someone takes a logarithm */
			if (rvec[ivec] == 0.f) {
				rvec[ivec] = twom24 * twom24;
			}
		}
		/*        Skipping to luxury.  As proposed by Martin Luscher. */
		++in24;
		if (in24 == 24) {
			in24 = 0;
			kount += nskip;
			i__2 = nskip;
			for (isk = 1; isk <= i__2; ++isk) {
				uni = seeds[j24 - 1] - seeds[i24 - 1] - carry;
				if (uni < 0.f) {
					uni += 1.f;
					carry = twom24;
				} else {
					carry = 0.f;
				}
				seeds[i24 - 1] = uni;
				i24 = next[i24 - 1];
				j24 = next[j24 - 1];
				/* L90: */
			}
		}
		/* L100: */
    }
    kount += *lenv;
    if (kount >= 1000000000) {
		++mkount;
		kount += -1000000000;
    }
    return 0;
	
	/*           Entry to input and float int seeds from previous run */
	
L_rluxin:
    twom24 = 1.f;
    for (i__ = 1; i__ <= 24; ++i__) {
		next[i__ - 1] = i__ - 1;
		/* L195: */
		twom24 *= .5f;
    }
    next[0] = 24;
    twom12 = twom24 * 4096.f;
    printf(" Full initialization of Ranlux with 25 integers: \n");
	for(i__ = 0; i__< 5; ++i__) {
		printf("%d  %d  %d  %d  %d \n", isdext[1+5*i__], isdext[2+5*i__], isdext[3+5*i__], isdext[4+5*i__], isdext[5+5*i__]);
    }
    for (i__ = 1; i__ <= 24; ++i__) {
		seeds[i__ - 1] = (float) isdext[i__] * twom24;
		/* L200: */
    }
    carry = 0.f;
    if (isdext[25] < 0) {
		carry = twom24;
    }
    isd = abs(isdext[25]);
    i24 = isd % 100;
    isd /= 100;
    j24 = isd % 100;
    isd /= 100;
    in24 = isd % 100;
    isd /= 100;
    luxlev = isd;
    if (luxlev <= 4) {
		nskip = ndskip[luxlev];
		printf("ranlux luxury level set by rluxin to: %d\n", luxlev);
    } else if (luxlev >= 24) {
		nskip = luxlev - 24;
		printf("ranlux p-value set by rluxin to: %d\n", luxlev);
    } else {
		nskip = ndskip[4];
		printf("ranlux ILLEGAL LUXURY rluxin: %d\n", luxlev);
		luxlev = 4;
    }
    inseed = -1;
	notyet = 0;
    return 0;
	
	/*                    Entry to ouput seeds as ints */
	
L_rluxut:
    for (i__ = 1; i__ <= 24; ++i__) {
		isdext[i__] = (int) (seeds[i__ - 1] * 4096.f * 4096.f);
		/* L300: */
    }
    isdext[25] = i24 + j24 * 100 + in24 * 10000 + luxlev * 1000000;
    if (carry > 0.f) {
		isdext[25] = -isdext[25];
    }
    return 0;
	
	/*                    Entry to output the "convenient" restart point */
	
L_rluxat:
    *lout = luxlev;
    *inout = inseed;
    *k1 = kount;
    *k2 = mkount;
    return 0;
	
	/*                    Entry to initialize from one or three ints */
	
L_rluxgo:
    if (*lux < 0) {
		luxlev = 3;
    } else if (*lux <= 4) {
		luxlev = *lux;
    } else if (*lux < 24 || *lux > 2000) {
		luxlev = 4;
		printf("ranlux ILLEGAL LUXURY rluxgo: %d\n", *lux);
    } else {
		luxlev = *lux;
		for (ilx = 0; ilx <= 4; ++ilx) {
			if (*lux == ndskip[ilx] + 24) {
				luxlev = ilx;
			}
			/* L310: */
		}
    }
    if (luxlev <= 4) {
		nskip = ndskip[luxlev];
		i__1 = nskip + 24;
		printf("ranlux luxury level set by rluxgo: %d, p = %d\n", luxlev, i__1);
    } else {
		nskip = luxlev - 24;
		printf("ranlux luxury level set by rluxgo to: %d\n", luxlev);
    }
    in24 = 0;
    if (*ins < 0) {
		printf("Illegal initialization by RLUXGO, negative input seed\n");
    }
    if (*ins > 0) {
		jseed = *ins;
		printf("ranlux initialized by rluxgo from seeds: %d, %d, %d \n", jseed, *k1, *k2);
    } else {
		jseed = 314159265;
		printf("ranlux initialized by rluxgo from default seed \n");
    }
    inseed = jseed;
    notyet = 0;
    twom24 = 1.f;
    for (i__ = 1; i__ <= 24; ++i__) {
		twom24 *= .5f;
		k = jseed / 53668;
		jseed = (jseed - k * 53668) * 40014 - k * 12211;
		if (jseed < 0) {
			jseed += 2147483563;
		}
		iseeds[i__ - 1] = jseed % 16777216;
		/* L325: */
    }
    twom12 = twom24 * 4096.f;
    for (i__ = 1; i__ <= 24; ++i__) {
		seeds[i__ - 1] = (float) iseeds[i__ - 1] * twom24;
		next[i__ - 1] = i__ - 1;
		/* L350: */
    }
    next[0] = 24;
    i24 = 24;
    j24 = 10;
    carry = 0.f;
    if (seeds[23] == 0.f) {
		carry = twom24;
    }
	/*        If restarting at a break point, skip K1 + IGIGA*K2 */
	/*        Note that this is the number of numbers delivered to */
	/*        the user PLUS the number skipped (if luxury .GT. 0). */
    kount = *k1;
    mkount = *k2;
    if (*k1 + *k2 != 0) {
		i__1 = *k2 + 1;
		for (iouter = 1; iouter <= i__1; ++iouter) {
			inner = 1000000000;
			if (iouter == *k2 + 1) {
				inner = *k1;
			}
			i__2 = inner;
			for (isk = 1; isk <= i__2; ++isk) {
				uni = seeds[j24 - 1] - seeds[i24 - 1] - carry;
				if (uni < 0.f) {
					uni += 1.f;
					carry = twom24;
				} else {
					carry = 0.f;
				}
				seeds[i24 - 1] = uni;
				i24 = next[i24 - 1];
				j24 = next[j24 - 1];
				/* L450: */
			}
			/* L500: */
		}
		/*         Get the right value of IN24 by direct calculation */
		in24 = kount % (nskip + 24);
		if (mkount > 0) {
			izip = 1000000000 % (nskip + 24);
			izip2 = mkount * izip + in24;
			in24 = izip2 % (nskip + 24);
		}
		/*       Now IN24 had better be between zero and 23 inclusive */
		if (in24 > 23) {
			printf("Error in RESTARTING with RLUXGO: \n");
			printf("The values %d, %d, %d cannot occur at luxury level %d\n", *ins, *k1, *k2, luxlev);
			in24 = 0;
		}
    }
    return 0;
} /* ranlux_ */

/* Subroutine */ int ranlux_(float *rvec, int *lenv)
{
    return ranlux_0_(0, rvec, lenv, (int *)0, (int *)0, (int *)0, 
					 (int *)0, (int *)0, (int *)0, (int *)0);
}

/* Subroutine */ int rluxin_(int *isdext)
{
    return ranlux_0_(1, (float *)0, (int *)0, isdext, (int *)0, (
																 int *)0, (int *)0, (int *)0, (int *)0, (int *)
					 0);
}

/* Subroutine */ int rluxut_(int *isdext)
{
    return ranlux_0_(2, (float *)0, (int *)0, isdext, (int *)0, (
																 int *)0, (int *)0, (int *)0, (int *)0, (int *)
					 0);
}

/* Subroutine */ int rluxat_(int *lout, int *inout, int *k1, 
							 int *k2)
{
    return ranlux_0_(3, (float *)0, (int *)0, (int *)0, lout, inout, 
					 k1, k2, (int *)0, (int *)0);
}

/* Subroutine */ int rluxgo_(int *lux, int *ins, int *k1, int 
							 *k2)
{
    return ranlux_0_(4, (float *)0, (int *)0, (int *)0, (int *)0, (
																   int *)0, k1, k2, lux, ins);
}

/* Subroutine */ int rnpssn_0_(int n__, float *amu, int *n, int *ierr, float *amx)
{
    /* Initialized data */
	
    static float amu0 = -12345.67f;
    static float amax = 88.f;
	
	static int c__1 = 1;
	
    /* Builtin functions */
	
    /* Local variables */
    static int j;
    static float p, r__, emu;
    extern /* Subroutine */ int ranlux_(float *, int *), rnormx_(float *, 
																 int *, int (*)(float *, int *));
	
    switch(n__) {
		case 1: goto L_rnpset;
	}
	
    *ierr = 0;
    if (*amu <= 0.f) {
		*ierr = 1;
		j = 0;
    } else if (*amu > amax) {
		rnormx_(&r__, &c__1, ranlux_);
		j = r__ * sqrt(*amu) + *amu + .5f;
    } else {
		if (*amu != amu0) {
			amu0 = *amu;
			emu = exp(-(*amu));
		}
		p = 1.f;
		j = -1;
	L1:
		++j;
		ranlux_(&r__, &c__1);
		p *= r__;
		if (p > emu) {
			goto L1;
		}
    }
	/* PN */
    if (j < 0) {
		printf("RNPSSN: Warning: J<0; J = %d\n", j);
		printf("Correction: J=0\n");
		printf("Increase AMAX value!");
	    j = 0;
    }
	/* PN */
    *n = j;
    return 0;
	
L_rnpset:
    amax = fminf(*amx,88.f);
	printf("+++++ CERN V136 RNPSSN: SWITCH TO NORMAL APPROXIMATION FOR AMU > %f \n", amax);
    return 0;
} /* rnpssn_ */

/* Subroutine */ int rnpssn_(float *amu, int *n, int *ierr)
{
    return rnpssn_0_(0, amu, n, ierr, (float *)0);
}

/* Subroutine */ int rnpset_(float *amx)
{
    return rnpssn_0_(1, (float *)0, (int *)0, (int *)0, amx);
}

/* Subroutine */ int rnormx_(float *devias, int *ndev, int (*routin)(float *, int *))
{
    /* Initialized data */
	
    static float s = .449871f;
    static float t = -.386595f;
    static float a = .196f;
    static float b = .25472f;
    static float r1 = .27597f;
    static float r2 = .27846f;
	
	static int c__2 = 2;
	
    /* System generated locals */
    int i__1;
    float r__1, r__2;
	
	
    /* Local variables */
    static float q, u[2], v, x, y;
    static int idev;
    static float deviat;
	
	/*        Generator of a vector of independent Gaussian-distributed */
	/*        (pseudo-)random numbers, of mean zero and variance one, */
	/*        making use of a uniform pseudo-random generator (RANMAR). */
	/*        The algorithm for converting uniform numbers to Gaussian */
	/*        is that of "Ratio of Uniforms with Quadratic Bounds."  The */
	/*        method is in principle exact (apart from rounding errors), */
	/*        and is based on the variant published by Joseph Leva in */
	/*        ACM TOMS vol. 18(1992), page 449 for the method and 454 for */
	/*        the Fortran algorithm (ACM No. 712). */
	/*        It requires at least 2 and on average 2.74 uniform deviates */
	/*        per Gaussian (normal) deviate. */
	/*   WARNING -- The uniform generator should not produce exact zeroes, */
	/*   since the pair (0.0, 0.5) provokes a floating point exception. */
    /* Parameter adjustments */
    --devias;
	
    /* Function Body */
	/*         generate pair of uniform deviates */
    i__1 = *ndev;
    for (idev = 1; idev <= i__1; ++idev) {
	L50:
		(*routin)(u, &c__2);
		v = (u[1] - .5f) * 1.7156f;
		x = u[0] - s;
		y = fabsf(v) - t;
		/* Computing 2nd power */
		r__1 = x;
		q = r__1 * r__1 + y * (a * y - b * x);
		/*           accept P if inside inner ellipse */
		if (q < r1) {
			goto L100;
		}
		/*           reject P if outside outer ellipse */
		if (q > r2) {
			goto L50;
		}
		/*           reject P if outside acceptance region */
		/* Computing 2nd power */
		r__1 = v;
		/* Computing 2nd power */
		r__2 = u[0];
		if (r__1 * r__1 > log(u[0]) * -4.f * (r__2 * r__2)) {
			goto L50;
		}
		/*           ratio of P's coordinates is normal deviate */
	L100:
		deviat = v / u[0];
		/* L200: */
		devias[idev] = deviat;
    }
    return 0;
} /* rnormx_ */




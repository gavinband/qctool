
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

/*
* strtod.cpp
*
* This is a version of strtod() modified to be safer and work with the genfile::string_utils library.
* It originated in a version from Ruby's implementation.
* The original copyright notice follows.
*/

/* 
 * strtod.c --
 *
 *	Source code for the "strtod" library procedure.
 *
 * Copyright (c) 1988-1993 The Regents of the University of California.
 * Copyright (c) 1994 Sun Microsystems, Inc.
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any purpose and without
 * fee is hereby granted, provided that the above copyright
 * notice appear in all copies.	 The University of California
 * makes no representations about the suitability of this
 * software for any purpose.  It is provided "as is" without
 * express or implied warranty.
 *
 */

/*
 * Modifications to strtod for slices of non-null terminated (C++-style) strings.
 * Copyright (c) 2011 Gavin Band
*/

#include <limits>
#include "genfile/string_utils/strtod.hpp"
#include "genfile/string_utils/slice.hpp"
#include "genfile/string_utils/string_utils.hpp"

namespace genfile {
	namespace string_utils {
		namespace {
			static int maxExponent = 511;	/* Largest possible base 10 exponent.  Any
							 * exponent larger than this will already
							 * produce underflow or overflow, so there's
							 * no need to worry about additional digits.
							 */
			static double powersOf10[] = {	/* Table giving binary powers of 10.  Entry */
				10.,			/* is 10^2^i.  Used to convert decimal */
				100.,			/* exponents into floating-point numbers. */
				1.0e4,
				1.0e8,
				1.0e16,
				1.0e32,
				1.0e64,
				1.0e128,
				1.0e256
			};
		}

		double strtod( slice const& string ) {
			/*
			 *----------------------------------------------------------------------
			 *
			 * strtod --
			 *
			 *	This procedure converts a floating-point number from an ASCII
			 *	decimal representation to internal double-precision format.
			 *
			 * Results:
			 *	The return value is the double-precision floating-point
			 *	representation of the characters in string.	 Unlike the
			 *	usual strtod() implementations, initial whitespace and
			 *	trailing characters are not allowed.
			 *	A StringConversionError is thrown on error.
			 *
			 * Arguments:
			 *
			 * a_slice: A decimal ASCII floating-point number.
			 * Must have form "-I.FE-X", where I is the
			 * integer part of the mantissa, F is the
			 * fractional part of the mantissa, and X
			 * is the exponent.	 Either of the signs
			 * may be "+", "-", or omitted.	 Either I
			 * or F may be omitted, or both.  The decimal
			 * point isn't necessary unless F is present.
			 * The "E" may actually be an "e".	E and X
			 * may both be omitted (but not just one).
			 */
			{
				bool sign, expSign = false;
				double fraction, dblExp, *d;
				register int c;
				int exp = 0;		/* Exponent read from "EX" field. */
				int fracExp = 0;		/* Exponent that derives from the fractional
							 * part.  Under normal circumstatnces, it is
							 * the negative of the number of digits in F.
							 * However, if I is very long, the last digits
							 * of I get dropped (otherwise a long I with a
							 * large negative exponent could cause an
							 * unnecessary overflow on I alone).  In this
							 * case, fracExp is incremented one for each
							 * dropped digit. */
				int mantSize;		/* Number of digits in mantissa. */
				int decPt;			/* Number of mantissa digits BEFORE decimal
							 * point. */
				const char *pExp;		/* Temporarily holds location of exponent
							 * in string. */

				register const char *p = &string[0];
				register const char * const end_p = &string[0] + string.size() ;
			
				if( p == end_p ) {
					throw StringConversionError() ;
				}

				if (*p == '-') {
					sign = true;
					p += 1;
				} else {
					if (*p == '+') {
						p += 1;
					}
					sign = false;
				}

				if( p == end_p ) {
					throw StringConversionError() ;
				}

				if(
					(end_p - p) == 3
					&& ( *p == 'n' || *p == 'N' )
					&& ( *(p+1) == 'a' || *(p+1) == 'A' )
					&& ( *(p+2) == 'n' || *(p+2) == 'N' )
				) {
					return std::numeric_limits< double >::quiet_NaN() ;
				}

				if(
					(end_p - p) == 3
					&& ( *p == 'i' || *p == 'I' )
					&& ( *(p+1) == 'n' || *(p+1) == 'N' )
					&& ( *(p+2) == 'f' || *(p+2) == 'F' )
				) {
					return (( sign ) ? -1.0 : 1.0 ) * std::numeric_limits< double >::infinity() ;
				}
				

				/*
				 * Count the number of digits in the mantissa (including the decimal
				 * point), and also locate the decimal point.
				 */

				decPt = -1;
				for (mantSize = 0; ; ++mantSize, ++p ) {
					if( p == end_p ) {
						break ;
					}
					else  {
						c = *p;
						if (!isdigit(c)) {
							if ((c != '.') || (decPt >= 0)) {
								break;
							}
							decPt = mantSize;
						}
					}
				}

				/*
				 * Now suck up the digits in the mantissa.	Use two integers to
				 * collect 9 digits each (this is faster than using floating-point).
				 * If the mantissa has more than 18 digits, ignore the extras, since
				 * they can't affect the value anyway.
				 */

				pExp  = p;
				p -= mantSize;
				if (decPt < 0) {
					decPt = mantSize;
				} else {
					mantSize -= 1;			/* One of the digits was the point. */
				}
				if (mantSize > 18) {
					fracExp = decPt - 18;
					mantSize = 18;
				} else {
					fracExp = decPt - mantSize;
				}
				if (mantSize == 0) {
					fraction = 0.0;
					p = &string[0] ;
					goto done;
				} else {
					int frac1, frac2;
					frac1 = 0;
					for ( ; mantSize > 9; --mantSize ) {
						c = *p;
						p += 1;
						if (c == '.') {
						c = *p;
						p += 1;
						}
						frac1 = 10*frac1 + (c - '0');
					}
					frac2 = 0;
					for (; mantSize > 0; --mantSize ) {
						c = *p;
						p += 1;
						if (c == '.') {
						c = *p;
						p += 1;
						}
						frac2 = 10*frac2 + (c - '0');
					}
					fraction = (1.0e9 * frac1) + frac2;
				}

				/*
				 * Skim off the exponent.
				 */

				p = pExp;
				if ( (p != end_p ) && ((*p == 'E') || (*p == 'e'))) {
					p += 1;
					if( p == end_p ) {
						throw StringConversionError() ;
					}
					if (*p == '-') {
						expSign = true;
						p += 1;
					} else {
						if (*p == '+') {
						p += 1;
						}
						expSign = false;
					}
					if( p == end_p ) {
						throw StringConversionError() ;
					}
					for( ; p != end_p && isdigit(*p); ++p ) {
						exp = exp * 10 + (*p - '0');
					}
				}
				if (expSign) {
					exp = fracExp - exp;
				} else {
					exp = fracExp + exp;
				}

				/*
				 * Generate a floating-point number that represents the exponent.
				 * Do this by processing the exponent one bit at a time to combine
				 * many powers of 2 of 10. Then combine the exponent with the
				 * fraction.
				 */
				if (exp < 0) {
					expSign = true;
					exp = -exp;
				} else {
					expSign = false;
				}
				if (exp > maxExponent) {
					exp = maxExponent;
					// errno = ERANGE;
				}
				dblExp = 1.0;
				for (d = powersOf10; exp != 0; exp >>= 1, d += 1) {
					if (exp & 01) {
						dblExp *= *d;
					}
				}
				if (expSign) {
					fraction /= dblExp;
				} else {
					fraction *= dblExp;
				}

			done:
				if( p != end_p ) {
					throw StringConversionError() ;
				}

				if (sign) {
					return -fraction;
				}
				return fraction;
			}
		}
	}
}

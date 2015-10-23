#ifndef PHYSICAL_CONSTANTS_H
#define PHYSICAL_CONSTANTS_H

/**
 * Namespace containing physical constants
 *
 * All members of this namespace should be given in SI standard units,
 * i.e. in terms of kg, m, s, K, mol, A, cd.
 */
namespace Phys
{

/**
  Ideal gas constant in SI standard units (J mol^-1 K^-1)

  Source:               http://physics.nist.gov/cgi-bin/cuu/Value?r
  Date visited:         2015-04-17
  Standard uncertainty: 0.000 0075 J mol^-1 K^-1
*/
const double R = 8.3144621;


/**
 * Molar masses of certain elements and chemical compounds
 */
namespace MolMass
{
	/**
	 * Water
	 *
	 * Source:       http://www.ciaaw.org/pubs/TSAW2013_xls.xls
	 * Date visited: 2015-05-28
	 *
	 * According to the IUPAC report the molar mass of O is in the range [15.999 03, 15.999 77] g/mol
	 * and the molar mass of H is in the range [1.007 84, 1.008 11] g/mol
	 */
	const double Water = 0.018016; ///< kg mol^-1

	/**
	 * N_2
	 *
	 * Source:       http://www.ciaaw.org/pubs/TSAW2013_xls.xls
	 * Date visited: 2015-05-28
	 *
	 * According to the IUPAC report the molar mass of N is in the range [14.006 43, 14.007 28] g/mol
	 */
	const double N2    = 0.028013; ///< kg mol^-1

	/**
	 * O_2
	 *
	 * Source:       http://www.ciaaw.org/pubs/TSAW2013_xls.xls
	 * Date visited: 2015-05-28
	 *
	 * According to the IUPAC report the molar mass of O is in the range [15.999 03, 15.999 77] g/mol
	 */
	const double O2    = 0.032;    ///< kg mol^-1
}

}

#endif // PHYSICAL_CONSTANTS_H

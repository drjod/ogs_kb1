
#include "ogs_display.h"

#include <cstdarg>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "Configure.h"


/**************************************************************************/
/* ROCKFLOW - Funktion: DisplayStartMsg
 */
/* Aufgabe:
   Gibt Eroeffnungsbildschirm aus
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   - void -
 */
/* Ergebnis:
   - void -
 */
/* Programmaenderungen:
   03/1994     MSR        Erste Version
   01/2010     NW         Automatic centering of the version information
 */
/**************************************************************************/
void DisplayStartMsg ( void )
{
	int i, pad_len;
	char buf[128];

	char version[] = OGS_VERSION;

	printf("\n");
	printf("     ####################################################\n");
	printf("     #                                                  #\n");
	printf("     #    ###    ###    ###        #    #  ####     #   #\n");
	printf("     #   #   #  #   #  #   #       #   #   #   #   ##   #\n");
	printf("     #   #   #  #       #          #  #    #   #  # #   #\n");
	printf("     #   #   #  #        #    ###  ###     ####     #   #\n");
	printf("     #   #   #  #  ##     #        #  #    #   #    #   #\n");
	printf("     #   #   #  #   #  #   #       #   #   #   #    #   #\n");
	printf("     #    ###    ###    ###        #    #  ####     #   #\n");
	printf("     #                                                  #\n");
	printf("     #                                                  #\n");
	printf("     ####################################################\n\n");
	printf("     %s\n\n", version);
	printf("     File name (without extension): ");
}

/**************************************************************************/
/* ROCKFLOW - Funktion: DisplayEndMsg
 */
/* Aufgabe:
   Gibt Programm-Abspann aus.
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   - void -
 */
/* Ergebnis:
   - void -
 */
/* Programmaenderungen:
   03/1994     MSR        Erste Version
 */
/**************************************************************************/
void DisplayEndMsg ( void )
{
	printf("\n          Programm beendet!\n\n\n");
}

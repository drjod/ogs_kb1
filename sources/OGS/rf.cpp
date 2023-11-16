/**************************************************************************/
/* ROCKFLOW - Modul: rf.c
 */
/* Aufgabe:
   ROCKFLOW-FEM - Hauptprogramm
 */
/* Programmaenderungen:
   07/1996     MSR        Erste Version
   06/1998     AH         Konfigurationsdatei
   08/1999     OK         RF-FEM Applikation
   10/1999     AH         Systemzeit

   last modified: OK 14.12.1999
 */
/**************************************************************************/

/**
 * the preprocessor directive RFW_FRACTURE is only useable until version 4.11 of OGS
 * */
#include "Configure.h"
#include "logger.h"
#include <stdexcept>
Logger logger = Logger::get_instance();


#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) || \
        defined(USE_MPI_GEMS) || defined(USE_MPI_KRC) 
#include "par_ddc.h"
#include <mpi.h>
#endif
#ifdef LIS
#include "lis.h"
#include <omp.h>
#endif

#include "logger.h"

#include "BuildInfo.h"

/* Preprozessor-Definitionen */
#include "makros.h"
#include "display.h"
#include "memory.h"
#include "ogs_display.h"
#define TEST
/* Benutzte Module */
#include "break.h"
#include "timer.h"
//16.12.2008. WW #include "rf_apl.h"
#include "FileTools.h"
#include "files0.h"
#ifdef SUPERCOMPUTER
// kg44 test for buffered outputh
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#endif
#include "problem.h"


/* Deklarationen */
int main ( int argc, char* argv[] );
void ShowSwitches ( void );
// LB,string FileName; //WW
// LB,string FilePath; //23.02.2009. WW
// ------  12.09.2007 WW:
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) || \
        defined(USE_MPI_GEMS) || defined(USE_MPI_KRC) 
double elapsed_time_mpi;
// ------
#endif


// Use PETSc. WW
#ifdef USE_PETSC
#include "petscksp.h"
#ifdef USEPETSC34
#include "petsctime.h"
#endif
#endif

// declaration of defaultOutputPath
#include "rf_out_new.h"




/* Definitionen */

/**************************************************************************/
/* ROCKFLOW - Funktion: main
 */
/* Aufgabe:
   Hauptprogramm
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E int argc: Anzahl der Kommandozeilenparameter (incl. Programmname)
   E char *argv[]: Zeiger auf Feld der argc Kommandozeilenparameter
 */
/* Ergebnis:
   Fehlerfreie Bearbeitung: Exit-Code 0
 */
/* Programmaenderungen:
   07/1996     MSR        Erste Version
   08/1999     OK         RF-FEM Applikation
 */
/**************************************************************************/
int main ( int argc, char* argv[] )
{
	/* parse command line arguments */
	std::string anArg;
	std::string modelRoot;

	for( int i = 1; i < argc; i++ )
	{
		anArg = std::string( argv[i] );
		if( anArg == "--help" || anArg == "-h")
		{
			std::cout << "Usage: ogs [MODEL_ROOT] [OPTIONS]\n"
			          << "Where OPTIONS are:\n"
			          << "  -h [--help]               print this message and exit\n"
			          << "  -b [--build-info]         print build info and exit\n"
			          << "  --output-directory DIR    put output files into DIR\n"
			          << "  --version                 print ogs version and exit" << "\n";
			continue;
		}
		if( anArg == "--build-info" || anArg == "-b" )
		{
			std::cout << "ogs version: " << OGS_VERSION << "\n"
			          << "ogs date: " << OGS_DATE << "\n";
#ifdef CMAKE_CMD_ARGS
			std::cout << "cmake command line arguments: " << CMAKE_CMD_ARGS << "\n";
#endif // CMAKE_CMD_ARGS
#ifdef GIT_COMMIT_INFO
			std::cout << "git commit info: " << GIT_COMMIT_INFO << "\n";
#endif // GIT_COMMIT_INFO
#ifdef SVN_REVISION
			std::cout << "subversion info: " << SVN_REVISION << "\n";
#endif // SVN_REVISION
#ifdef BUILD_TIMESTAMP
			std::cout << "build timestamp: " << BUILD_TIMESTAMP << "\n";
#endif // BUILD_TIMESTAMP
			continue;
		}
		if( anArg == "--version" )
		{
			std::cout << OGS_VERSION << "\n";
			continue;
		}
		if( anArg == "--model-root" || anArg == "-m" )
		{
			if (i+1 >= argc) {
				std::cerr << "Error: Parameter " << anArg << " needs an additional argument\n";
				std::exit(EXIT_FAILURE);
			}
			modelRoot = std::string( argv[++i] );
			continue;
		}
		if (anArg == "--output-directory")
		{
			if (i+1 >= argc) {
				std::cerr << "Error: Parameter " << anArg << " needs an additional argument\n";
				std::exit(EXIT_FAILURE);
			}
			std::string path = argv[++i];

			if (! path.empty()) defaultOutputPath = path;
			continue;
		}
		// anything left over must be the model root, unless already found
		if ( modelRoot == "" )
			modelRoot = std::string( argv[i] );
	} // end of parse argc loop

	if( argc > 1 && modelRoot == "" ) // non-interactive mode and no model given
		exit(0);             // e.g. just wanted the build info

	char* dateiname(NULL);
#ifdef SUPERCOMPUTER
// *********************************************************************
// buffered output ... important for performance on cray
// (unbuffered output is limited to 10 bytes per second)
// georg.kosakowski@psi.ch 11.10.2007

	char buf[1024 * 1024];
	int bsize;

	bsize = 1024 * 1024; // question: what happens if buffer is full?
	                     // according to documentation the buffer is flushed when full.
	                     // If we have a lot of output, increasing buffer is usefull.
	if(bsize > 0)
//        bufstd = malloc(bsize);
		setvbuf(stdout, buf, _IOFBF, bsize);
	//**********************************************************************
#endif
/*---------- MPI Initialization ----------------------------------*/
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) || \
	defined(USE_MPI_GEMS) || defined(USE_MPI_KRC) 
	printf("Before MPI_Init\n");
#if defined(USE_MPI_GEMS)
	int prov;
	MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED, &prov);
#else
	MPI_Init(&argc,&argv);
#endif
	MPI_Barrier (MPI_COMM_WORLD); // 12.09.2007 WW
	elapsed_time_mpi = -MPI_Wtime(); // 12.09.2007 WW
	MPI_Comm_size(MPI_COMM_WORLD,&mysize);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	std::cout << "After MPI_Init myrank = " << myrank << '\n';
	time_ele_paral = 0.0;
#endif
/*---------- MPI Initialization ----------------------------------*/


#ifdef USE_PETSC
	int rank, r_size;
	PetscLogDouble v1,v2;
	char help[] = "OGS with PETSc \n";
	//PetscInitialize(argc, argv, help);
	PetscInitialize(&argc,&argv,(char *)0,help);
	//kg44 quick fix to compile PETSC with version PETSCV3.4
#ifdef USEPETSC34
       PetscTime(&v1);
#else
       PetscGetTime(&v1);
#endif
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &r_size);
	PetscSynchronizedPrintf(PETSC_COMM_WORLD, "===\nUse PETSc solver");
	PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Number of CPUs: %d, rank: %d\n", r_size, rank);
#endif




/*---------- LIS solver -----------------------------------------*/
#ifdef LIS
	//Initialization of the lis solver.
	lis_initialize(&argc, &argv);
#endif
/*========================================================================*/
/* Kommunikation mit Betriebssystem */
	/* Timer fuer Gesamtzeit starten */
#ifdef TESTTIME
	TStartTimer(0);
#endif
	/* Intro ausgeben */
#if defined(USE_MPI) //WW
	if(myrank == 0)
#endif
#ifdef USE_PETSC
        if(rank == 0 )
#endif

	DisplayStartMsg();
	/* Speicherverwaltung initialisieren */
	if (!InitMemoryTest())
	{
		DisplayErrorMsg("Fehler: Speicherprotokoll kann nicht erstellt werden!");
		DisplayErrorMsg("        Programm vorzeitig beendet!");
		return 1; // LB changed from 0 to 1 because 0 is indicating success
	}
	if( argc == 1 )               // interactive mode

		dateiname = ReadString();
	else                         // non-interactive mode
	{
		if ( argc == 2 )     // a model root was supplied
		{
			dateiname = (char*) Malloc((int)strlen(argv[1]) + 1);
			dateiname = strcpy(dateiname,argv[1]);
		}
		else                // several args supplied
		if( modelRoot != "")
		{
			dateiname = (char*) Malloc( (int) modelRoot.size() + 1 );
			dateiname = strcpy( dateiname, modelRoot.c_str() );
		}
		DisplayMsgLn(dateiname);
	}
	//WW  DisplayMsgLn("");
	//WW  DisplayMsgLn("");
	// ----------23.02.2009. WW-----------------

	// LB Check if file exists
	std::string tmpFilename = dateiname;
	tmpFilename.append(".pcs");
	if(!IsFileExisting(tmpFilename))
	{
		std::cout << " Error: Cannot find file " << dateiname << "\n";
		return 1;
	}

	// If no option is given, output files are placed in the same directory as the input files
	if (defaultOutputPath.empty()) defaultOutputPath = pathDirname(std::string(dateiname));

	FileName = dateiname;
	size_t indexChWin, indexChLinux;
	indexChWin = indexChLinux = 0;
	indexChWin = FileName.find_last_of('\\');
	indexChLinux = FileName.find_last_of('/');
	//
	if(indexChWin != std::string::npos)
		FilePath = FileName.substr(0,indexChWin) + "\\";
	else if(indexChLinux != std::string::npos)
		FilePath = FileName.substr(0,indexChLinux) + "/";

	logger.initialize(FilePath);
	// ---------------------------WW
	Problem* aproblem = new Problem(dateiname);

	try
	{

#ifdef USE_PETSC
		aproblem->setRankandSize(rank, r_size);
#endif
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) || defined(USE_MPI_GEMS)  || defined(USE_MPI_KRC)
		aproblem->setRankandSize(myrank, mysize);
#endif

		aproblem->Euler_TimeDiscretize();
		if(aproblem->PrintTimes())
		  if(ClockTimeVec.size()>0){
			ClockTimeVec[0]->PrintTimes();  //CB time
			DestroyClockTime();
		  }
		delete aproblem;
		aproblem = NULL;
	}
	catch (const std::runtime_error& re)
	{
			std::cout << re.what() << std::endl;
			std::cerr << re.what() << std::endl;
	}

#ifdef TESTTIME
#if defined(USE_MPI)
     if(myrank == 0)
#endif
#if defined(USE_PETSC) 
     if(rank == 0)
#endif
	std::cout << "Simulation time: " << TGetTimer(0) << "s" << "\n";
#endif
	/* Abspann ausgeben */
/*--------- MPI Finalize ------------------*/
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) || defined(USE_MPI_KRC)
	elapsed_time_mpi += MPI_Wtime(); // 12.09.2007 WW
	std::cout << "\n *** Total CPU time of parallel modeling: " << elapsed_time_mpi <<
	"\n";                                                                          //WW
	// Count CPU time of post time loop WW
	MPI_Finalize();
#endif
/*--------- MPI Finalize ------------------*/
/*--------- LIS Finalize ------------------*/
#ifdef LIS
	lis_finalize();
#endif
/*--------- LIS Finalize ------------------*/

	free(dateiname);

#ifdef USE_PETSC
	//kg44 quick fix to compile PETSC with version PETSCV3.4
#ifdef USEPETSC34
       PetscTime(&v2);
#else
       PetscGetTime(&v2);
#endif


   PetscPrintf(PETSC_COMM_WORLD,"\t\n>>Total elapsed time by using PETSC:%f s\n",v2-v1);

   PetscFinalize();
#endif

	return 0;
}

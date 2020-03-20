#!/bin/bash

#############################################################################################  
# To generate ogs binaries with intel compiler (mpi, openmp, mkl) on:
#    1. Lockstedt GPI server - Eclipse IDE
#    2. RZ CLUSTER - Eclipse IDE, PETSC
#    3. NEC cluster - PETSC            
#                                                                    by JOD 10/2015  
#    parameter:
#       $1: path to ogs folder or empty (script steps into build folder for cmake and make)
#       $2: configurationSELECTED 
#                   0: OGS_FEM
#                   1: OGS_FEM_SP
#                   2: OGS_FEM_MKL (sequential)
#                   3: OGS_FEM_MPI
#                   4: OGS_FEM_MPI_KRC
#                   5: OGS_FEN_PETSC (parallel)
#                            (more can be added in 1./3. below)
#       $3: BUILD_CONFIGURATION   [Debug, Release]   No Debug for NEC 
#                                 if BUILD_CONFIGURATION preselected
#       $4: for BUILD_flag, if $4 = build BUILD_flag is set 1 (and build executed), else 0
#
#############################################################################################  
#  
# USER GUIDE: 
#     Put in OGS folder (same level as sources and libs)  
#     cd into this folder and type ./compileInKiel.sh  
#     as an option, call script from somewhere else and pass ogs folder as argument $1
# 


#############################################################################################
# 0. output on console and in compilation.log
#
# parameters:
#     $1 ["ERROR", "WARNING", "INFO"] message type 
#     $2 "message"
#

printMessage()
{
    case $1 in
        "ERROR")
            tput setaf 1;  # red
            ;;
        "WARNING")
            tput setaf 5;  # cyan
            ;;
        "INFO")
            tput setaf 2;  # green
            ;;
        *)  
            ;;
    esac
    
    echo -e "\n"$1 - $2
    tput sgr0;  # reset color    

    echo -e $(date) $1 - $2 >> $OGS_FOLDER/compilation.log
}

############################################################################################# 
# 1. Initialization
#     first function called in main
#    parameter:
#        $1: path to ogs folder 
 
initialize()
{
    nCPUs=6      # number of CPUs for compilation (<= number of nodes on cluster login node or server)
    CALLEDFROM=$PWD         

    if [ -z $1 ]; then
        OGS_FOLDER=$CALLEDFROM # call from ogs folder
    else
        OGS_FOLDER=$1          # path passed as parameter into script
    fi
    
    printMessage "INFO" "Up to compile ${OGS_FOLDER##*/}"
    
    cConfigurations=(  # extend compiler table (3.) if you add code configurations
        "OGS_FEM"  
        "OGS_FEM_SP"  
        "OGS_FEM_MKL"   
        "OGS_FEM_MPI"  
        "OGS_FEM_MPI_KRC" 
        "OGS_FEM_PETSC" 
    )
    
    configurationSELECTED="" # [OGS_FEM, OGS_FEM_SP, ...]  
    for (( i=0; i<${#cConfigurations[@]}; i++ ))  
    do  
        if [ "${cConfigurations[i]}" == "$2" ]; then 
             configurationSELECTED=$i   
        fi        
    done     
    
     if [ "$3" == "Debug" ]; then
        BUILD_CONFIGURATION="Debug"
        IDE="ECLIPSE"             # [empty,ECLIPSE]
    elif [ "$3" == "Release" ]; then    
        BUILD_CONFIGURATION="Release"
        IDE=""
    else
        BUILD_CONFIGURATION=""    
        IDE=""
    fi    
    
    if [ "$4" == "build" ]; then
        BUILD_flag="1"
    else
        BUILD_flag="0"
    fi

    COMPILER_VERSION=""
 
    # paths to 
    ROOT_FOLDER=""    # where folder for code configurations will be placed 
                         # for specific BUILD_CONFIGURATION (Debug, Release) and COMPILER_VERSION 
                         # ($OGS_FOLDER/"Build_${BUILD_CONFIGURATION}_$COMPILER_VERSION")
    BUILD_FOLDER=""   # where folder for specific BUILD_CONFIGURATION will be placed 
                         #($ROOT_FOLDER/$cConfigurationSELECTED)
    SOFTWARE_FOLDER=""   # where intel folder are
    COMPOSER_ROOT=""   # not used with intel16 compiler on rz cluster
    MPI_ROOT=""
    ICC=""           # intel c compiler
    ICPC=""            # intel c++ compilerx
    MPIICC=""        # mpi c compiler
    MPIICPC=""        # mpi c++ compiler
}
 
############################################################################################# 
#  2. SET PATHS 
#
#  INTEL        : rzcluster (1502, composer_xe_2015.2.164),
#                 NEC cluster(15.0.3, composer_xe_2015.3.187)
#                 Lockstedt (13.1.0, composer_xe_2013.2.146)
#  INTEL mpi    : rzcluster (5.0.3.048), NEC cluster (4.1.1.036), Lockstedt
#  PETSC 3.5.3     : rzcluster (intel14)
#                  NEC cluster (intel14)
#  Requirements:
#    initialization of COMPILER_VERSION (1.)
#  Results:
#    paths to intel tools set (modules loaded)

setPaths()
{
    setPaths__host=$HOSTNAME 

    case ${setPaths__host:0:2} in 
        ca) # rzcluster 
            SOFTWARE_FOLDER="/home/Software"
            
            # COMPILER_VERSION="intel1402"   
            # COMPOSER_ROOT="$SOFTWARE_FOLDER/$COMPILER_VERSION/composer_xe_2013_sp1.2.144"     
            # MPI_ROOT="$SOFTWARE_FOLDER/$COMPILER_VERSION/impi/4.1.3.048"  
 
            # COMPILER_VERSION="intel1502"      
            # COMPOSER_ROOT="$SOFTWARE_FOLDER/$COMPILER_VERSION/composer_xe_2015.2.164"          
            # MPI_ROOT="$SOFTWARE_FOLDER/$COMPILER_VERSION/impi/5.0.3.048" 
            # module load $COMPILER_VERSION   
            
            COMPILER_VERSION="intel19.0.4"      
            COMPOSER_ROOT="$SOFTWARE_FOLDER/intel/$COMPILER_VERSION/usr/compilers_and_libraries_2019/linux"          
		# "$SOFTWARE_FOLDER/$COMPILER_VERSION/compilers_and_libraries_2019.4.243/linux/mpi" 
            MPI_ROOT="$COMPOSER_ROOT/mpi"  

            ICPC=icpc
		# $COMPOSER_ROOT/bin/intel64/icpc"

            MPIICC="$MPI_ROOT/intel64/bin/mpiicc"  
            MPIICPC="$MPI_ROOT/intel64/bin/mpiicpc"  

	    module load cmake/3.15.4                
            module load intel/18.0.4
            module load intelmpi/18.0.4
            #module load intel16.0.0
            #module load intelmpi16.0.0
	    #module load petsc-3.7.5            
            #module load eclipse
            
            MKLROOT="$COMPOSER_ROOT/mkl"   
            export PATH=$PATH:$MKLROOT/lib/intel64
            export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKLROOT/lib/intel64 
            # . $MKLROOT/bin/intel64/mklvars_intel64.sh            
                ;;
        ne) # NEC cluster 
            SOFTWARE_FOLDER="/opt"    
            COMPOSER_ROOT="$SOFTWARE_FOLDER/intel/composer_xe_2015.3.187"
            MPI_ROOT="$SOFTWARE_FOLDER/intel/impi/4.1.1.036"
    
            ICC="$COMPOSER_ROOT/bin/intel64/icc"
            ICPC="$COMPOSER_ROOT/bin/intel64/icpc"
     
            MPIICC="$MPI_ROOT/intel64/bin/mpiicc"  
            MPIICPC="$MPI_ROOT/intel64/bin/mpiicpc"
    
            COMPILER_VERSION="intel15.0.3"
            module load $COMPILER_VERSION    
            module load petsc-3.5.3-intel
            
            MKLROOT="$COMPOSER_ROOT/mkl"                  
            export PATH=$PATH:$MKLROOT/lib/intel64
            export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKLROOT/lib/intel64 
            . $MKLROOT/bin/intel64/mklvars_intel64.sh
                ;;
        Lo) # GPI server 
            SOFTWARE_FOLDER="/opt"
            COMPOSER_ROOT="$SOFTWARE_FOLDER/intel/composer_xe_2013.2.146"
            MPI_ROOT="$SOFTWARE_FOLDER/openmpi"
        
            ICC="$COMPOSER_ROOT/bin/intel64/icc"
            ICPC="$COMPOSER_ROOT/bin/intel64/icpc"
     
            MPIICC="$MPI_ROOT/bin/mpicc"  
            MPIICPC="$MPI_ROOT/bin/mpicxx"        
        
            . $COMPOSER_ROOT/bin/compilervars.sh intel64
            . $COMPOSER_ROOT/mkl/bin/intel64/mklvars_intel64.sh
        
            COMPILER_VERSION="intel" 
                ;;
	je) # my Fujitsu
	    ICC="/usr/bin/gcc"
	    ICPC="/usr/bin/g++"

            export PATH=/usr/bin:$PATH
            export LD_LIBRARY_PATH=/usr/lib

	    COMPILER_VERSION="gnu"
                ;;
	am) # AMAK
	    ICC="/usr/bin/gcc"
	    ICPC="/usr/bin/g++"
	    
	    MPIICC="/usr/bin/mpicc"
	    MPIICPC="/usr/bin/mpicxx"

            export PATH=/usr/bin:$PATH
            export LD_LIBRARY_PATH=/usr/lib

	    . /opt/intel/compilers_and_libraries_2017.0.098/linux/mkl/bin/mklvars.sh intel64

	    export PETSC_DIR=/home/jens/petsc/petsc-3.7.6
	    export PETSC_ARCH=arch-linux2-c-debug

	    COMPILER_VERSION="gnu"
		;;

	co) # DOCKER CONTAINER
            ICC="/usr/bin/gcc"
            ICPC="/usr/bin/g++"
            MPIICC="/usr/local/bin/mpicc"
            MPIICPC="/usr/local/bin/mpicxx"

            export PATH=/usr/bin:/usr/local/bin:$PATH
            export LD_LIBRARY_PATH=/usr/lib:/usr/local/lib:$LD_LIBRARY_PATH

            export PETSC_DIR=/home/ogs_user/petsc-3.7.6
            export PETSC_ARCH=arch-linux2-c-debug
            COMPILER_VERSION="gnu"
                ;;
        *)
            printMessage "ERROR - Check HOSTNAME"
                ;;      
    esac

    printMessage "INFO" "Paths set for $HOSTNAME" 
}
  
############################################################################################# 
# 3. Compiler table
#    here you can add code configurations  
#    !!!!! MATCH LINES OF cConfigurations (1.) AND compilerTable 
#     Requirements:
#         paths set (2.)
#    Result:
#        Compiler are assigned to code configurations (from 1.)
#
 
setCompilerTable()
{ 
    compilerTable=( 
        #    -DPARALLEL_USE_OPENMP=     -DCMAKE_C_COMPILER=     -DCMAKE_CXX_COMPILER=          
            "OFF"                    "$ICC"                    "$ICPC"                    # OGS_FEM   
            "OFF"                    "$ICC"                    "$ICPC"                    # OGS_FEM_SP   
            "ON"                    "$ICC"                    "$ICPC"                    # OGS_FEM_MKL   
            "OFF"                    "$MPIICC"                "$MPIICPC"                # OGS_FEM_MPI  
            "OFF"                    "$MPIICC"                "$MPIICPC"                # OGS_FEM_MPI_KRC 
            "ON"                    "$MPIICC"                "$MPIICPC"                # OGS_FEM_PETSC                      
    )    
}  
    
#############################################################################################  
# 4. User input
#    supported configurations in 1. above 
#    Debug BUILD_CONFIGURATION for ECLIPSE IDE supported for rzcluster and Lokstedt 
#   Requirements:
#        cConfigurations list
#    Results:
#        variables set
#            configurationSELECTED     
#            BUILD_CONFIGURATION        
#            BUILD_flag                
#

selectConfiguration()  # code configuration from list cConfigurations
{      
    echo -e "Select (x for exit)"  
    for (( i=0; i<${#cConfigurations[@]}; i++ ))  
    do  
        echo -e "\t$i: ${cConfigurations[$i]}"  
    done  
    echo -e "\ta: all"  
    read -n1 configurationSELECTED  
    
    if [ $configurationSELECTED != "x" ]; then    
        # exception handling - restart if input error
        if echo $configurationSELECTED | egrep -q '^[0-9]+$'; then 
            if [ "$configurationSELECTED" -lt 0 ] || [ "$configurationSELECTED" -ge ${#cConfigurations[@]} ]; then
                printMessage "ERROR" "Number out of range - Restart"                
                main
            fi
        else
            if [ $configurationSELECTED != "a" ]; then
              printMessage "ERROR" "Input neither a number nor \"a\" to select all - Restart"                
              main
            fi
        fi
    fi
}  

selectBuild()  
{  
    selectBuild__cInput=""  # used as local variable
    selectBuild__host=$HOSTNAME
    
    # configuration
    if [ "${selectBuild__host:0:2}" == "rz" ] || [ "${selectBuild__host:0:2}" == "Lo" ] || [ "${selectBuild__host:0:2}" == "am" ]; then
        # Eclipse exists on RZ cluster, amak and Lokstedt
        echo -e "\n[d]ebug or [r]elease?" 
        read -n1 selectBuild__cInput  
        if [ "$selectBuild__cInput" == "d" ]; then  
            BUILD_CONFIGURATION="Debug"
            IDE="ECLIPSE"
        elif [ "$selectBuild__cInput" == "r" ]; then   
            BUILD_CONFIGURATION="Release"        
        else
            printMessage "ERROR" "Take \"d\" or \"r\" - Restart"                
            main
        fi 
    else
        BUILD_CONFIGURATION="Release"     
    fi    
        
    # flag
    echo -e "\n--------------------------------------------------\n"
    echo -e "\nCreate Build Files ([y]es or [n]o)?"  
    read -n1 selectBuild__cInput  
    if [ "$selectBuild__cInput" == "y" ]; then  
       BUILD_flag=1  
    elif [ "$selectBuild__cInput" == "n" ]; then  
       BUILD_flag=0
    else
        printMessage "ERROR" "Take \"y\" or \"n\" - Restart"        
        main
    fi  
}  

#############################################################################################  
# 5. linking
#     do cmake
#   parameter:
#       $1: main__configurationNDX (6.)  
#     Requirements:
#        compiler paths, compiler table (2.)
#        configurations for code (cConfigurationSELECTED, 1.) and build (BUILD_CONFIGURATION, 4.)    
#     Result:
#        Build directories exist
#

build()
{
    rm -rf $BUILD_FOLDER  # remove old build
    mkdir $BUILD_FOLDER  
    cd $BUILD_FOLDER   # step into build folder for cmake
   
    # local variables
    OPENMP=${compilerTable[(($1 * 3))]}          
    build__COMPILER_C=${compilerTable[(($1 * 3 + 1))]} 
    build__COMPILER_CXX=${compilerTable[(($1 * 3 + 2))]}

    printMessage "INFO" "Building files - Debugger $build__COMPILER_C $build__COMPILER_CXX"
    if [ "$IDE" == "ECLIPSE" ]; then  # only difference is GENERATOR_OPTION -G
        cmake ../../sources -G "Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=$BUILD_CONFIGURATION -D$cConfigurationSELECTED=ON -DPARALLEL_USE_OPENMP=${compilerTable[(($1 * 3))]} -DCMAKE_C_COMPILER=$build__COMPILER_C  -DCMAKE_CXX_COMPILER=$build__COMPILER_CXX                       
    else

        if [ "$build__COMPILER_C" == "" ]; then
           cmake ../../sources -DCMAKE_BUILD_TYPE=$BUILD_CONFIGURATION -D$cConfigurationSELECTED=ON -DPARALLEL_USE_OPENMP=$OPENMP -DCMAKE_CXX_COMPILER=$build__COMPILER_CXX                      
        else
           cmake ../../sources -DCMAKE_BUILD_TYPE=$BUILD_CONFIGURATION -D$cConfigurationSELECTED=ON -DPARALLEL_USE_OPENMP=$OPENMP -DCMAKE_C_COMPILER=$build__COMPILER_C  -DCMAKE_CXX_COMPILER=$build__COMPILER_CXX                      
        fi
    fi
}

#############################################################################################
# 6. main function
#     Calls: 
#        select functions (3.) 
#        build (4.) for cmake if BUILD_flag = 1 or if Build folder does not exist
#         make for compilation
#     Requirements:
#         --- (here script starts)
#     Result:
#         binaries (renamed)
#    parameter:
#        $1: path to ogs folder
#

main()
{
    # config - variables and paths
    initialize $1 $2 $3 $4               
    setPaths
    setCompilerTable
    
    # user input
    if [ -z $2 ]; then # configuration not preselected (by argument)
        selectConfiguration    
    fi
    
    if [ "$configurationSELECTED" != "x" ]; then # else exit
        if [ -z $3 ]; then # build flag not previously set (by argument)
            selectBuild 
        fi
        # config main loop
        ROOT_FOLDER="$OGS_FOLDER/Build_${BUILD_CONFIGURATION}_$COMPILER_VERSION"
        mkdir -p $ROOT_FOLDER
        # loop over all configurations (from list in 1.)
        for (( main__configurationNDX=0; main__configurationNDX<${#cConfigurations[@]}; main__configurationNDX++ ))  
        do  
            # either one or all can be selected
            if [ "$main__configurationNDX" == "$configurationSELECTED" ] || [ "$configurationSELECTED" == "a" ]; then  
                # pre-processing
                cConfigurationSELECTED=${cConfigurations[main__configurationNDX]}   
                printMessage "INFO" "GENERATING $cConfigurationSELECTED $BUILD_CONFIGURATION $IDE" 
                BUILD_FOLDER=$ROOT_FOLDER/$cConfigurationSELECTED
                # build
                if [ "$BUILD_flag" -eq 1 ]; then  
                    build $main__configurationNDX
                else
                    if [ ! -d "$BUILD_FOLDER" ]; then
                        printMessage "WARNING" "Build folder does not exist - Building it now"
                        build $main__configurationNDX
                    fi
                    cd $BUILD_FOLDER # step into build folder for make
                fi  

                # compile
                printMessage "INFO" "Compiling"
                make -j $nCPUs    

		cd ../..  # step into ROOT_FOLDER
                # post-processing
                if [ -e $BUILD_FOLDER/bin/ogs ]; then            
                    mv $BUILD_FOLDER/bin/ogs $BUILD_FOLDER/bin/ogs_$cConfigurationSELECTED     # rename
                    printMessage "INFO" "Binaries ogs_$cConfigurationSELECTED generated"                
                else
                    printMessage "WARNING" "No binaries generated"
                fi
            fi  
        done   

        cd $CALLEDFROM      # back to initial folder
    
        if [ -z $2 ]; then 
            main $1 "" "" # restart - than BUILD_CONFIGURATION, Build_flag and configuration always selected by user
        fi   # else BUILD_CONFIGURATION was preselected - exit now
    fi

} 

if [ "$configurationSELECTED" != "x" ]; then
    main $1 $2 $3 $4 # start
fi # else exit

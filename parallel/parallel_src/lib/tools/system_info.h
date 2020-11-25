/******************************************************************************
 * distributed_quality_metrics.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef SYSTEM_INFO_UAA1EX3T
#define SYSTEM_INFO_UAA1EX3T

#ifdef __GNUC__ //apple clang does not have the sys/* libraries

#include <mpi.h>
#include <string.h>

#include "sys/types.h"
#include "sys/sysinfo.h"
#include "sys/times.h"
#include "sys/vtimes.h"

//taken from 
//https://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
unsigned long long getFreeRam(const MPI_Comm communicator, double& myMem, bool printMessage=false){

    struct sysinfo memInfo;
    PEID rank;
    MPI_Comm_rank( communicator, &rank);
    const double kb = 1024.0;
    const double mb = kb*1024;
    [[maybe_unused]] const double gb = mb*1024;

    sysinfo (&memInfo);
    long long totalVirtualMem = memInfo.totalram;
    //Add other values in next statement to avoid int overflow on right hand side...
    totalVirtualMem += memInfo.totalswap;
    totalVirtualMem *= memInfo.mem_unit;

    long long totalPhysMem = memInfo.totalram;
    //Multiply in next statement to avoid int overflow on right hand side...
    totalPhysMem *= memInfo.mem_unit;

    long long physMemUsed = memInfo.totalram - memInfo.freeram;
    //Multiply in next statement to avoid int overflow on right hand side...
    physMemUsed *= memInfo.mem_unit;

    unsigned long long freeRam = memInfo.freeram;
    freeRam *= memInfo.mem_unit;

    unsigned long long sharedRam = memInfo.sharedram;
    sharedRam *= memInfo.mem_unit;    
    
    unsigned long long buffRam = memInfo.bufferram;
    buffRam *= memInfo.mem_unit; 

    auto parseLine = [](char* line){
        // This assumes that a digit will be found and the line ends in " Kb".
        int i = strlen(line);
        const char* p = line;
        while (*p <'0' || *p > '9') p++;
        line[i-3] = '\0';
        i = atoi(p);
        return i;
    };

    //Note: this value is in KB!
    auto getValue = [&](){ 
        FILE* file = fopen("/proc/self/status", "r");
        int result = -1;
        char line[128];

        while (fgets(line, 128, file) != NULL){
            if (strncmp(line, "VmRSS:", 6) == 0){
                result = parseLine(line);
                break;
            }
        }
        fclose(file);
        return result;
    };

    myMem = getValue()/kb;

    if( printMessage ){
        double maxMem;
        MPI_Reduce( &myMem, &maxMem, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
        if( rank==0 ){
                std::cout << "totalPhysMem: " << (totalPhysMem/mb) << 
                " MB, physMemUsed: " << physMemUsed/mb << 
                " MB, free ram: " << freeRam/mb << 
                " max mem used: " << maxMem << " MB" << std::endl;
        }

    }

    return freeRam;
}

#else //for other compilers

unsigned long getFreeRam(){
    std::cout<< "Memory usage for non-gcc compilers is not supported " <<std::endl;
}

#endif //check for apple clang

#endif /* end of include guard: SYSTEM_INFO_UAA1EX3T */

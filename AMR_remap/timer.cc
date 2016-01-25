/* LANL Copyright 2015
 *
 *
 *  Authors: 
 *        Gerald Collom        gcollom@lanl.gov
 *        Colin Redman        credman@lanl.gov
 */

// Timer routines to simplify main loop

#include <stdlib.h>
#include <sys/time.h>

#include "timer.h"

void cpu_timer_start(struct timeval *tstart_cpu){
   gettimeofday(tstart_cpu, NULL);
}

double cpu_timer_stop(struct timeval tstart_cpu){
   double result;
   struct timeval tstop_cpu, tresult;

   gettimeofday(&tstop_cpu, NULL);
   tresult.tv_sec = tstop_cpu.tv_sec - tstart_cpu.tv_sec;
   tresult.tv_usec = tstop_cpu.tv_usec - tstart_cpu.tv_usec;
   result = (double)tresult.tv_sec + (double)tresult.tv_usec*1.0e-6;
   return(result);
}


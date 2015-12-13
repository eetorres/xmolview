//--------------------------------------------------------------------
//
//  FILE:  ctimer.h
//
//  Class for time testing
//
//  created: Sat Feb  4 12:37:18 MST 2012
//  updated:
//
//  Copyrigth 2002 by Edmanuel Torres
//                    eetorres@gmail.com
//--------------------------------------------------------------------

// Insert a declaration for system_time here.

#ifndef _CTIME_H_
#define _CTIME_H_

#include <ctime>
#include <iostream>

typedef clock_t system_time;
//typedef time_t system_time;

class CTimer
{
  public:
    // Start and stop the timer
    void start(void){ start_time=clock();};
    //void start(void){ time(&start_time);};
    void stop(void){  stop_time=clock();};
    //void stop(void){  time(&stop_time);};

    // Compute the elapsed time (in seconds)
    double get_elapsed_time () const {
      return  difftime(stop_time,start_time )/CLOCKS_PER_SEC;
      //return  (double)difftime(stop_time,start_time );
    };
    // Display the elapsed time (in seconds)
    void show() const {
      double t =  difftime(stop_time,start_time )/CLOCKS_PER_SEC;
      std::cout<<" elapsed time = "<<t<<" s"<<std::endl;
      //return  (double)difftime(stop_time,start_time );
    };

  private:
    system_time start_time,   // Time that the timer was started
                stop_time;    // Time that the timer was stopped
};

#endif

//END
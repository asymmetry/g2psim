// this is just a header file to turn on or off 
// the debug information

#ifndef _GLOBALDEBUGER_H_
#define _GLOBALDEBUGER_H_ 1

#include <iostream>
#include <string>
#include <stdio.h>

#define STOP4DEBUG Stop4Debug();

extern int  Global_Debug_Level;
extern int  iKeepGoing;
extern int  Global_Skip_Counter;

template <class T>
void ECHO(T i)
{
    std::cout<<">>>debug position "<<i<<std::endl;
}

template <class T>
void ECHO(const char *name, T value)
{
    std::cout<<">>> "<<name<<" = "<<value<<" <<<"<<std::endl;
}


template <class T>
void SetGlobalDebugLevel(const char *caller, T val)
{
    Global_Debug_Level=int(val);
    std::cout<<">>> "<<caller<<" set Global_Debug_Level to "<<Global_Debug_Level<<" <<<"<<std::endl;
}

extern int Stop4Debug(int iNewEvent=0);
extern void PrintDebugMenu();

#endif //_GLOBALDEBUGER_H_

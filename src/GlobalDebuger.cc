// this is just a header file to turn on or off 
// the debug information
// #define ECHO(i) printf(">>>debug position %g\n",(float)(i));

#include "GlobalDebuger.hh"
#include <stdlib.h>
using namespace std;

int  Global_Debug_Level=1;
int  iKeepGoing=0;	//-1 exit; 0 false; 1 true
int  Global_Skip_Counter=0;

void PrintDebugMenu()
{  
	cout<<endl;
	cout<<"==================================================="<<endl;
	cout<<"GlobalDebug: The help menu of global debug system is:"<<endl;
	cout<<"\t h|H: print help menu and current status"<<endl
		<<"\t l|L #: change the debug level, if the level==-999, exit immediately"<<endl
		<<"\t j|J #: process # of events without stopping or showing this msg"<<endl
		<<"\t k|K: keep running without showing this msg till the end"<<endl
		<<"\t q|Q: quit now"<<endl;
	cout<<"---------------------------------------------------"<<endl;
	cout<<"Current debug information: Global_Debug_Level="<<Global_Debug_Level
		<<" iKeepGoing="<<((iKeepGoing)? "true" : "false")<<endl;
	cout<<"==================================================="<<endl;
	cout<<endl;
}	

int Stop4Debug(int iNewEvent)
{
	if(Global_Skip_Counter>0) {Global_Skip_Counter--;return 0;}
	if(!iKeepGoing) 
	{ 
		char foo;
		if(iNewEvent) cout<<"=======================Event End  ============================"<<endl;
		else cout<<"=======================DEBUG BLOCK End  ============================"<<endl;

		cout<<"GlobalDebug: Press k to skip asking, L# change debug level, any other key to continue:";
		scanf("%c",&foo);
		if(foo=='k' || foo=='K')  iKeepGoing=1;
		else if(foo=='h' || foo=='H')  PrintDebugMenu();
		else if(foo=='l' || foo=='L')  
		{
			scanf("%d",&Global_Debug_Level);
			if(Global_Debug_Level==-999) exit(Global_Debug_Level);
			cout<<"==================================================="<<endl;
			cout<<"GlobalDebug: Changing global debug level to "<<Global_Debug_Level<<endl;
			cout<<"==================================================="<<endl;
		}
		else if(foo=='q' || foo=='Q')  
		{
			iKeepGoing=-1;
		}

		if(iNewEvent) cout<<"=======================Event Start ============================"<<endl;
		else cout<<"=======================DEBUG BLOCK Start============================"<<endl;
	}
	return iKeepGoing;
}


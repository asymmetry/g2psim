#include<stdlib.h>
#include<stdio.h>
#include<math.h>

extern "C" {
void epc_(int *, int *, int *, double *, double *, double *, double *);
}

void photon(int tgt_Z,int tgt_N,double E1,double PTP,double THP,double *xs)
{
  int i,j;
  double p0,p1,p2;
  int thbin,Z,N,PART;
  double th,sinth,sinthmax,sinthrange,domega;
  double deltap,sig;
  double tempsig,temprate,theta,E;
  double mpi0,pi;
  double a,b,c;
  
  mpi0=134.9766;
  pi=3.141592654;
  thbin=2000;
  
  E=E1;Z=tgt_Z;N=tgt_N;PART=111;
  p0=PTP;
  if(mpi0/2.0/p0<1)
    sinthmax=mpi0/2/p0;
  else
    sinthmax=1;
  sinthrange=sinthmax/thbin;
  tempsig=0;temprate=0;
    
  sinth=sinthrange/2.0;
  th=asin(sinth);
  a=mpi0*mpi0-4*p0*p0*sin(th)*sin(th);
  b=mpi0*mpi0/(2*p0*sin(th)*sin(th));
  c=mpi0*cos(th)/(2*p0*sin(th)*sin(th));
  p1=sqrt((b+c*sqrt(a))*(b+c*sqrt(a))-mpi0*mpi0);
  a=mpi0*mpi0-4*(p0-0.1)*(p0-0.1)*sin(th)*sin(th);
  b=mpi0*mpi0/(2*(p0-0.1)*sin(th)*sin(th));
  c=mpi0*cos(th)/(2*(p0-0.1)*sin(th)*sin(th));
  if((b+c*sqrt(a))*(b+c*sqrt(a))-mpi0*mpi0<0)
    p2=p1;
  else
    p2=sqrt((b+c*sqrt(a))*(b+c*sqrt(a))-mpi0*mpi0);
  deltap=p2-p1;

  theta=0+THP;
  epc_(&PART,&Z,&N,&E,&p0,&theta,&tempsig);
  domega=asin(sinth+sinthrange/2.0)*asin(sinth+sinthrange/2.0);
  temprate=temprate+tempsig*domega*deltap;
  //for(i=1;i<thbin;i++)
  for(i=1837;i<1838;i++)
  {
    sinth=sinthrange*(i+0.5);
    th=asin(sinth);
    a=mpi0*mpi0-4*p0*p0*sin(th)*sin(th);
    b=mpi0*mpi0/(2*p0*sin(th)*sin(th));
    c=mpi0*cos(th)/(2*p0*sin(th)*sin(th));
    p1=sqrt((b+c*sqrt(a))*(b+c*sqrt(a))-mpi0*mpi0);
    a=mpi0*mpi0-4*(p0-0.1)*(p0-0.1)*sin(th)*sin(th);
    b=mpi0*mpi0/(2*(p0-0.1)*sin(th)*sin(th));
    c=mpi0*cos(th)/(2*(p0-0.1)*sin(th)*sin(th));
    if((b-c*sqrt(a))*(b-c*sqrt(a))-mpi0*mpi0<0)
      p2=p1;
    else
      p2=sqrt((b+c*sqrt(a))*(b+c*sqrt(a))-mpi0*mpi0);
    deltap=p2-p1;

//    printf("%lf\t%lf\t%lf\t%lf\t%lf\n",th,b,c,p2,deltap);
    domega=(asin(sinth+sinthrange/2.0)*asin(sinth+sinthrange/2.0)- asin(sinth-sinthrange/2.0)*asin(sinth-sinthrange/2.0))/2.0;
    theta=THP+th*180.0/pi;
    epc_(&PART,&Z,&N,&E,&p1,&theta,&tempsig);
    temprate=temprate+tempsig*domega*deltap;
//    printf("%d\t%lf\t%lf\t%lf\n",i,tempsig,domega,deltap);
    theta=THP-th*180.0/pi;
    epc_(&PART,&Z,&N,&E,&p1,&theta,&tempsig);
    temprate=temprate+tempsig*domega*deltap;
//    printf("%d\t%lf\t%lf\t%lf\n",i,THP,theta,th*180.0/pi);
  }
  sig=temprate/(2*pi)/0.1;
  *xs=sig;
}
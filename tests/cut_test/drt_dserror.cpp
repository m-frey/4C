/*!----------------------------------------------------------------------
\file drt_dserror.cpp

\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/

#include "drt_dserror.H"
#include <mpi.h>

#ifdef THROWELEMENTERRORS

#include "drt_globalproblem.H"
#include "../drt_io/io_control.H"


static bool elementcall = false;
static int elementerrorcount = 0;


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void EnterElementLoop()
{
  if (elementerrorcount!=0)
    dserror("inconsistent error system: %d", elementerrorcount);
  elementcall = true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ExitElementLoop()
{
  elementcall = false;
  if (elementerrorcount!=0)
    dserror("%d element errors occured",elementerrorcount);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ElementError(int ele, const std::string& err)
{
  if (not elementcall)
    dserror("element error outside of element loop");
  elementerrorcount += 1;
  fprintf(DRT::Problem::Instance()->ErrorFile()->Handle(),"element %d error: %s\n",ele,err.c_str());
  fflush(DRT::Problem::Instance()->ErrorFile()->Handle());
}


#endif

static int         latest_line = -1;
static std::string latest_file = "{dserror_func call without prototype}";
/*----------------------------------------------------------------------*
 |  assert function                                          mwgee 11/06|
 | used by macro dsassert in dserror.H                                  |
 *----------------------------------------------------------------------*/
void cpp_dsassert_func(const std::string file, const int line, const bool test, const std::string text)
{
#ifdef DEBUG
  if (!test)
  {
    latest_file = file;
    latest_line = line;
    cpp_dserror_func(text);
  }
#endif
  return;
} /* end of dsassert_func */

/*----------------------------------------------------------------------*
 |  assert function                                          mwgee 11/06|
 | used by macro dsassert in dserror.H                                  |
 *----------------------------------------------------------------------*/
extern "C"
void cpp_dsassert_func(const char* file, const int line, const int test, const char* text)
{
#ifdef DEBUG
  if (!test)
  {
    latest_file = file;
    latest_line = line;
    cpp_dserror_func(text);
  }
#endif
  return;
} /* end of dsassert_func */

/*----------------------------------------------------------------------*
 |  set file and line                                        mwgee 11/06|
 | used by macro dsassert and dserror in dserror.H                      |
 *----------------------------------------------------------------------*/
void cpp_dslatest(const std::string file, const int line)
{
  latest_file = file;
  latest_line = line;
}

/*----------------------------------------------------------------------*
 |  set file and line                                        mwgee 11/06|
 | used by macro dsassert and dserror in dserror.H                      |
 *----------------------------------------------------------------------*/
extern "C"
void cpp_dslatest(const char* file, const int line)
{
  latest_file = file;
  latest_line = line;
}

/*----------------------------------------------------------------------*
 | error function                                            mwgee 11/06|
 | used by macro dsassert in dserror.H                                  |
 *----------------------------------------------------------------------*/
extern "C"
void cpp_dserror_func(const char* text, ...)
{
  int myrank=0;

  char errbuf[2048];

  va_list ap;
  va_start(ap, text);

  sprintf(errbuf,"PROC %d ERROR in %s, line %i:\n",myrank,latest_file.c_str(),latest_line);
  vsprintf(&errbuf[strlen(errbuf)],text,ap);

  va_end(ap);

#ifdef THROWELEMENTERRORS
  if (elementcall)
  {
    throw std::string(errbuf);
  }
  else
#endif
  {
    char line[] = "=========================================================================\n";
    printf("\n\n");
    printf(line);
    printf(errbuf);
    printf("\n");
    printf(line);
    printf("\n\n");
    fflush(stdout);

#ifdef DSERROR_DUMP
    abort();
#endif

#ifdef PARALLEL
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
#else
    exit(1);
#endif
  }
} /* end of dserror_func */

/*----------------------------------------------------------------------*
 | error function                                            mwgee 11/06|
 | used by macro dsassert in dserror.H                                  |
 *----------------------------------------------------------------------*/
void cpp_dserror_func(const std::string text, ...)
{
  int myrank=0;

  char errbuf[2048];

  va_list ap;
  va_start(ap, text);

  sprintf(errbuf,"PROC %d ERROR in %s, line %i:\n",myrank,latest_file.c_str(),latest_line);
  vsprintf(&errbuf[strlen(errbuf)],text.c_str(),ap);

  va_end(ap);

#ifdef THROWELEMENTERRORS
  if (elementcall)
  {
    throw std::string(errbuf);
  }
  else
#endif
  {
    char line[] = "=========================================================================\n";
    printf("\n\n");
    printf(line);
    printf(errbuf);
    printf("\n");
    printf(line);
    printf("\n\n");
    fflush(stdout);

#ifdef DSERROR_DUMP
    abort();
#endif

#ifdef PARALLEL
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
#else
    exit(0);
#endif
  }
} /* end of dserror_func */

/*!---------------------------------------------------------------------
\file
\brief file reading routines

---------------------------------------------------------------------*/

#include "../headers/standardtypes.h"
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;                      


/*! 
\addtogroup FRSYSTEM 
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in input_control_global.c
and the type is in standardtypes.h                                                  
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>

*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*!---------------------------------------------------------------------
\brief init file reading system                                              

<pre>                                                        m.gee 8/00 
The routine open the inputfile allfiles.input on proc 0 and
reads it twice:
-Only counting lines, which are no comment (comment: //)
-Copying file to allfiles.input_file without comments
-Broadcasting allfiles.input_file to all procs
-Sets the fr-systems pointers used by the fr routines on all procs
</pre>
\return void                                                

------------------------------------------------------------------------*/
void frinit()
{
int     i=0;
int     linecount=0;
char   *remarkpointer;
allfiles.numcol=MAXNUMCOL;
if (par.myrank==0)
{
/*--------------------make a copy of the input file without commentars */
   rewind(allfiles.in_input);
/*----------------------------------------------------------read a line */
   while ( feof(allfiles.in_input) == 0 )
   {
      if ( (fgets(allfiles.line,499,allfiles.in_input) == NULL) &&
           (feof(allfiles.in_input) == 0) )
      {
         dserror("An error occured reading a line from input file");
      }
/* ----------------------------------check whether its a commentar line */
      if (strncmp(allfiles.line,"//",2)!=0)
      {
         linecount++; 
      }
   }
}
/*--------------------------------------------broadcast number of lines */
#ifdef PARALLEL 
if (par.nprocs>1)
{
   MPI_Bcast(&linecount,1,MPI_INT,0,MPI_COMM_WORLD);
}
#endif 
allfiles.numrows=linecount;
/*--------------------------------------allocate space for copy of file */
allfiles.input_file_hook=(char*)CALLOC(linecount*(allfiles.numcol),sizeof(char));
if (allfiles.input_file_hook==NULL) dserror("Allocation of memory failed");

allfiles.input_file=(char**)CALLOC(linecount,sizeof(char*));
if (allfiles.input_file==NULL) dserror("Allocation of memory failed");

for (i=0; i<linecount; i++)
{
   allfiles.input_file[i] = &(allfiles.input_file_hook[i*(allfiles.numcol)]);
}
/*-------------------know read the input file into allfiles.input_file */
if (par.myrank==0)
{
   rewind(allfiles.in_input);
   for (i=0; i<linecount; i++)
   {
      read:
      if (fgets(allfiles.input_file[i],499,allfiles.in_input)==NULL)
         dserror("An error occured reading a line from input file");
/*----------------------------------check whether its a commentar line */
      if ((strncmp(allfiles.input_file[i],"//",2)==0))
      {
         goto read;
      }
/*----------------------check whether there is a commentar in the line */
      remarkpointer=strchr(allfiles.input_file[i],'/');
      if (remarkpointer!=NULL)
      {
         strcpy(remarkpointer,"\n\0");
      }
   }
}
/*--------------------------------broadcast the copy of the input file */
#ifdef PARALLEL 
if (par.nprocs>1)
{
      MPI_Bcast(allfiles.input_file_hook,
                allfiles.numrows*allfiles.numcol,
                MPI_CHAR,
                0,
                MPI_COMM_WORLD);
}
#endif
/*---------------------give a copy of the "cleaned" input file on .err */
if (par.myrank==0)
{
   fprintf(allfiles.out_err,"===========================================\n");
   fprintf(allfiles.out_err,"broadcasted copy of input file:            \n");
   fprintf(allfiles.out_err,"===========================================\n");
   for (i=0; i<linecount; i++)
   {
      fprintf(allfiles.out_err,"%s",allfiles.input_file[i]);
   }
   fprintf(allfiles.out_err,"===========================================\n");
   fprintf(allfiles.out_err,"end of broadcasted copy of input file      \n");
   fprintf(allfiles.out_err,"===========================================\n");
   fflush(allfiles.out_err);
/*---------------------------close input file, 'cause no longer needed */
   fclose(allfiles.in_input);
}
/*--------------------------------set fr-system to begin of input_file */
allfiles.actrow=0;
allfiles.actplace=&(allfiles.input_file[0][0]);


return;
} /* end of frinit */




/*!---------------------------------------------------------------------
\brief rewind the copy of the input file                                              

<pre>                                                        m.gee 8/00 
sets pointer allfiles.actplace to allfiles.input_file[0][0]
</pre>
\return void   
                                             
------------------------------------------------------------------------*/
void frrewind()
{
allfiles.actrow=0;
allfiles.actplace=&(allfiles.input_file[0][0]);
return;
} /* end of frrewind */





/*!---------------------------------------------------------------------
\brief find a character string                                             

<pre>                                                        m.gee 8/00 
searches for a given character string in allfiles.input_file
and sets allfiles.actplace to it
If string is not found it terminates the programm with an error message
</pre>
\param string   char[]   (i)   string to search for in input                                
\return void                                               

------------------------------------------------------------------------*/
void frfind(char string[])
{
char message[100];
int  i=0;

#ifdef DEBUG 
dstrc_enter("frfind");
#endif

frrewind();
while ( strstr(allfiles.input_file[i],string) == NULL )
{
   if ( i==allfiles.numrows )
   {
      
      sprintf(message,"frfind:  String %s is not in input file",string);
      dserror(message);
   }
   i++;
}
allfiles.actrow=i;
allfiles.actplace=strstr(allfiles.input_file[i],string);

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of frfind */





/*!---------------------------------------------------------------------
\brief sets a pointer to the next line in the copy of input                                             

<pre>                                                        m.gee 8/00 
sets allfiles.actplace to the next line in allfiles.input_file
</pre>
\return void                                               

------------------------------------------------------------------------*/
void frread()
{

#ifdef DEBUG 
dstrc_enter("frread");
#endif

allfiles.actrow+=1;
if (allfiles.actrow>=allfiles.numrows) dserror("Can't read line, end of input_file reached");
allfiles.actplace=&(allfiles.input_file[allfiles.actrow][0]);

#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of frread */


/*!---------------------------------------------------------------------
\brief reads n integers from input_file                                             

<pre>                                                        m.gee 4/01 
reads n integers from input_file
</pre>
\param string   char[]   (i)   string to search for in actual line
\param var      int*     (o)   adress of field to hold values read    
\param num      int      (i)   number of values to read 
\param ierr     int*     (o)   =0  keyword not found / =1 values read
\return void                                               
\sa frdouble_n() , frint()

------------------------------------------------------------------------*/
void frint_n(char string[],int *var,int num, int *ierr)
{
char message[100];
int  i=0;
char *foundit = NULL;
#ifdef DEBUG 
dstrc_enter("frint");
#endif
/*----------------------------------------------------------------------*/
foundit = strstr(allfiles.input_file[allfiles.actrow],string);
if (foundit != NULL)
{
   foundit += strlen(string);
   for (i=0; i<num; i++)
   {
      foundit=strpbrk(foundit,"-1234567890");
      var[i] = strtol(foundit,&foundit,10);
   }
   *ierr=1;
}
else
{
*ierr=0;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of frint */



/*!---------------------------------------------------------------------
\brief reads an integer from input_file                                            

<pre>                                                        m.gee 8/00 
reads an integer from input_file starting in allfiles.actrow
Searches for the keyword and reads the first thing consisting of -1234567890 behind it
</pre>
\param string   char[]   (i)   string to search for in actual line
\param var      int*     (o)   adress of variable to hold value read    
\param ierr     int*     (o)   =0  keyword not found / =1 values read
\return void                                               
\sa frdouble_n() , frint_n() , frdouble() , frchar()

------------------------------------------------------------------------*/
void frint(char string[],int *var, int *ierr)
{
char message[100];
int  i=0;
char *foundit = NULL;
#ifdef DEBUG 
dstrc_enter("frint");
#endif
/*----------------------------------------------------------------------*/
foundit = strstr(allfiles.input_file[allfiles.actrow],string);
if (foundit != NULL)
{
   foundit += strlen(string);
   foundit=strpbrk(foundit,"-1234567890");
   i=sscanf(foundit," %d ",var);
   if (i!=1 || i==EOF)
   {
      sprintf(message,"frint:  an error occured reading %s",string);
      dserror(message);
   }
   else
   {
   *ierr=1;
   }
}
else
{
*ierr=0;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of frint */


/*!---------------------------------------------------------------------
\brief reads n doubles from input_file                                            

<pre>                                                        m.gee 8/00 
reads n doubles from input_file starting in allfiles.actrow
Searches for the keyword and reads the first thing consisting of -.1234567890 behind it. 
All values to be read must then be continous
</pre>
\param string   char[]   (i)   string to search for in actual line
\param var      double*  (o)   adress of field to hold values read    
\param ierr     int*     (o)   =0  keyword not found / =1 values read
\return void                                               
\sa frdouble() , frint_n() , frint() , frchar()

------------------------------------------------------------------------*/
void frdouble_n(char string[],double *var,int num, int *ierr)
{
int  i;
char *foundit = NULL;
char *end;
#ifdef DEBUG 
dstrc_enter("frdouble");
#endif
/*----------------------------------------------------------------------*/
foundit = strstr(allfiles.input_file[allfiles.actrow],string);
if (foundit != NULL) 
{
   foundit += strlen(string);
   for (i=0; i<num; i++)
   {
      foundit=strpbrk(foundit,"-.1234567890");
      var[i] = strtod(foundit,&foundit);
   }
   *ierr=1;
}
else
{
*ierr=0;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of frdouble */



/*!---------------------------------------------------------------------
\brief reads a double from input_file                                            

<pre>                                                        m.gee 8/00 
reads a double from input_file starting in allfiles.actrow
Searches for the keyword and reads the first thing consisting of -.1234567890 behind it.
</pre>
\param string   char[]   (i)   string to search for in actual line
\param var      double*  (o)   adress of variable to hold value read    
\param ierr     int*     (o)   =0  keyword not found / =1 value read
\return void                                               
\sa frdouble_n() , frint_n() , frint() , frchar()

------------------------------------------------------------------------*/
void frdouble(char string[],double *var, int *ierr)
{
char *foundit = NULL;
char *end;
#ifdef DEBUG 
dstrc_enter("frdouble");
#endif
/*----------------------------------------------------------------------*/
foundit = strstr(allfiles.input_file[allfiles.actrow],string);
if (foundit != NULL) 
{
   foundit += strlen(string);
   foundit=strpbrk(foundit,"-.1234567890");
   *var = strtod(foundit,&end);
   *ierr=1;
}
else
{
*ierr=0;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of frdouble */





/*!---------------------------------------------------------------------
\brief reads a charstring from input file                                            

<pre>                                                        m.gee 8/00 
searches for keyword string and reads continous charstring behind it
user must assure, that the given charpointer space is long enough
to hold the string
</pre>
\param string   char[]   (i)   string to search for in actual line
\param var      char*    (o)   adress of variable to hold string read    
\param ierr     int*     (o)   =0  keyword not found / =1 value read
\return void                                               
\sa frdouble_n() , frint_n() , frint() 

------------------------------------------------------------------------*/
void frchar(char string[],char *var, int *ierr)
{
int  i=0;
char *foundit = NULL;
char *end;

#ifdef DEBUG 
dstrc_enter("frchar");
#endif

foundit = strstr(allfiles.input_file[allfiles.actrow],string);
if (foundit != NULL) 
{
   foundit += strlen(string);
   sscanf(foundit," %s ",var);
   *ierr=1;
}
else
{
*ierr=0;
}

#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of frchar */



/*!---------------------------------------------------------------------
\brief checks for a keyword in actual line                                            

<pre>                                                        m.gee 8/00 
checks for a keyword in actual line of input_file
</pre>
\param string   char[]   (i)   string to search for in actual line
\param ierr     int*     (o)   =0  keyword not found / =1 value read
\return void                                               
\sa frdouble_n() , frint_n() , frint() , frchar()

------------------------------------------------------------------------*/
void frchk(char string[], int *ierr)
{
int  i=0;
char *foundit = NULL;
char *end;

foundit = strstr(allfiles.input_file[allfiles.actrow],string);
if (foundit != NULL) *ierr=1;
else *ierr=0;
return;
} /* end of frchk */



/*!---------------------------------------------------------------------
\brief close and delete input file copy                                            

<pre>                                                        m.gee 4/01 
close and delete input file copy. All memory allocated to the
fr-system is freed
</pre>
\warning Nothing can be read after call to this routine
\return void                                               
\sa frdouble_n() , frint_n() , frint() , frchar() , frchk()

------------------------------------------------------------------------*/
void frend()
{
/*----------------------------------------------------------------------*/
FREE(allfiles.input_file);
FREE(allfiles.input_file_hook);
allfiles.input_file=NULL;
allfiles.input_file_hook=NULL;
allfiles.actplace = NULL;
allfiles.numrows=0;
allfiles.numcol=0;
allfiles.actrow=0;
/*----------------------------------------------------------------------*/
return;
} /* end of frend */

/*! @} (documentation module close)*/

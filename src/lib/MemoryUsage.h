#ifndef _MemoryUsage_included_
#define _MemoryUsage_included_

# include <unistd.h>
/*
void printMemoryUsage()
{
        int who=RUSAGE_SELF;
        struct rusage usage;
        int ret;

        ret=getrusage(who, &usage);
        cout <<endl;
        //cout <<"Memory usage= "<<usage.ru_utime<<endl;
        //cout <<"Memory usage= "<<usage.ru_stime<<endl;
        cout <<"Memory usage= "<<usage.ru_minflt<<endl;
        cout <<"Memory usage= "<<usage.ru_majflt<<endl;
        cout <<"Memory usage= "<<usage.ru_nswap<<endl;
        cout <<endl;
}
*/
void printMemoryUsage()
{
        pid_t pid=getpid();
        string mem, procFile;
        vector<string> lines;

        procFile="/proc/"+IntToStr(pid)+"/stat";
        ReadLines(procFile, lines, "procFile");
        mem=GetWord2(lines[0], 23);
        cout <<"Memory Usage= "<<mem<<endl;
}
#endif

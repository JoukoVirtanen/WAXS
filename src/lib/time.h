#ifndef _time_included_
#define _time_included_

Real calcTimeDiff(timeval start, timeval end)
{
        long seconds, nseconds;

        seconds=end.tv_sec-start.tv_sec;
        nseconds=end.tv_usec-start.tv_usec;
        return seconds+nseconds/1000000.0;
}

#endif

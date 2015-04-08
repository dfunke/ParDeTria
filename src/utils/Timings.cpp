#include "Timings.h"

#include <unistd.h>
#include <time.h>

#include "version.h"

ExperimentRun::ExperimentRun() {

    // output datetime, hostname and git commit

    char datetime[64];
    time_t tnow = time(NULL);
    strftime(datetime,sizeof(datetime),"%Y-%m-%d %H:%M:%S", localtime(&tnow));
    addTrait("datetime", datetime);

    char hostname[128];
    gethostname(hostname, sizeof(hostname));
    addTrait("host", hostname);

    addTrait("git-rev", GIT_COMMIT);

}


# Changing the source code

In order to perform more customised simulations, the source code of RAMSES can be changed/extended.

There are two different ways in which this task can be managed:

1. Create an appropriately called git branch and all the development changes will be reflected in the code commits. Great variety of resources can be found on this, but see e.g. [this brief introduction](https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging). This is the **recommended** way of patching RAMSES.
2. Create a patch directory and store there your changed files. This is presently **disfavoured**, but will be discussed here due to historical reasons.

# Patch directory
After storing all the changed source files in the directory of choice you have to pass this path to Makefile - adjust the `PATCH` variable (make sure there is no trailing space!).

During the compilation, `ramses/utils/scripts/cr_write_patch.sh` will be invoked. This script prepares a Fortran source code file which contains all the patches as a long string. This is written in each RAMSES output directory for later reference (`patches.txt`).

**WARNING**: Changes in RAMSES source code that are uncommitted or not in patch directory are not being tracked!
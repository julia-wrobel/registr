## Test environments
* local R installation, R 3.6.2
* ubuntu 16.04 (on travis-ci), R 3.6.2
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.


## CRAN submission 1.0

Documentation fixes were requested in an email received by Jelena Saf on March 6, 2020. Bullets below detail changes that were made to address these requests.

* Year was added to references in Description field of DESCRIPTION file
* Description field of DESCRIPTION file was lengthened to provide a full paragraph of information
* Author field of DESCRIPTION file was changed to include Erin McDonnell and Jihui Lee was removed (Jihui was accidentally referenced as an author in another file but was not actually a contributor to the package)
* Removed \dontrun{} from examples that can be checked relatively quickly
* Replaced \dontrun{} with \donttest{} for slow examples. The main user-facing functions all have examples that are tested, the \donttest{} examples just provide additional documentation for user on how to use the software for a large real-world data set.
* Have adjusted .Rd files to comply with CRAN standards
* Have added \value section to all .Rd files that are not data files to explain the function results
  



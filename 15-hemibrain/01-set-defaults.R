# This script assumed that you have run the file "15-hemibrain/00-setup.R"

## First, we need to set up some default, like which dataset you want to look at.
if (!require("usethis")) install.packages("usethis")

## To make life easier, you can then edit your R.environ file to contain
## information about the neuPrint server you want to speak with, your token
## and the dataset hosted by that server, that you want to read. For
## detailed instructions use:
?neuprint_login

## In brief:
### You will very likely want to set the following environment variables in
### your .Renviron file. This
### file is read by R on startup. In this way the neuprintr package will
### automatically login to your preferred neuPrint server. Note that
### environment variables will also be inherited by child R sessions. This
### means for example that they will be available when running knitr reports,
### tests or R CMD Check from RStudio. In order to edit your R.profile or
### R.environ files easily and directly, try using
### usethis::edit_r_environ() and usethis::edit_r_profile()
### Copy and paste the lines below there to access the HemiBrain data.
### You will need to use your own token, once you sign in to neuPrint using a google account.
### The one given is a non-functional example.
neuprint_server = "https://neuprint.janelia.org"
neuprint_token = "asBatEsiOIJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJlbWFpbCI6ImFsZXhhbmRlci5zaGFrZWVsLmJhdGVzQGdtYWlsLmNvbSIsImxldmVsIjoicmVhZHdyaXRlIiwiaW1hZ2UtdXJsIjoiaHR0cHM7Ly9saDQuZ29vZ2xldXNlcmNvbnRlbnQuY29tLy1QeFVrTFZtbHdmcy9BQUFBQUFBQUFBDD9BQUFBQUFBQUFBQS9BQ0hpM3JleFZMeEI4Nl9FT1asb0dyMnV0QjJBcFJSZlBRL21vL3Bob3RvLapwZz9zej01MCIsImV4cCI6MTczMjc1MjU2HH0.jhh1nMDBPl5A1HYKcszXM518NZeAhZG9jKy3hzVOWEU"
neuprint_dataset = "hemibrain:v1.0"
usethis::edit_r_environ()

## Once you have updated your R environ file, you will need to restart the session.

## Then see if a simple function works for you:
available.datasets = neuprint_datasets()
available.datasets
### Did you get something??

## Great! Let's see how you are connecting:
conn = neuprint_login()
### After successful login, the neuprint_connection object will
### contain a cookie field that includes a sessionid that is required
### for subsequent GET/POST operations using the package httr. When
### Cache=TRUE (the default) the open connection object is cached and
### will be used when EITHER neuprint_login is called with enough
### information to indicate that the same server is desired OR (when no
### information about the server is passed to neuprint_login) the last
### opened connection will be used. A new connection can be made using
### Force = TRUE, which is advisable as a first call for debugging if
### you are having issues querying the server.


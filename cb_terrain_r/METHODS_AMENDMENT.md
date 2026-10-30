# Proposed STAR Methods amendment (author review)

The analyses were independently reproduced in R (version 4.3.3), using jsonlite
(version 1.8.8) for metadata import. Linear models and conditional logistic
models were implemented using explicit design matrices and cluster-robust
covariance calculations to match the archived Python specifications. Numerical
agreement was assessed against archived Python outputs with an absolute
tolerance of 10^-6. Circular-shift analyses used identical archived lag schedules
in both implementations; wild-cluster tests enumerated all 128 sign assignments
for seven bats. Undefined cluster-robust uncertainty for single-bat fits was
reported as missing. The code, environment information and validation report are
available at [INSERT IMMUTABLE CODE URL FROM PUBLICATION_STATUS.md]. The software
archive DOI and data accession remain to be assigned/confirmed.

Before insertion: retain the existing Python package/version information because
Python was used for the original analysis. Confirm author names, license, data
availability, original terrain processing and device time-zone metadata. Do not
state that raw-DEM preprocessing or final map rendering has been reproduced in R.

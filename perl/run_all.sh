#!/bin/sh

perl pg_01_add_nsamples.pl
perl pg_02_add_nreads.pl
perl pg_03_remove_variables.pl
perl pg_04_remove_samples.pl

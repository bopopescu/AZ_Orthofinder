#!/bin/bash
module load blast
module load mcl
module unload python; module load python/64_2.7.3
module unload perl; module load perl
export PERL5LIB=/opt/az/perlfoundation/perl/5.16/lib/site_perl/5.16.0/x86_64-linux-thread-multi-ld/:$PERL5LIB
setenv PERL5LIB /opt/az/perlfoundation/perl/5.16/lib/site_perl/5.16.0/x86_64-linux-thread-multi-ld/:$PERL5LIB

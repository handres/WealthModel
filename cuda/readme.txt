COMPILATION

  /usr/local/cuda/bin/nvcc -O3 -o wealthgpu wealthgpu.cu

EXECUTION

 ./wealthgpu npersons niterations initwealth deltawealth ndw beta

As you can guess,

 ./wealthgpu 100000 100 1.0 0.01 1024 1.36

will run 100,000 persons for 100 iterations
with each person’s initial wealth equal to 1.0 and
delta wealth packet per drop of 0.01,
1024 packet drops per iteration, and beta exponent of 1.36.

If initwealth is specified as a negative value, say, -10.0,
then its negative (i.e., the absolute value) will be used as
the mean of uniformly randomly distributed initial wealths among persons.

ENVIRONMENT VARIABLES

There are several environment variables that can control the behavior.

 * DEBUG
   export DEBUG=d will set debug level to the specified integer level 'd'.

   Specify DEBUG=0 for the fastest runs.
   Specify larger integers for increasingly detailed output.

 * VERIFY
   export VERIFY=TRUE|FALSE will turn on/off error checks at runtime.

   Specify VERIFY=FALSE for the fastest runs.

 * WRITERESULTS

   export WRITERESULTS=TRUE|FALSE will turn on/off periodic writing of wealths to files (see also NWRITES).

 * NWRITES

   export NWRITES=w will write (across w equally spaced intervals) the
   individual wealths into binary files timestamped by number of iterations.
   Obviously, 1<=p<=niterations.

 * NPRINTS

   export NPRINTS=p will print p status updates of execution.
   Obviously, 1<=p<=niterations.

 * NORMALIZE

   export NORMALIZE=TRUE|FALSE will normalize exponentiated wealths to [0,1], or use actual values.

Kalyan Perumalla - July 20, 2016

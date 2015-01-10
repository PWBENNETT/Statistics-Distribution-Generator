Statistics-Distribution-Generator
=================================

A means of composing random number generators together using simple operators

Install using the CPAN interface of your choice, or download from CPAN, and then

    perl Makefile.PL
    make
    make install

If you want to risk random breakage (this is a random number generator after
all), you may run `make test` before `make install`. The test suite is
imperfect, however, precisely because it's hard to predict what random numbers
will do, so you _may_ run into test "failures" that are merely the result of
something very unlikely happening, rather than something being broken.

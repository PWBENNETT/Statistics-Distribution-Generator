#!perl -T

use Test::More tests => 1;

BEGIN {
    use_ok( 'Statistics::Distribution::Generator' );
}

diag( "Testing Statistics::Distribution::Generator $Statistics::Distribution::Generator::VERSION, Perl $], $^X" );

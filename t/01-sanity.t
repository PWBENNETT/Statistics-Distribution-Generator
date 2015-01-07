#!perl -T

use Test::More tests => 5;

use Statistics::Descriptive;
use Statistics::Distribution::Generator qw( :all );

my $num_tests = $ENV{ SDG_TESTS } || 100_000;
my $accuracy = defined($ENV{ SDG_ACC }) ? $ENV{ SDG_ACC } : 0.01;

diag('Set $ENV{ SDG_TESTS } to control the number of iterations inside each sanity test') unless defined $ENV{ SDG_TESTS };
diag('Set $ENV{ SDG_ACC } to control the desired accuracy of each sanity test') unless defined $ENV{ SDG_ACC };

{
    my $gaussian = gaussian(0, 1);
    my $s = Statistics::Descriptive::Sparse->new();
    for (1 .. $num_tests) {
        $s->add_data($gaussian);
    }
    ok(
        $s->standard_deviation > (1 - $accuracy)
        && $s->standard_deviation < (1 + $accuracy)
    );
    ok(
        $s->mean > -$accuracy
        && $s->mean < $accuracy
    );
}

{
    my $cointoss = supplied(0) | supplied(1);
    my $s = Statistics::Descriptive::Sparse->new();
    for (1 .. $num_tests) {
        $s->add_data($cointoss);
    }
    ok(
        $s->mean > (0.5 - $accuracy)
        && $s->mean < (0.5 + $accuracy)
    );
}

{
    my $coinA = supplied(0) | supplied(1);
    my $coinB = supplied(0) | supplied(1);
    my $twocoins = $coinA x $coinB;
    my $s = Statistics::Descriptive::Sparse->new();
    for (1 .. $num_tests) {
        $s->add_data($_) for @$twocoins;
    }
    ok(
        $s->mean > (0.5 - $accuracy)
        && $s->mean < (0.5 + $accuracy)
    );
}

{
    my $d3 = supplied(1) | supplied(2) | supplied(3);
    my $s = Statistics::Descriptive::Sparse->new();
    for (1 .. $num_tests) {
        $s->add_data($d3);
    }
    ok(
        $s->mean > (2 - $accuracy)
        && $s->mean < (2 + $accuracy)
    );
}


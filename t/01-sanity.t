#!perl -T

use Test::More tests => 3;

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

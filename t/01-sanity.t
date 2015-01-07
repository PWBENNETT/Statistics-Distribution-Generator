#!perl -T

use Test::More tests => 2;

use Statistics::Descriptive;
use Statistics::Distribution::Generator qw( :all );

{
    my $g = gaussian(0, 1);
    my $s = Statistics::Descriptive::Full->new();
    for (1 .. 10000) {
        $s->add_data($g);
    }
    ok(
        $s->standard_deviation > 0.5
        && $s->standard_deviation < 1.5
    );
    ok(
        $s->mean > -0.5
        && $s->mean < 0.5
    );
}

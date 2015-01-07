#!perl -T

use Test::More tests => 1;

use Statistics::Descriptive;
use Statistics::Distribution::Generator qw( :all );

{
    my $g = gaussian(0, 1);
    my $s = Statistics::Descriptive::Full->new();
    for (1 .. 1000) {
        $s->add_data($g + 0);
    }
    ok('All good');
}

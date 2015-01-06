package Statistics::Distribution::Generator;

use 5.018;
use utf8;
use overload (
    '0+' => '_render',
    '""' => '_render',
    '|' => '_add_alternative',
    'x' => '_add_dimension',
    fallback => 0,
);

use List::AllUtils qw( reduce );
use Exporter qw( import );

sub gaussian ($$);
sub uniform ($$);
sub logistic ();
sub supplied ($);
sub gamma ($$);

our $VERSION = '0.001';

our @EXPORT_OK = qw( gaussian uniform logistic supplied gamma );
our %EXPORT_TAGS = (':all' => \@EXPORT_OK);

our $pi = 3.14159265358979323846264338327950288419716939937510;
our $two_pi = 2 * $pi;
our $e = exp 1;

sub _render {
    my $self = shift;
    if ($self->{ alts }) {
        my $accum = reduce { $a + $b } map { $_->{ weight } // 1 } @{$self->{ alts }};
        my $n = rand() * $accum;
        my $alt;
        for $alt (@{$self->{ alts }}) {
            if ($n >= ($alt->{ weight } // 1)) {
                $n -= $alt_weight;
                last if $n <= 0;
            }
        }
        return $alt->_render;
    }
    elsif ($self->{ dims }) {
        my @rv;
        for my $dimension (@{$self->{ dims }}) {
            push @rv, $dimension->_render;
        }
        return \@rv;
    }
    else {
        die "Something horrible has happened";
    }
}

sub gaussian ($$) {
    my ($mean, $sigma) = @_;
    $mean //= 0;
    $sigma //= 1;
    return bless { mean => $mean, sigma => $sigma }, 'Statistics::Distribution::Generator::gaussian';
}

sub uniform ($$) {
    my ($min, $max) = @_;
    $min //= 0;
    $max //= 1;
    return bless { min => $min, max => $max }, 'Statistics::Distribution::Generator::uniform';
}

sub logistic () {
    return bless { }, 'Statistics::Distribution::Generator::logistic';
}

sub supplied ($) {
    my ($iv) = @_;
    my $rv;
    if (ref $iv eq 'CODE') {
        $rv = { code => $iv };
    }
    else {
        $rv = { code => sub { return $iv } };
    }
    return bless $rv, 'Statistics::Distribution::Generator::supplied';
}

sub gamma ($$) {
    my ($order, $scale) = map { $_ // 1 } @_;
    return bless {
        order => $order,
        scale => $scale,
        norder => int($order),
    }, 'Statistics::Distribution::Generator::gamma';
}

sub _rand_nonzero {
    my $rv;
    1 while (!$rv = rand);
    return $rv;
}

sub _gamma_int {
    my $order = shift;
    if ($order < 12){
        my $prod = 1;
        for (my $i=0; $i<$order; $i++){
            $prod *= _rand_nonzero();
        }
        return -log($prod);
    }
    else {
        return _gamma_large_int($order);
    }
}

sub _tan { sin($_[0]) / cos($_[0]); }

sub _gamma_large_int {
    my $order = shift;
    my $sqrt = sqrt(2 * $order - 1);
    my ($x,$y,$v);
    do {
        do {
            $y = _tan($pi * rand);
            $x = $sqrt * $y + $order - 1;
        } while ($x <= 0);
        $v = rand;
    } while ($v > (1 + $y * $y) * exp(($order - 1) * log($x / ($order - 1)) - $sqrt * $y));
    return $x;
}

sub _gamma_frac {
    my $order = shift;
    my $p = $e / ($order + $e);
    my ($q, $x, $u, $v);
    do {
        $u = rand;
        $v = _rand_nonzero();
        if ($u < $p){
            $x = exp((1 / $order) * log($v));
            $q = exp(-$x);
        }
        else {
            $x = 1 - log($v);
            $q = exp(($order - 1) * log($x));
        }
    } while (rand >= $q);
    return $x;
}

sub _add_alternative {
    my ($lhs, $rhs, $swapped) = @_;
    ($lhs, $rhs) = ($rhs, $lhs) if $swapped;
    $rhs = supplied $rhs unless ref($rhs) =~ /^Statistics::Distribution::Generator/;
    my $self
        = ref($lhs) eq 'Statistics::Distribution::Generator'
        ? { %$lhs }
        : { alts => [ $lhs ] }
        ;
    bless $self, 'Statistics::Distribution::Generator';
    push @{$self->{ alts }}, $rhs;
    return $self;
}

sub _add_dimension {
    my ($lhs, $rhs, $swapped) = @_;
    ($lhs, $rhs) = ($rhs, $lhs) if $swapped;
    $rhs = supplied $rhs unless ref($rhs) =~ /^Statistics::Distribution::Generator/;
    my $self
        = ref($lhs) eq 'Statistics::Distribution::Generator'
        ? { %$lhs }
        : { dims => [ $lhs ] }
        ;
    bless $self, 'Statistics::Distribution::Generator';
    push @{$self->{ dims }}, $rhs;
    return $self;
}

sub Statistics::Distribution::Generator::gaussian::_render {
    my $self = shift;
    my $U = rand;
    my $V = rand;
    return $self->{ mean } + (sqrt(-2 * log $U) * cos($two_pi * $V) * $self->{ sigma });
}

sub Statistics::Distribution::Generator::uniform::_render {
    my $self = shift;
    return ($self->{ max } - $self->{ min }) * rand() + $self->{ min };
}

sub Statistics::Distribution::Generator::logistic::_render {
    my $self = shift;
    return 1 / log rand;
}

sub Statistics::Distribution::Generator::supplied::_render {
    my $self = shift;
    return $self->{ code }->();
}

sub Statistics::Distribution::Generator::gamma::_render {
    my $self = shift;
    my $rv;
    if ($self->{ order } == $self->{ norder }) {
        $rv = $self->{ scale } * _gamma_int($self->{ norder });
    }
    elsif ($self->{ norder } == 0) {
        $rv = $self->{ scale } * _gamma_frac($self->{ order });
    }
    else {
        $rv = $self->{ scale } * (_gamma_int($self->{ norder }) + _gamma_frac($self->{ norder } - $self->{ order }));
    }
    return $rv;
}

1;

__END__

=head1 NAME

Statistics::Distribution::Generator - A way to compose complicated probability functions

=head1 VERSION

Version 0.001

=head1 SYNOPSIS

    use Statistics::Distribution::Generator qw( :all );
    my $g = gaussian 3, 1;
    say $g; # something probably between -3 and 9, but probably about 2 .. 4-ish
    my $cloud = (gaussian 0, 1 x gaussian 0, 1 x gaussian 0, 1);
    say @$cloud; # a 3D vector almost certainly within (+/- 6, +/- 6, +/- 6) and probably within (+/- 3, +/- 3, +/- 3)
    my $combo = (gaussian 100, 15 | uniform 0, 200); # one answer with an equal chance of being picked from either distribution

=head1 DESCRIPTION

This module allows you to bake together multiple "simple" probability distributions into a more complex random number generator.

=head1 EXPORTABLE FUNCTIONS

=over

=item gaussian MEAN, SIGMA

Gaussian Normal Distribution

=back

=over

=item uniform MIN, MAX

A uniform distribution, with equal chance of any n where MIN Z<<>= n Z<<> MAX

=back

=over

=item logistic

Standard Logistic Distribution

=back

=over

=item supplied VALUE

=item supplied CALLBACK

Allows the caller to supply either a constant VALUE which will always be returned, or a coderef CALLBACK that may use any algorithm you like to generate a random number

=back

=item gamma ORDER, SCALE

Gamma Distribution

The distribution function is

    p(x) dx = {1 \over \Gamma(a) b^a} x^{a-1} e^{-x/b} dx
    for x > 0.

=back

=head1 OVERLOADED OPERATORS

=over

=item x

Allows you to compose multi-dimensional random vectors.

    $randvect = $foo x $bar x $baz; # generate a three-dimensional vector

=back

=over

=item |

Allows you to pick a single (optionally weighted) generator from some set of generators.

    $cointoss = supplied 0 | supplied 1; # fair 50:50 result of either 0 or 1

=back

=head1 AUTHOR

The main body of this work is by Paul W Bennett E<paul.w.bennett@gmail.com>

The idea of composing probabilities together comes from a paper by B<TODO: CITE THE PAPER HERE>

The implementation of the Gamma Distribution is by Nigel Wetters Gourlay.

=head1 CAVEATS

Almost no error checking is done. Garbage in I<will> result in garbage out.

=head1 TODO

Build in more probability density functions.

=head1 LICENSE

Artistic 2.0

=cut

package Statistics::Distribution::Generator;

use 5.018;
use utf8;
use overload (
    '0+' => '_render',
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


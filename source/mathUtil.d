import std.typetuple;
import std.range;
import std.traits;
import std.conv : to;
import std.algorithm : map;
import std.range : iota;

public import std.math : PI, SQRT2;
public import core.stdc.math;
public import std.complex;

static import std.math;

//alias stdMath = std.math;
alias stdMath = core.stdc.math;

version(unittest) import std.stdio;

auto linSpace(T = double, T0, T1)(T0 start, T1 stop, size_t nPoints)
{
    T stopT = stop;
    T startT = start;

    T step;
    if(nPoints == 1)
    {
        step = 0; //avoid nan
    }
    else
    {
        step = (stopT - startT)/(nPoints - 1);
    }
    //stopT += step;
    //return iota(dstart, dstop, step).takeExactly(nPoints);
    return iota(nPoints).map!((n) => startT + n*step);
}

auto expi(double y)
{
    return complex(stdMath.cos(y), stdMath.sin(y));
}

auto complex(T0, T1)(T0 r, T1 i)
{
    return std.complex.complex(r, i);
}

auto complex(T)(T c)
{
    static if(is(T : creal))
        return std.complex.complex(c.re, c.im);
    else static if(is(T : Complex!X, X))
        return c;
    else
        return std.complex.complex(c);
}


/+ //ldc can't handle this;
enum bool isStdComplex(T) = is(T == Complex!X, X);
enum bool isBuiltinComplex(T) = is(T == cfloat) || is(T == cdouble) || is(T == creal);
enum bool isComplex(T) = isStdComplex!T || isBuiltinComplex!T;
+/

template isStdComplex(T)
{
    static if(is(T == Complex!X, X))
        enum isStdComplex = true;
    else
        enum isStdComplex = false;
}


template isBuiltinComplex(T)
{
    static if(is(T == cfloat) || is(T == cdouble) || is(T == creal))
        enum isBuiltinComplex = true;
    else
        enum isBuiltinComplex = false;
}

template isComplex(T)
{
    enum isComplex = isStdComplex!T || isBuiltinComplex!T;
}

auto exp(T)(T c)
{
    static if(is(T : Complex!X, X))
        return complex(stdMath.exp(c.re) * expi(c.im));
    else
        return stdMath.exp(c);
}
 
bool approxEqual(T, U, V)(T lhs, U rhs, V maxRelDiff, V maxAbsDiff = 1e-05)
{
    static if(isComplex!T && isComplex!U)
    {
        return std.math.approxEqual(lhs.re, rhs.re, maxRelDiff, maxAbsDiff) &&
            approxEqual(lhs.im, rhs.im, maxRelDiff, maxAbsDiff);
    }
    else
    {
        return std.math.approxEqual(lhs, rhs, maxRelDiff, maxAbsDiff);
    }
}

bool approxEqual(T, U)(T lhs, U rhs)
{
    return .approxEqual(lhs, rhs, 1e-2, 1e-5);
}

auto Re(T)(T v)
    if(isComplex!T)
{
    return v.re;
}
auto Im(T)(T v)
    if(isComplex!T)
{
    return v.im;
}

private alias ComplexMathFunNames = TypeTuple!("cos", "sin", "sqrt");
mixin({
        string res;
        foreach(name; ComplexMathFunNames)
        {
            res ~=
`auto ` ~ name ~ `(T)(T c)
{
    static if(is(T : Complex!X, X))
        return std.complex.` ~ name ~ `(c);
    else
        return stdMath.` ~ name ~ `(c);
}
`;
        }
        return res;
        }()
     );

auto abs(T)(T c)
{
    static if(is(T : Complex!X, X))
        return std.complex.abs(c);
    else
        return std.math.abs(c);
}


private alias BuiltinTypes = TypeTuple!(bool, byte, ubyte, short, ushort, int, uint,
        long, ulong, float, double, real, char, wchar, dchar);

static if(__traits(compiles, mixin(`double(3)`)))
{
    private string builtinConstructCompat()
    {
        string s;
        foreach(T; BuiltinTypes)
        {
            s ~= "alias " ~ T.stringof ~ "_ = " ~ T.stringof ~ ";\n";
        }
        return s;
    }
}
else
{
    pragma(msg, "no uniform construction for builtin types");
    private string builtinConstructCompat()
    {
        string s;
        foreach(T; BuiltinTypes)
        {
            s ~= "auto " ~ T.stringof ~ "_(Q)(Q v) { return v.to!" ~ T.stringof ~ "; }\n";
        }
        return s;
    }
}
mixin(builtinConstructCompat());

auto reshape(R, D ...)(R r, D dims)
    if(isInputRange!r)
{
    return r.take(dims[0]);
}


ulong nextHammingNumber(T)(T targetNum)
    if(isIntegral!T || isFloatingPoint!T)
{
    ulong target;
    static if(isIntegral!T)
        target = targetNum;
    else static if(isFloatingPoint!T)
        target = ceil(targetNum).to!ulong;

    if (target <= 6)
    {
        return target;
    }

    // Quickly check if it's already a power of 2
    if (!(target & (target-1)))
    {
        return target;
    }

    auto match = ulong.max; // Anything found will be smaller
    ulong p5 = 1;
    while (p5 < target)
    {
        ulong p35 = p5;
        while (p35 < target)
        {
            // Ceiling integer division, avoiding conversion to float
            // (quotient = ceil(real(target) / p35))
            auto quotient = iDivCeiling(target, p35);

            // Quickly find next power of 2 >= quotient
            auto p2 = nextPow2(quotient);

            auto N = p2 * p35;
            if (N == target)
                return N;
            else if (N < match)
                match = N;
            p35 *= 3;
            if (p35 == target)
                return p35;
        }
        if (p35 < match)
            match = p35;
        p5 *= 5;
        if (p5 == target)
            return p5;
    }
    if (p5 < match)
        match = p5;
    return match;
}

auto nextPow2(T)(T v)
{
    ulong p2 = 2;
    while(p2 < v)
    {
        p2 <<= 1;
    }
    return p2;
}

unittest
{
    assert(nextPow2(5) == 8, "result " ~  nextPow2(5).to!string ~ " != 8");
    assert(nextPow2(14.5) == 16, "result " ~  nextPow2(14.5).to!string ~ " != 16");
    assert(nextPow2(4) == 4, "result " ~  nextPow2(4).to!string ~ " != 4");
    assert(nextPow2(4UL) == 4, "result " ~  nextPow2(4UL).to!string ~ " != 4");
}

auto nextMultipleOf(T0, T1)(T0 base, T1 multiplier)
{
    if(base % multiplier)
    {
        return base + multiplier - (base % multiplier);
    }
    return base;
}


auto iDivCeiling(ulong a, ulong b)
{
    return 1 + ((a - 1) / b);
}

auto dftFreqs(ulong n, double d = 1.0)
{
    auto step = 1 / (d * n);
    if(n & 1)
    {
        auto res = chain(
                linSpace(0, (n-1)/(2*d*n), (n/2)+1),
                linSpace((1-long_(n))/(2*d*n), -1/(d*n), n/2));
        return res;
    }
    else
    {
        auto res = chain(
                linSpace(0, ((n/2) - 1)/(d * n), n/2),
                linSpace(-0.5/d, -1/(n*d), n/2));
        return res;
    }
}

template ComplexComponentType(T)
{
    static if(is(T == creal))        alias ComplexComponentType = real;
    else static if(is(T == cdouble))      alias ComplexComponentType = double;
    else static if(is(T == cfloat))       alias ComplexComponentType = float;
    else static if(is(T : Complex!X, X)) alias ComplexComponentType = X;
}

unittest
{
    static assert(is(ComplexComponentType!creal == real));
    static assert(is(ComplexComponentType!cdouble == double));
    static assert(is(ComplexComponentType!cfloat == float));
    static assert(is(ComplexComponentType!(Complex!double) == double));
}

template BuiltinComplex(T)
{
    static if(is(T == Complex!X, X))
        alias E = X;
    else static if(is(T == Complex!(X)[], X))
        alias E = X;
    else
        static assert(false, "There is no builtin complex type corresponding to " ~ T.stringof);

    static if(is(X == real))
    {
        alias BuiltinComplex = creal;
    }
    else static if(is(X == double))
    {
        alias BuiltinComplex = cdouble;
    }
    else static if(is(X == float))
    {
        alias BuiltinComplex = cfloat;
    }
}

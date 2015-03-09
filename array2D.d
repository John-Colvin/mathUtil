import std.algorithm : map, joiner;
import std.range : ElementType, chunks, zip, takeExactly,
       hasSlicing, isInputRange, isForwardRange,
       isBidirectionalRange, isRandomAccessRange,
       hasLength, lockstep, stride, iota;
import std.traits : Unqual, isIterable;
import std.conv: text;
import std.stdio : writeln, write;
import std.array : array;

version(unittest)
{
    import std.algorithm : equal, copy;
    import std.range : take, repeat, retro;
}

struct Array2D(T) 
{
    alias E = ElementType!T;
    T data;
    Slice1D[2] slicing;
    size_t nI, nJ;

    this(T data, size_t nI, size_t nJ)
    {
        this.data = data;
        this.nI = nI;
        this.nJ = nJ;
        slicing[0].e = nI;
        slicing[1].e = nJ;
    }

    private this(T data, Slice1D[2] slicing, size_t nI, size_t nJ)
    {
        this.data = data;
        this.slicing = slicing;
        this.nI = nI;
        this.nJ = nJ;
    }

    private auto baseIdx() @property
    {
        return slicing[0].b*nJ + slicing[1].b;
    }

    private auto endIdx() @property
    {
        return slicing[0].e*nJ + slicing[1].e;
    }

/+
    private auto iOrigin() @property
    {
        return data[(i + slicing[0].b)*nJ .. $];
    }

    private auto iRange()
    {
        return data[(i + slicing[0].b)*nJ .. slicing[0].e*nJ]
    }
+/

    //GRRR THIS SHOULD WORK!!
   // alias opSlice(uint i) = Slice1D;

    auto opSlice(size_t i)(size_t b, size_t e)
    {
        return Slice1D(b, e);
    }

    auto opDollar(size_t i)() @property
    {
        return slicing[i].e - slicing[i].b;
    }
/+
    auto opIndex()
    {
        return this;
    }
+/
    auto opIndex(size_t i, size_t j)
    {
        return data[(i + slicing[0].b)*nJ + slicing[1].b + j];
    }

    auto opIndex(size_t i, Slice1D s)
    {
        return data[baseIdx + i*nJ + s.b .. baseIdx + i*nJ + s.e]; 
    }

    auto opIndex(Slice1D s, size_t j)
    {
        return data[baseIdx + s.b*nJ + j .. $][0 .. (s.e - 1 - s.b)*nJ + 1].stride(nJ);
    }

    auto opIndex(Slice1D s0, Slice1D s1)
    {
        Slice1D[2] tmp;
        tmp[0] = s0;
        tmp[1] = s1;
        return typeof(this)(data, tmp, nI, nJ);
    }
    
    void opIndexAssign(Q)(Q rhs)
        if(is(Q : Array2D!P, P))
    in
    {
        assert(rhs.shape == shape);
    }
    body
    {
        size_t i=0;
        debug writeln("opIndexAssign another Array2D");

        auto bE = this.byElement;
        foreach(r; rhs.byElement)
        {
            bE.front = r;
            bE.popFront();
        }
    }

    void opIndexAssign(Q : E)(Q rhs)
    {
        foreach(ref el; this.byElement)
            el = rhs;
    }

    void opIndexAssign(Q : E)(Q v, size_t i, size_t j)
    {
        debug writeln("opIndexAssign single element");
        data[baseIdx + i*nJ + j] = v;
    }

    void opIndexAssign(Q, I0, I1)(Q v, I0 i0, I1 i1)
    {
        debug writeln("opIndexAssign slices");
        auto sliced = this[i0, i1];
        static if(is(typeof(sliced) : Array2D!P, P))
            sliced[] = v;
        else static if(__traits(compiles, {sliced[] = v;}))
            sliced[] = v;
        else static if(__traits(compiles, {sliced[0..$] = v;} ))
            sliced[0..$] = v;
        else static if(isIterable!(typeof(sliced)) && isIterable!Q)
            foreach(ref dest, src; zip(sliced, v))
                dest = src;
        else static if(isIterable!(typeof(sliced)))
            foreach(ref el; sliced)
                el = v;
    }

    auto opUnary(string op)()
    {
        return array2D(data.map!(op~'a'), slicing, nI, nJ);
    }
    
    auto opBinary(string op, Q)(Q rhs)
        if(is(Q : Array2D!P, P))
    in
    {
        assert(rhs.shape == shape);
    }
    body
    {
        auto newData = zip(this.data, rhs.data).map!("a[0]"~op~"a[1]");
        return array2D(newData, slicing, nI, nJ);
    }

    auto byElement()
    out
    {
        static assert(hasSlicing!(Unqual!(typeof(__result))));
    }
    body
    {
        debug writeln("byElement");
        return ByElement(data,
                slicing[0].b,
                slicing[1].b,
                slicing[0].e - slicing[0].b,
                slicing[1].e - slicing[1].b,
                nJ,
                0,
                0);
    }

    struct ByElement
    {
        alias E = ElementType!T;
        T data;
        size_t originI, originJ, sliceNI, sliceNJ, nJ;
        size_t frontOffset, backOffset;

        this(T data, size_t originI, size_t originJ, size_t sliceNI, size_t sliceNJ,
                size_t nJ, size_t frontOffset, size_t backOffset)
        {
            debug writeln("ByElement.this");
            this.data = data;
            this.originI = originI;
            this.originJ = originJ;
            this.sliceNI = sliceNI;
            this.sliceNJ = sliceNJ;
            this.nJ = nJ;
            this.frontOffset = frontOffset;
            this.backOffset = backOffset;
        }

        auto ref front() @property
        in
        {
            assert(!empty);
        }
        body
        {
            debug writeln("front ", originI*nJ + originJ + frontOffset, ' ', data.length);
            return data[originI*nJ + originJ + frontOffset];
        }

        auto ref back() @property
        in
        {
            assert(!empty);
        }
        body
        {
            debug writeln("back");
            return data[(originI+sliceNI)*nJ + originJ + sliceNJ - backOffset];
        }

        private void debugSpew() const
        {
            foreach(i, el; this.tupleof[1..$])
                    write(__traits(allMembers, typeof(this))[i+2], ": ", el,"  ");
            write('\n');
        }

        bool empty() @property
        {
            debug writeln("empty");
            return !(sliceNI > 1
                || (sliceNI == 1 && frontOffset + backOffset < sliceNJ));
        }

        void popFront()
        in
        {
            assert(!empty);
        }
        body
        {
            debug writeln("popFront", frontOffset, ' ', originI, ' ', sliceNI, ' ', empty);
            ++frontOffset;
            if(frontOffset == sliceNJ)
            {
                frontOffset = 0;
                ++originI;
                --sliceNI;
            }
            debug writeln(frontOffset, ' ', originI, ' ', sliceNI, ' ', empty);
        }

        void popBack()
        in
        {
            assert(!empty);
        }
        body
        {
            debug writeln("popBack");
            ++backOffset;
            if(backOffset == sliceNJ)
            {
                backOffset = 0;
                --sliceNI;
            }
        }

        invariant()
        {
            debug(spew) debugSpew();
            assert(frontOffset < sliceNJ, text(frontOffset));
            assert(backOffset < sliceNJ, text(backOffset));
            assert(sliceNJ <= nJ);
        }

        private size_t linearIdx(size_t n)
        {
            return ((n+frontOffset)/sliceNJ + originI)*nJ
                + originJ + (n + frontOffset)%sliceNJ;
        }

        auto ref opIndex(size_t n)
        in
        {
            assert(n < length);
        }
        body
        {
            debug writeln("opIndex");
            return data[linearIdx(n)];
        }

        auto opIndex()
        {
            return this;
        }

        auto opIndex(Slice1D s)
        in
        {
            assert(s.b <= length, 
                    "s.b(" ~ text(s.b) ~ ") > length(" ~ text(length) ~ ")");
            assert(s.e <= length,
                    "s.e(" ~ text(s.e) ~ ") > length(" ~ text(length) ~ ")");
            assert(s.b <= s.e,
                    "s.b(" ~ text(s.b) ~ ") > s.e(" ~ text(s.e) ~ ")");
        }
        body
        {
            debug writeln("opIndex slice ", s);
            auto iShift = (s.b + frontOffset)/sliceNJ;
            auto iBackShift = (sliceNI*sliceNJ - frontOffset - s.e)/sliceNJ;
            return typeof(this)(data,
                    originI + iShift,
                    originJ,
                    sliceNI - iShift - iBackShift,
                    sliceNJ,
                    nJ,
                    (s.b + frontOffset)%sliceNJ,
                    sliceNJ - (s.e + frontOffset - 1)%sliceNJ - 1); 
        }

        auto opSlice(size_t dim)(size_t b, size_t e)
        {
            debug writeln("opSlice ", b, ' ', e);
            return Slice1D(b, e);
        }

        auto opDollar()
        {
            return length;
        }

        auto length() @property
        out
        {
            assert(__result != 0 || empty);
        }
        body 
        {
            debug writeln("length");
            auto tmp = backOffset;
            scope(exit) assert(backOffset == tmp);
            return sliceNI * sliceNJ - frontOffset - backOffset;
        }

        auto save() @property
        {
            return this;
        }

        static assert(isInputRange!(typeof(this)));
        static assert(isForwardRange!(typeof(this)));
        static assert(isBidirectionalRange!(typeof(this)));
        static assert(isRandomAccessRange!(typeof(this)));
    }

    auto byI()
    {
        return iota(nI).map!(i => this[i,0..$]);
    }

    auto byJ()
    {
        return iota(nJ).map!(j => this[0..$,j]);
    }

    size_t[2] shape() @property
    {
        return [slicing[0].e - slicing[0].b, slicing[1].e - slicing[1].b];
    }

    size_t size() @property
    {
        return (slicing[0].e - slicing[0].b) * (slicing[1].e - slicing[1].b);
    }
}

struct Slice1D
{
    size_t b, e;
}

private auto array2D(T)(T data, Slice1D[2] s, size_t nI, size_t nJ)
{
    return Array2D!T(data, s, nI, nJ);
}

auto array2D(T)(T data, size_t nI, size_t nJ)
{
    return Array2D!T(data, nI, nJ);
}

auto array2D(T)(size_t nI, size_t nJ)
{
    alias E = ElementType!T;
    return Array2D!T(new E[nI*nJ], nI, nJ);
}

unittest
{
    auto a2d = array2D!(long[])(4, 3);
    assert(a2d.data.length == 4*3);
    assert(a2d.size == 4*3);
    assert(a2d.shape == [4, 3]);

    auto byEl = a2d.byElement;
    assert(byEl.length == 4*3);

    a2d[] = 0;
    assert(a2d.byElement.equal(0.repeat(a2d.size)));

    iota(a2d.size).copy(a2d.data);

    assert(a2d.byElement.equal(iota(a2d.size)));
    assert(a2d.byElement[].equal(iota(a2d.size)));
    assert(a2d.byElement[0..$].equal(iota(a2d.size)));
    
    iota(a2d.size).map!"a*2".copy(a2d.byElement);
    assert(a2d.byElement[4..7].equal([8,10,12]), a2d.byElement[4..7].text);
    assert(a2d[1,0..$] == [6,8,10], a2d[1,0..$].text);

    a2d[1,0..$] = -1;
    assert(a2d.byElement.equal([0,2,4,-1,-1,-1,12,14,16,18,20,22]));

    a2d[] = 0;
    [1,2,1].copy(a2d[1..4, 1]);

    auto tmp0 = a2d*a2d;
    assert(tmp0.byElement.equal([
            0,0,0,
            0,1,0,
            0,4,0,
            0,1,0]));

    a2d[] = a2d*a2d;
    assert(a2d.byElement.equal([0,0,0,
                                0,1,0,
                                0,4,0,
                                0,1,0]));

    assert(a2d[3..$,2].equal([0]));

    auto b = a2d*-a2d;

    assert(b.byElement.equal([
            0, 0, 0,
            0,-1, 0,
            0,-16,0,
            0,-1, 0]), tmp0.byElement.text);
}

version(unittest)
{
    void main(){}
}
else
{
auto foo(Array2D!(float[]) a2d, size_t i0, size_t i1, size_t i2, float v)
{
    version(LDC) pragma(LDC_never_inline);
//    a2d = a2d[0..$/2, 0..$-1];
    a2d[i0 .. i1, i2] = v;

    import std.datetime : StopWatch;
    StopWatch sw;
    sw.start();
    //a2d[] = a2d*a2d;
    a2d[] = a2d - a2d;
    sw.peek().hnsecs.writeln;

    return a2d;
}

void main()
{
    auto a2d = foo(array2D!(float[])(1000, 1000), 2, 5, 3, 5.3f);
    //foreach(row; a2d.byElement.chunks(9))
    //    writeln(row);
}
}

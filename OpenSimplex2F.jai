
#import "Math";

#scope_file

PSIZE :: 2048;
PMASK :: 2047;
N2    :: 0.01001634121365712;
N3    :: 0.030485933181293584;
N4    :: 0.009202377986303158;

#scope_export

LatticePoint2D :: struct {
    xsv: int;
    ysv: int;
    dx: float64;
    dy: float64;
}

lattice_point_2d :: (xsv: int, ysv: int) -> LatticePoint2D {
    res : LatticePoint2D = ---;
    ssv : float64 = (xsv + ysv) * -0.211324865405187;
    res.xsv = xsv;
    res.ysv = ysv;
    res.dx  = -cast(float64)xsv - ssv; 
    res.dy  = -cast(float64)ysv - ssv;
    return res;
}

#scope_file
POINT_2D_ARRAY :: LatticePoint2D.[
    #run lattice_point_2d(1, 0),
    #run lattice_point_2d(0, 0),
    #run lattice_point_2d(1, 1),
    #run lattice_point_2d(0, 1),
];
#scope_export

LatticePoint3D :: struct {
    dxr: float64;
    dyr: float64;
    dzr: float64;
    xrv: int;
    yrv: int;
    zrv: int;
    nextOnFailure: int;
    nextOnSuccess: int;
}

lattice_point_3d :: (xrv: int, yrv: int, zrv: int, lattice: int) -> LatticePoint3D {
    res : LatticePoint3D = ---;
    res.dxr = -cast(float64)xrv + lattice * 0.5;
    res.dyr = -cast(float64)yrv + lattice * 0.5;
    res.dzr = -cast(float64)zrv + lattice * 0.5;
    res.xrv = xrv + lattice * 1024;
    res.yrv = yrv + lattice * 1024;
    res.zrv = zrv + lattice * 1024;
    res.nextOnFailure = -1;
    res.nextOnSuccess = -1;
    return res;
}

#scope_file

POINT_3D_ARRAY :: #run () -> [8*8]LatticePoint3D {
    res : [8*8]LatticePoint3D = ---;

    for i: 0..8-1 {
        i1, j1, k1, i2, j2, k2: int;
        i1 = (i >> 0) & 1;
        j1 = (i >> 1) & 1;
        k1 = (i >> 2) & 1;
        i2 = i1 ^ 1;
        j2 = j1 ^ 1;
        k2 = k1 ^ 1;

        // The two points within this octant, one from each of the two cubic half-lattices.
        q0 := i + 0*8;
        q1 := i + 1*8;
        q2 := i + 2*8;
        q3 := i + 3*8;
        q4 := i + 4*8;
        q5 := i + 5*8;
        q6 := i + 6*8;
        q7 := i + 7*8;

        res[q0] = lattice_point_3d(i1, j1, k1, 0);
        res[q1] = lattice_point_3d(i1+i2, j1+j2, k1+k2, 1);

        // Each single step away on the first half-lattice.
        res[q2] = lattice_point_3d(i1^1, j1, k1, 0);
        res[q3] = lattice_point_3d(i1, j1^1, k1, 0);
        res[q4] = lattice_point_3d(i1, j1, k1^1, 0);

        // Each single step away on the second half-lattice.
        res[q5] = lattice_point_3d(i1+(i2^1), j1+j2, k1+k2, 1);
        res[q6] = lattice_point_3d(i1+i2, j1+(j2^1), k1+k2, 1);
        res[q7] = lattice_point_3d(i1+i2, j1+j2, k1+(k2^1), 1);

        // First two are guaranteed.
        res[q0].nextOnFailure = q1;
        res[q0].nextOnSuccess = q1;
        res[q1].nextOnFailure = q2;
        res[q1].nextOnSuccess = q2;

        // Once we find one of the first half-lattice, the rest are out.
        // In addition, knowing c2 rules out c5.
        res[q2].nextOnFailure = q3;
        res[q2].nextOnSuccess = q6;
        res[q3].nextOnFailure = q4;
        res[q3].nextOnSuccess = q5;
        res[q4].nextOnFailure = q4;
        res[q4].nextOnSuccess = q5;

        // Once we find one on the second half-lattice, the rest are out.
        res[q5].nextOnFailure = q6;
        res[q5].nextOnSuccess = -1;
        res[q6].nextOnFailure = q7;
        res[q6].nextOnSuccess = -1;
        res[q7].nextOnFailure = -1;
        res[q7].nextOnSuccess = -1;
    }

    return res;    
}();

#scope_export

LatticePoint4D :: struct {
    xsv, ysv, zsv, wsv: int;
    dx, dy, dz, dw: float64;
    xsi, ysi, zsi, wsi: float64;
    ssiDelta: float64;
}

lattice_point_4d :: (xsv: int, ysv: int, zsv: int, wsv: int) -> LatticePoint4D {
    res : LatticePoint4D = ---;
    res.xsv = xsv + 409;
    res.ysv = ysv + 409;
    res.zsv = zsv + 409;
    res.wsv = wsv + 409;
    ssv := (xsv + ysv + zsv + wsv) * cast(float64)0.309016994374947;
    res.dx = -cast(float64)xsv - ssv;
    res.dy = -cast(float64)ysv - ssv;
    res.dz = -cast(float64)zsv - ssv;
    res.dw = -cast(float64)wsv - ssv;
    res.xsi = 0.2 - cast(float64)xsv;
    res.ysi = 0.2 - cast(float64)ysv;
    res.zsi = 0.2 - cast(float64)zsv;
    res.wsi = 0.2 - cast(float64)wsv;
    res.ssiDelta = (0.8 - cast(float64)(xsv + ysv + zsv + wsv)) * 0.309016994374947;
    return res;
}

#scope_file
POINT_4D_ARRAY :: #run () -> [16]LatticePoint4D {
    res : [16]LatticePoint4D = ---;
    for i: 0..16-1 res[i] = lattice_point_4d( (i>>0)&1, (i>>1)&1, (i>>2)&1, (i>>3)&1 );
    return res;    
}();
#scope_export

Grad2 :: struct {
    dx, dy: float64;
}
Grad3 :: struct {
    dx, dy, dz: float64;
}
Grad4 :: struct {
    dx, dy, dz, dw: float64;
}


fast_floor :: (x: float64) -> int {
    xi := cast(int)x;
    return ifx x < cast(float64)xi then xi-1 else xi;
}



OpenSimplex :: struct {
    perm: [PSIZE]int;
    permGrad2: [PSIZE]Grad2;
    permGrad3: [PSIZE]Grad3;
    permGrad4: [PSIZE]Grad4;
}

open_simplex :: (seed: int) -> OpenSimplex {
    res : OpenSimplex = ---;
    source : [PSIZE]int = ---;
    for *source it[0] = it_index;

    for i: 0..PSIZE-1 {
        seed = seed * 6364136223846793005 + 1442695040888963407;
        r := (seed + 31) % (i+1);
        if r < 0 then r += (i+1);

        res.perm[i] = source[r];
        res.permGrad2[i] = GRADIENTS_2D[res.perm[i]];
        res.permGrad3[i] = GRADIENTS_3D[res.perm[i]];
        res.permGrad4[i] = GRADIENTS_4D[res.perm[i]];
        source[r] = source[i];
    }

    return res;
}



noise2_base :: (os: OpenSimplex, xs: float64, ys: float64) -> float64 {
    value : float64 = 0.0;

    // Get base points and offsets
    xsb := fast_floor(xs);
    ysb := fast_floor(ys);
    xsi := xs - xsb;
    ysi := ys - ysb;

    // Index to point list
    index := cast(int)((ysi-xsi)/2 + 1);

    ssi := (xsi+ysi) * -0.211324865405187;
    xi := xsi + ssi;
    yi := ysi + ssi;

    // Point contributions
    for i: 0..3-1 {
        c := POINT_2D_ARRAY[index + i];

        dx := xi + c.dx;
        dy := yi + c.dy;
        attn := 0.5 - dx*dx - dy*dy;
        if attn <= 0 continue;

        pxm := (xsb + c.xsv) & PMASK;
        pym := (ysb + c.ysv) & PMASK;
        grad := os.permGrad2[os.perm[pxm] ^ pym];
        extrapolation := grad.dx * dx + grad.dy * dy;

        attn *= attn;
        value += attn * attn * extrapolation;
    }

    return value;
}


/* 
    2D Simplex noise, standard lattice orientation 
*/
noise2 :: (os: OpenSimplex, x: float64, y: float64) -> float64 {
    // Get points for A2* lattice
    s := 0.366025403784439 * (x + y);
    xs := x + s;
    ys := y + s;

    return noise2_base(os, xs, ys);
}

/*
    2D Simplex noise, with Y pointing down the main diagonal.
    Might be better for a 2D sandbox style game, where Y is vertical.
    Probably slightly less optimal for heightmaps or continent maps.
*/
noise2_x_before_y :: (os: OpenSimplex, x: float64, y: float64) -> float64 {
    // Skew transform and rotation baked into one.
    x2 := x * 0.7071067811865476;
    y2 := y * 1.224744871380249;
    return noise2_base(os, y2 + x2, y2 - x2);
}

/*
    Generate overlapping cubic lattices for 3D re-oriented BCC noise.
    Lookup table implementation inspired by DigitalShadow.
    It was actually faster to narrow down the points in the loop itself,
    than to build up the index with enough info to isolate 4 points.
*/
noise3_bcc :: (os: OpenSimplex, xr: float64, yr: float64, zr: float64) -> float64 {
    // Get base and offsets inside cube of first lattice.
    xrb := fast_floor(xr);
    yrb := fast_floor(yr);
    zrb := fast_floor(zr);
    xri := xr - xrb;
    yri := yr - yrb;
    zri := zr - zrb;

    // Identify which octant of the cube we're in. This determines which cell
    // in the other cubic lattice we're in, and also narrows down one point on each.
    xht := cast(int)(xri + 0.5);
    yht := cast(int)(yri + 0.5);
    zht := cast(int)(zri + 0.5);
    index := (xht << 0) | (yht << 1) | (zht << 2);

    // Point contributions
    value : float64 = 0.0;
    c := POINT_3D_ARRAY[index];
    while true {
        dxr := xri + c.dxr;
        dyr := yri + c.dyr;
        dzr := zri + c.dzr;
        attn := 0.5 - dxr * dxr - dyr * dyr - dzr * dzr;

        if attn < 0 then index = c.nextOnFailure;
        else {
            pxm := (xrb + c.xrv) & PMASK;
            pym := (yrb + c.yrv) & PMASK;
            pzm := (zrb + c.zrv) & PMASK;

            grad := os.permGrad3[os.perm[os.perm[pxm] ^ pym] ^ pzm];
            extrapolation := grad.dx * dxr + grad.dy * dyr + grad.dz * dzr;

            attn *= attn;
            value += attn * attn * extrapolation;
            index = c.nextOnSuccess;
        }

        if index == -1 break;
        c = POINT_3D_ARRAY[index];
    }

    return value;
}


/*
    3D re-oriented 4-point BBC noise, classic orientation.
    Proper substititue for 3D Simplex in light of Forbidden Formulae.
    Use noise3_xy_before_z or noise3_xz_before_y instead, wherever appropriate.
*/
noise3_classic :: (os: OpenSimplex, x: float64, y: float64, z: float64) -> float64 {
    // Re-orient the cubic lattices via rotation, to produce the expected look on cardinal planar slices.
    // If texturing objects that don't tend to have cardinal plane faces, you could even remove this.
    // Orthonormal rotation. Not a skew transform.
    r := (2.0/3.0) * (x+y+z);
    xr := r-x;
    yr := r-y;
    zr := r-z;
    return noise3_bcc(os, xr, yr, zr);
}

/*
    3D re-oriented 4-point BBC noise, with better visual isotropy in (X,Y).
    Recommended for 3D terrain and time-varied animations.
    The Z coordinate should always be the "different" coordinate in your use case.
    If Y is vertical in world coordinates, call noise3_xy_before_z(x,z,Y) or use noise3_xz_before_y.
    If Z is vertical in world coordinates, call noise3_xy_before_z(x,y,Z).
    For a time-varied animation, call noise3_xy_before_z(x,y,T).
*/
noise3_xy_before_z :: (os: OpenSimplex, x: float64, y: float64, z: float64) -> float64 {
    // Re-orient the cubic lattices without skewing, to make X and Y triangular like 2D.
    // Orthonormal rotation. Not a skew transform.
    xy := x + y;
    s2 := xy * -0.211324865405187;
    zz := z * 0.577350269189626;
    xr := x + s2 - zz;
    yr := y + s2 - zz;
    zr := xy * 0.577350269189626 + zz;

    // Evaluate both lattices to form a BCC lattice.
    return noise3_bcc(os, xr, yr, zr);
}


/*
    3D re-oriented 4-point BCC noise, with better visual isotropy in (X,Z).
    Recommended for 3D terrain and time-varied animations.
    The Y coordinate should always be the "different" coordinate in your use case.
    If Y is vertical in world coordinates, call noise3_xz_before_y(x,Y,z).
    If Z is vertical in world coordinates, call noise3_xz_before_Y(x,Z,y) or use noise3_xy_before_z.
    For a time-varied animation, call noise3_xz_before_y(x,T,y) or use noise3_xy_before_z.
*/
noise3_xz_before_y :: (os: OpenSimplex, x: float64, y: float64, z: float64) -> float64 {
    // Re-orient the cubic lattices without skewing, to make X and Z triangular like 2D.
    // Orthonormal rotations. Not a skew transform.
    xz := x + z;
    s2 := xz * -0.211324865405187;
    yy := y * 0.577350269189626;
    xr := x + s2 - yy;
    zr := z + s2 - yy;
    yr := xz * 0.577350269189626 + yy;

    // Evaluate both lattices to form a BCC lattice.
    return noise3_bcc(os, xr, yr, zr);
}


/*
    4D OpenSimplex2F noise base.
    Current implementations not fully optimized by lookup tables.
    But still comes out slightly ahead of Gustavson's Simplex in tests.
*/
noise4_base :: (os: OpenSimplex, xs: float64, ys: float64, zs: float64, ws: float64) -> float64 {
    value : float64 = 0.0;

    // Get base points and offsets
    xsb := fast_floor(xs);
    ysb := fast_floor(ys);
    zsb := fast_floor(zs);
    wsb := fast_floor(ws);
    xsi := xs - xsb;
    ysi := ys - ysb;
    zsi := zs - zsb;
    wsi := ws - wsb;

    // If we're in the lower half, flip, so we can repeat the code for the upper half.
    // We'll flip back later.
    si_sum := xsi + ysi + zsi + wsi;
    ssi := si_sum * 0.309016994374947; // Prep for vertex contributions.
    in_lower_half := si_sum < 2;
    if in_lower_half {
        xsi = 1 - xsi;
        ysi = 1 - ysi;
        zsi = 1 - zsi;
        wsi = 1 - wsi;
        si_sum = 4 - si_sum;
    }

    // Consider opposing vertex pairs of the octahedron formed by the central cross-section of the stretched tesseract
    aabb := xsi + ysi - zsi - wsi;
    abab := xsi - ysi + zsi - wsi;
    abba := xsi - ysi - zsi + wsi;
    aabb_score := Math.abs(aabb);
    abab_score := Math.abs(abab);
    abba_score := Math.abs(abba);

    // Find the closest point on the stretched tesseract as if it were the upper half
    vertex_index : int = ---;
    via : int = ---;
    vib : int = ---;
    asi: float64 = ---;
    bsi : float64 = ---;

    if aabb_score > abab_score && aabb_score > abba_score {
        if aabb > 0 {
            asi = zsi;
            bsi = wsi;
            vertex_index = 0b0011;
            via = 0b0111;
            vib = 0b1011;
        } else {
            asi = xsi;
            bsi = ysi;
            vertex_index = 0b1100;
            via = 0b1101;
            vib = 0b1110;
        }
    } else if abab_score > abba_score {
        if abab > 0 {
            asi = ysi;
            bsi = wsi;
            vertex_index = 0b0101;
            via = 0b0111;
            vib = 0b1101;
        } else {
            asi = xsi;
            bsi = zsi;
            vertex_index = 0b1010;
            via = 0b1011;
            vib = 0b1110;
        }
    } else {
        if abba > 0 {
            asi = ysi;
            bsi = zsi;
            vertex_index = 0b1001;
            via = 0b1011;
            vib = 0b1101;
        } else {
            asi = xsi;
            bsi = wsi;
            vertex_index = 0b0110;
            via = 0b0111;
            vib = 0b1110;
        }
    }

    if bsi > asi {
        via = vib;
        tmp := bsi;
        bsi = asi;
        asi = tmp;
    }

    if si_sum + asi > 3 {
        vertex_index = via;
        if si_sum + bsi > 4 then vertex_index = 0b1111;
    }

    // Now flip back if we're actually in the lower half.
    if in_lower_half {
        xsi = 1 - xsi;
        ysi = 1 - ysi;
        zsi = 1 - zsi;
        wsi = 1 - wsi;
        vertex_index ^= 0b1111;
    }

    // Five points to add, total, from five copies of the A4 lattice.
    for i: 0..5-1 {
        c := POINT_4D_ARRAY[vertex_index];
        xsb += c.xsv;
        ysb += c.ysv;
        zsb += c.zsv;
        wsb += c.wsv;

        xi := xsi + ssi;
        yi := ysi + ssi;
        zi := zsi + ssi;
        wi := wsi + ssi;

        dx := xi + c.dx;
        dy := yi + c.dy;
        dz := zi + c.dz;
        dw := wi + c.dw;

        attn := 0.5 - dx*dx - dy*dy - dz*dz - dw*dw;
        if attn > 0 {
            pxm := xsb & PMASK;
            pym := ysb & PMASK;
            pzm := zsb & PMASK;
            pwm := wsb & PMASK;

            grad := os.permGrad4[os.perm[os.perm[os.perm[pxm]^pym]^pzm]^pwm];
            ramped := grad.dx * dx + grad.dy * dy + grad.dz * dz + grad.dw * dw;

            attn *= attn;
            value += attn * attn * ramped;
        }

        // Maybe this helps the compiler/JVM/LLVM/etc. know we can end the loop here. Maybe not.
        if i == 4 break;

        // Update the relative skewed coordinates to reference the vertex we just added.
        // Rather, reference its counterpart on the lattice copy that is shifted down by
        // the vector (-0.2, -0.2, -0.2, -0.2)
        xsi += c.xsi;
        ysi += c.ysi;
        zsi += c.zsi;
        wsi += c.wsi;
        ssi += c.ssiDelta;

        // Next point is the closest vertex on the 4-simplex whose base vertex is the aforementioned vertex.
        score0 := 1.0 + ssi * (-1.0/0.309016994374947); // Seems slightly faster than 1.0-xsi-ysi-zsi-wsi
        vertex_index = 0b0000;
        if xsi >= ysi && xsi >= zsi && xsi >= wsi && xsi >= score0 then vertex_index = 0b0001;
        if ysi >= xsi && ysi >= zsi && ysi >= wsi && ysi >= score0 then vertex_index = 0b0010;
        if zsi >= xsi && zsi >= ysi && zsi >= wsi && zsi >= score0 then vertex_index = 0b0100;
        if wsi >= xsi && wsi >= ysi && wsi >= zsi && wsi >= score0 then vertex_index = 0b1000;
    }

    return value;
}


/*
    4D OpenSimplex2F noise, classic lattice orientation.
*/
noise4_classic :: (os: OpenSimplex, x: float64, y: float64, z: float64, w: float64) -> float64 {
    // Get points for A4 lattice
    s := -0.138196601125011 * (x + y + z + w);
    xs := x + s;
    ys := y + s;
    zs := z + s;
    ws := w + s;

    return noise4_base(os, xs, ys, zs, ws);
}

/*
    4D OpenSimplex3F noise, with XY and ZW forming orthogonal triangular-based planes.
    Recommended for 3D terrain, where X and Y (or Z and W) are horizontal.
    Recommended for noise(x, y, sin(time), cos(time)) trick.
*/
noise4_xy_before_zw :: (os: OpenSimplex, x: float64, y: float64, z: float64, w: float64) -> float64 {
    s2 := (x + y) * -0.178275657951399372 + (z + w) * 0.215623393288842828;
    t2 := (z + w) * -0.403949762580207112 + (x + y) * -0.375199083010075342;
    xs := x + s2;
    ys := y + s2;
    zs := z + t2;
    ws := w + t2;

    return noise4_base(os, xs, ys, zs, ws);
}

/*
    4D OpenSimplex2F noise, with XZ and YW forming orthogonal triangular-based planes.
    Recommended for 3D terrain, where X and Z (or Y and W) are horizontal.
*/
noise4_xz_before_yw :: (os: OpenSimplex, x: float64, y: float64, z: float64, w: float64) -> float64 {
    s2 := (x + z) * -0.178275657951399372 + (y + w) * 0.215623393288842828;
    t2 := (y + w) * -0.403949762580207112 + (x + z) * -0.375199083010075342;
    xs := x + s2;
    ys := y + t2;
    zs := z + s2;
    ws := w + t2;

    return noise4_base(os, xs, ys, zs, ws);
}

/*
    4D OpenSimplex3F noise, with XYZ oriented like noise3_classic,
    and W for an extra degree of freedom. W repeats eventually.
    Recommended for time-varied animations which texture a 3D object (W = time)
*/
noise4_xyz_before_w :: (os: OpenSimplex, x: float64, y: float64, z: float64, w: float64) -> float64 {
    xyz := x + y + z;
    ww := w * 0.2236067977499788;
    s2 := xyz * -0.16666666666666666 + ww;
    xs := x + s2;
    ys := y + s2;
    zs := z + s2;
    ws := -0.5 * xyz + ww;

    return noise4_base(os, xs, ys, zs, ws);
}



#scope_file
GRADIENTS_2D :: #run ()->[PSIZE]Grad2 {
    res : [PSIZE]Grad2 = ---;

    arr :: Grad2.[
        Grad2.{0.130526192220052  /N2, 0.99144486137381   /N2},
        Grad2.{0.38268343236509   /N2, 0.923879532511287  /N2},
        Grad2.{0.608761429008721  /N2, 0.793353340291235  /N2},
        Grad2.{0.793353340291235  /N2, 0.608761429008721  /N2},
        Grad2.{0.923879532511287  /N2, 0.38268343236509   /N2},
        Grad2.{0.99144486137381   /N2, 0.130526192220051  /N2},
        Grad2.{0.99144486137381   /N2, -0.130526192220051 /N2},
        Grad2.{0.923879532511287  /N2, -0.38268343236509  /N2},
        Grad2.{0.793353340291235  /N2, -0.60876142900872  /N2},
        Grad2.{0.608761429008721  /N2, -0.793353340291235 /N2},
        Grad2.{0.38268343236509   /N2, -0.923879532511287 /N2},
        Grad2.{0.130526192220052  /N2, -0.99144486137381  /N2},
        Grad2.{-0.130526192220052 /N2, -0.99144486137381  /N2},
        Grad2.{-0.38268343236509  /N2, -0.923879532511287 /N2},
        Grad2.{-0.608761429008721 /N2, -0.793353340291235 /N2},
        Grad2.{-0.793353340291235 /N2, -0.608761429008721 /N2},
        Grad2.{-0.923879532511287 /N2, -0.38268343236509  /N2},
        Grad2.{-0.99144486137381  /N2, -0.130526192220052 /N2},
        Grad2.{-0.99144486137381  /N2, 0.130526192220051  /N2},
        Grad2.{-0.923879532511287 /N2, 0.38268343236509   /N2},
        Grad2.{-0.793353340291235 /N2, 0.608761429008721  /N2},
        Grad2.{-0.608761429008721 /N2, 0.793353340291235  /N2},
        Grad2.{-0.38268343236509  /N2, 0.923879532511287  /N2},
        Grad2.{-0.130526192220052 /N2, 0.99144486137381   /N2},
    ];

    for i: 0..PSIZE-1 res[i] = arr[i%24];

    return res;    
}();

GRADIENTS_3D :: #run ()->[PSIZE]Grad3 {
    res : [PSIZE]Grad3 = ---;

    arr :: Grad3.[
        Grad3.{-2.22474487139      /N3 , -2.22474487139      /N3 , -1.0                /N3 },
        Grad3.{-2.22474487139      /N3 , -2.22474487139      /N3 ,  1.0                /N3 },
        Grad3.{-3.0862664687972017 /N3 , -1.1721513422464978 /N3 ,  0.0                /N3 },
        Grad3.{-1.1721513422464978 /N3 , -3.0862664687972017 /N3 ,  0.0                /N3 },
        Grad3.{-2.22474487139      /N3 , -1.0                /N3 , -2.22474487139      /N3 },
        Grad3.{-2.22474487139      /N3 ,  1.0                /N3 , -2.22474487139      /N3 },
        Grad3.{-1.1721513422464978 /N3 ,  0.0                /N3 , -3.0862664687972017 /N3 },
        Grad3.{-3.0862664687972017 /N3 ,  0.0                /N3 , -1.1721513422464978 /N3 },
        Grad3.{-2.22474487139      /N3 , -1.0                /N3 ,  2.22474487139      /N3 },
        Grad3.{-2.22474487139      /N3 ,  1.0                /N3 ,  2.22474487139      /N3 },
        Grad3.{-3.0862664687972017 /N3 ,  0.0                /N3 ,  1.1721513422464978 /N3 },
        Grad3.{-1.1721513422464978 /N3 ,  0.0                /N3 ,  3.0862664687972017 /N3 },
        Grad3.{-2.22474487139      /N3 ,  2.22474487139      /N3 , -1.0                /N3 },
        Grad3.{-2.22474487139      /N3 ,  2.22474487139      /N3 ,  1.0                /N3 },
        Grad3.{-1.1721513422464978 /N3 ,  3.0862664687972017 /N3 ,  0.0                /N3 },
        Grad3.{-3.0862664687972017 /N3 ,  1.1721513422464978 /N3 ,  0.0                /N3 },
        Grad3.{-1.0                /N3 , -2.22474487139      /N3 , -2.22474487139      /N3 },
        Grad3.{ 1.0                /N3 , -2.22474487139      /N3 , -2.22474487139      /N3 },
        Grad3.{ 0.0                /N3 , -3.0862664687972017 /N3 , -1.1721513422464978 /N3 },
        Grad3.{ 0.0                /N3 , -1.1721513422464978 /N3 , -3.0862664687972017 /N3 },
        Grad3.{-1.0                /N3 , -2.22474487139      /N3 ,  2.22474487139      /N3 },
        Grad3.{ 1.0                /N3 , -2.22474487139      /N3 ,  2.22474487139      /N3 },
        Grad3.{ 0.0                /N3 , -1.1721513422464978 /N3 ,  3.0862664687972017 /N3 },
        Grad3.{ 0.0                /N3 , -3.0862664687972017 /N3 ,  1.1721513422464978 /N3 },
        Grad3.{-1.0                /N3 ,  2.22474487139      /N3 , -2.22474487139      /N3 },
        Grad3.{ 1.0                /N3 ,  2.22474487139      /N3 , -2.22474487139      /N3 },
        Grad3.{ 0.0                /N3 ,  1.1721513422464978 /N3 , -3.0862664687972017 /N3 },
        Grad3.{ 0.0                /N3 ,  3.0862664687972017 /N3 , -1.1721513422464978 /N3 },
        Grad3.{-1.0                /N3 ,  2.22474487139      /N3 ,  2.22474487139      /N3 },
        Grad3.{ 1.0                /N3 ,  2.22474487139      /N3 ,  2.22474487139      /N3 },
        Grad3.{ 0.0                /N3 ,  3.0862664687972017 /N3 ,  1.1721513422464978 /N3 },
        Grad3.{ 0.0                /N3 ,  1.1721513422464978 /N3 ,  3.0862664687972017 /N3 },
        Grad3.{ 2.22474487139      /N3 , -2.22474487139      /N3 , -1.0                /N3 },
        Grad3.{ 2.22474487139      /N3 , -2.22474487139      /N3 ,  1.0                /N3 },
        Grad3.{ 1.1721513422464978 /N3 , -3.0862664687972017 /N3 ,  0.0                /N3 },
        Grad3.{ 3.0862664687972017 /N3 , -1.1721513422464978 /N3 ,  0.0                /N3 },
        Grad3.{ 2.22474487139      /N3 , -1.0                /N3 , -2.22474487139      /N3 },
        Grad3.{ 2.22474487139      /N3 ,  1.0                /N3 , -2.22474487139      /N3 },
        Grad3.{ 3.0862664687972017 /N3 ,  0.0                /N3 , -1.1721513422464978 /N3 },
        Grad3.{ 1.1721513422464978 /N3 ,  0.0                /N3 , -3.0862664687972017 /N3 },
        Grad3.{ 2.22474487139      /N3 , -1.0                /N3 ,  2.22474487139      /N3 },
        Grad3.{ 2.22474487139      /N3 ,  1.0                /N3 ,  2.22474487139      /N3 },
        Grad3.{ 1.1721513422464978 /N3 ,  0.0                /N3 ,  3.0862664687972017 /N3 },
        Grad3.{ 3.0862664687972017 /N3 ,  0.0                /N3 ,  1.1721513422464978 /N3 },
        Grad3.{ 2.22474487139      /N3 ,  2.22474487139      /N3 , -1.0                /N3 },
        Grad3.{ 2.22474487139      /N3 ,  2.22474487139      /N3 ,  1.0                /N3 },
        Grad3.{ 3.0862664687972017 /N3 ,  1.1721513422464978 /N3 ,  0.0                /N3 },
        Grad3.{ 1.1721513422464978 /N3 ,  3.0862664687972017 /N3 ,  0.0                /N3 },
    ];

    for i: 0..PSIZE-1 res[i] = arr[i%48];

    return res;    
}();


GRADIENTS_4D :: #run ()->[PSIZE]Grad4 {
    res : [PSIZE]Grad4 = ---;

    arr :: Grad4.[
        Grad4.{-0.753341017856078 /N4, -0.37968289875261624 /N4, -0.37968289875261624 /N4, -0.37968289875261624},
        Grad4.{-0.7821684431180708 /N4, -0.4321472685365301 /N4, -0.4321472685365301 /N4, 0.12128480194602098},
        Grad4.{-0.7821684431180708 /N4, -0.4321472685365301 /N4, 0.12128480194602098 /N4, -0.4321472685365301},
        Grad4.{-0.7821684431180708 /N4, 0.12128480194602098 /N4, -0.4321472685365301 /N4, -0.4321472685365301},
        Grad4.{-0.8586508742123365 /N4, -0.508629699630796 /N4, 0.044802370851755174 /N4, 0.044802370851755174},
        Grad4.{-0.8586508742123365 /N4, 0.044802370851755174 /N4, -0.508629699630796 /N4, 0.044802370851755174},
        Grad4.{-0.8586508742123365 /N4, 0.044802370851755174 /N4, 0.044802370851755174 /N4, -0.508629699630796},
        Grad4.{-0.9982828964265062 /N4, -0.03381941603233842 /N4, -0.03381941603233842 /N4, -0.03381941603233842},
        Grad4.{-0.37968289875261624 /N4, -0.753341017856078 /N4, -0.37968289875261624 /N4, -0.37968289875261624},
        Grad4.{-0.4321472685365301 /N4, -0.7821684431180708 /N4, -0.4321472685365301 /N4, 0.12128480194602098},
        Grad4.{-0.4321472685365301 /N4, -0.7821684431180708 /N4, 0.12128480194602098 /N4, -0.4321472685365301},
        Grad4.{0.12128480194602098 /N4, -0.7821684431180708 /N4, -0.4321472685365301 /N4, -0.4321472685365301},
        Grad4.{-0.508629699630796 /N4, -0.8586508742123365 /N4, 0.044802370851755174 /N4, 0.044802370851755174},
        Grad4.{0.044802370851755174 /N4, -0.8586508742123365 /N4, -0.508629699630796 /N4, 0.044802370851755174},
        Grad4.{0.044802370851755174 /N4, -0.8586508742123365 /N4, 0.044802370851755174 /N4, -0.508629699630796},
        Grad4.{-0.03381941603233842 /N4, -0.9982828964265062 /N4, -0.03381941603233842 /N4, -0.03381941603233842},
        Grad4.{-0.37968289875261624 /N4, -0.37968289875261624 /N4, -0.753341017856078 /N4, -0.37968289875261624},
        Grad4.{-0.4321472685365301 /N4, -0.4321472685365301 /N4, -0.7821684431180708 /N4, 0.12128480194602098},
        Grad4.{-0.4321472685365301 /N4, 0.12128480194602098 /N4, -0.7821684431180708 /N4, -0.4321472685365301},
        Grad4.{0.12128480194602098 /N4, -0.4321472685365301 /N4, -0.7821684431180708 /N4, -0.4321472685365301},
        Grad4.{-0.508629699630796 /N4, 0.044802370851755174 /N4, -0.8586508742123365 /N4, 0.044802370851755174},
        Grad4.{0.044802370851755174 /N4, -0.508629699630796 /N4, -0.8586508742123365 /N4, 0.044802370851755174},
        Grad4.{0.044802370851755174 /N4, 0.044802370851755174 /N4, -0.8586508742123365 /N4, -0.508629699630796},
        Grad4.{-0.03381941603233842 /N4, -0.03381941603233842 /N4, -0.9982828964265062 /N4, -0.03381941603233842},
        Grad4.{-0.37968289875261624 /N4, -0.37968289875261624 /N4, -0.37968289875261624 /N4, -0.753341017856078},
        Grad4.{-0.4321472685365301 /N4, -0.4321472685365301 /N4, 0.12128480194602098 /N4, -0.7821684431180708},
        Grad4.{-0.4321472685365301 /N4, 0.12128480194602098 /N4, -0.4321472685365301 /N4, -0.7821684431180708},
        Grad4.{0.12128480194602098 /N4, -0.4321472685365301 /N4, -0.4321472685365301 /N4, -0.7821684431180708},
        Grad4.{-0.508629699630796 /N4, 0.044802370851755174 /N4, 0.044802370851755174 /N4, -0.8586508742123365},
        Grad4.{0.044802370851755174 /N4, -0.508629699630796 /N4, 0.044802370851755174 /N4, -0.8586508742123365},
        Grad4.{0.044802370851755174 /N4, 0.044802370851755174 /N4, -0.508629699630796 /N4, -0.8586508742123365},
        Grad4.{-0.03381941603233842 /N4, -0.03381941603233842 /N4, -0.03381941603233842 /N4, -0.9982828964265062},
        Grad4.{-0.6740059517812944 /N4, -0.3239847771997537 /N4, -0.3239847771997537 /N4, 0.5794684678643381},
        Grad4.{-0.7504883828755602 /N4, -0.4004672082940195 /N4, 0.15296486218853164 /N4, 0.5029860367700724},
        Grad4.{-0.7504883828755602 /N4, 0.15296486218853164 /N4, -0.4004672082940195 /N4, 0.5029860367700724},
        Grad4.{-0.8828161875373585 /N4, 0.08164729285680945 /N4, 0.08164729285680945 /N4, 0.4553054119602712},
        Grad4.{-0.4553054119602712 /N4, -0.08164729285680945 /N4, -0.08164729285680945 /N4, 0.8828161875373585},
        Grad4.{-0.5029860367700724 /N4, -0.15296486218853164 /N4, 0.4004672082940195 /N4, 0.7504883828755602},
        Grad4.{-0.5029860367700724 /N4, 0.4004672082940195 /N4, -0.15296486218853164 /N4, 0.7504883828755602},
        Grad4.{-0.5794684678643381 /N4, 0.3239847771997537 /N4, 0.3239847771997537 /N4, 0.6740059517812944},
        Grad4.{-0.3239847771997537 /N4, -0.6740059517812944 /N4, -0.3239847771997537 /N4, 0.5794684678643381},
        Grad4.{-0.4004672082940195 /N4, -0.7504883828755602 /N4, 0.15296486218853164 /N4, 0.5029860367700724},
        Grad4.{0.15296486218853164 /N4, -0.7504883828755602 /N4, -0.4004672082940195 /N4, 0.5029860367700724},
        Grad4.{0.08164729285680945 /N4, -0.8828161875373585 /N4, 0.08164729285680945 /N4, 0.4553054119602712},
        Grad4.{-0.08164729285680945 /N4, -0.4553054119602712 /N4, -0.08164729285680945 /N4, 0.8828161875373585},
        Grad4.{-0.15296486218853164 /N4, -0.5029860367700724 /N4, 0.4004672082940195 /N4, 0.7504883828755602},
        Grad4.{0.4004672082940195 /N4, -0.5029860367700724 /N4, -0.15296486218853164 /N4, 0.7504883828755602},
        Grad4.{0.3239847771997537 /N4, -0.5794684678643381 /N4, 0.3239847771997537 /N4, 0.6740059517812944},
        Grad4.{-0.3239847771997537 /N4, -0.3239847771997537 /N4, -0.6740059517812944 /N4, 0.5794684678643381},
        Grad4.{-0.4004672082940195 /N4, 0.15296486218853164 /N4, -0.7504883828755602 /N4, 0.5029860367700724},
        Grad4.{0.15296486218853164 /N4, -0.4004672082940195 /N4, -0.7504883828755602 /N4, 0.5029860367700724},
        Grad4.{0.08164729285680945 /N4, 0.08164729285680945 /N4, -0.8828161875373585 /N4, 0.4553054119602712},
        Grad4.{-0.08164729285680945 /N4, -0.08164729285680945 /N4, -0.4553054119602712 /N4, 0.8828161875373585},
        Grad4.{-0.15296486218853164 /N4, 0.4004672082940195 /N4, -0.5029860367700724 /N4, 0.7504883828755602},
        Grad4.{0.4004672082940195 /N4, -0.15296486218853164 /N4, -0.5029860367700724 /N4, 0.7504883828755602},
        Grad4.{0.3239847771997537 /N4, 0.3239847771997537 /N4, -0.5794684678643381 /N4, 0.6740059517812944},
        Grad4.{-0.6740059517812944 /N4, -0.3239847771997537 /N4, 0.5794684678643381 /N4, -0.3239847771997537},
        Grad4.{-0.7504883828755602 /N4, -0.4004672082940195 /N4, 0.5029860367700724 /N4, 0.15296486218853164},
        Grad4.{-0.7504883828755602 /N4, 0.15296486218853164 /N4, 0.5029860367700724 /N4, -0.4004672082940195},
        Grad4.{-0.8828161875373585 /N4, 0.08164729285680945 /N4, 0.4553054119602712 /N4, 0.08164729285680945},
        Grad4.{-0.4553054119602712 /N4, -0.08164729285680945 /N4, 0.8828161875373585 /N4, -0.08164729285680945},
        Grad4.{-0.5029860367700724 /N4, -0.15296486218853164 /N4, 0.7504883828755602 /N4, 0.4004672082940195},
        Grad4.{-0.5029860367700724 /N4, 0.4004672082940195 /N4, 0.7504883828755602 /N4, -0.15296486218853164},
        Grad4.{-0.5794684678643381 /N4, 0.3239847771997537 /N4, 0.6740059517812944 /N4, 0.3239847771997537},
        Grad4.{-0.3239847771997537 /N4, -0.6740059517812944 /N4, 0.5794684678643381 /N4, -0.3239847771997537},
        Grad4.{-0.4004672082940195 /N4, -0.7504883828755602 /N4, 0.5029860367700724 /N4, 0.15296486218853164},
        Grad4.{0.15296486218853164 /N4, -0.7504883828755602 /N4, 0.5029860367700724 /N4, -0.4004672082940195},
        Grad4.{0.08164729285680945 /N4, -0.8828161875373585 /N4, 0.4553054119602712 /N4, 0.08164729285680945},
        Grad4.{-0.08164729285680945 /N4, -0.4553054119602712 /N4, 0.8828161875373585 /N4, -0.08164729285680945},
        Grad4.{-0.15296486218853164 /N4, -0.5029860367700724 /N4, 0.7504883828755602 /N4, 0.4004672082940195},
        Grad4.{0.4004672082940195 /N4, -0.5029860367700724 /N4, 0.7504883828755602 /N4, -0.15296486218853164},
        Grad4.{0.3239847771997537 /N4, -0.5794684678643381 /N4, 0.6740059517812944 /N4, 0.3239847771997537},
        Grad4.{-0.3239847771997537 /N4, -0.3239847771997537 /N4, 0.5794684678643381 /N4, -0.6740059517812944},
        Grad4.{-0.4004672082940195 /N4, 0.15296486218853164 /N4, 0.5029860367700724 /N4, -0.7504883828755602},
        Grad4.{0.15296486218853164 /N4, -0.4004672082940195 /N4, 0.5029860367700724 /N4, -0.7504883828755602},
        Grad4.{0.08164729285680945 /N4, 0.08164729285680945 /N4, 0.4553054119602712 /N4, -0.8828161875373585},
        Grad4.{-0.08164729285680945 /N4, -0.08164729285680945 /N4, 0.8828161875373585 /N4, -0.4553054119602712},
        Grad4.{-0.15296486218853164 /N4, 0.4004672082940195 /N4, 0.7504883828755602 /N4, -0.5029860367700724},
        Grad4.{0.4004672082940195 /N4, -0.15296486218853164 /N4, 0.7504883828755602 /N4, -0.5029860367700724},
        Grad4.{0.3239847771997537 /N4, 0.3239847771997537 /N4, 0.6740059517812944 /N4, -0.5794684678643381},
        Grad4.{-0.6740059517812944 /N4, 0.5794684678643381 /N4, -0.3239847771997537 /N4, -0.3239847771997537},
        Grad4.{-0.7504883828755602 /N4, 0.5029860367700724 /N4, -0.4004672082940195 /N4, 0.15296486218853164},
        Grad4.{-0.7504883828755602 /N4, 0.5029860367700724 /N4, 0.15296486218853164 /N4, -0.4004672082940195},
        Grad4.{-0.8828161875373585 /N4, 0.4553054119602712 /N4, 0.08164729285680945 /N4, 0.08164729285680945},
        Grad4.{-0.4553054119602712 /N4, 0.8828161875373585 /N4, -0.08164729285680945 /N4, -0.08164729285680945},
        Grad4.{-0.5029860367700724 /N4, 0.7504883828755602 /N4, -0.15296486218853164 /N4, 0.4004672082940195},
        Grad4.{-0.5029860367700724 /N4, 0.7504883828755602 /N4, 0.4004672082940195 /N4, -0.15296486218853164},
        Grad4.{-0.5794684678643381 /N4, 0.6740059517812944 /N4, 0.3239847771997537 /N4, 0.3239847771997537},
        Grad4.{-0.3239847771997537 /N4, 0.5794684678643381 /N4, -0.6740059517812944 /N4, -0.3239847771997537},
        Grad4.{-0.4004672082940195 /N4, 0.5029860367700724 /N4, -0.7504883828755602 /N4, 0.15296486218853164},
        Grad4.{0.15296486218853164 /N4, 0.5029860367700724 /N4, -0.7504883828755602 /N4, -0.4004672082940195},
        Grad4.{0.08164729285680945 /N4, 0.4553054119602712 /N4, -0.8828161875373585 /N4, 0.08164729285680945},
        Grad4.{-0.08164729285680945 /N4, 0.8828161875373585 /N4, -0.4553054119602712 /N4, -0.08164729285680945},
        Grad4.{-0.15296486218853164 /N4, 0.7504883828755602 /N4, -0.5029860367700724 /N4, 0.4004672082940195},
        Grad4.{0.4004672082940195 /N4, 0.7504883828755602 /N4, -0.5029860367700724 /N4, -0.15296486218853164},
        Grad4.{0.3239847771997537 /N4, 0.6740059517812944 /N4, -0.5794684678643381 /N4, 0.3239847771997537},
        Grad4.{-0.3239847771997537 /N4, 0.5794684678643381 /N4, -0.3239847771997537 /N4, -0.6740059517812944},
        Grad4.{-0.4004672082940195 /N4, 0.5029860367700724 /N4, 0.15296486218853164 /N4, -0.7504883828755602},
        Grad4.{0.15296486218853164 /N4, 0.5029860367700724 /N4, -0.4004672082940195 /N4, -0.7504883828755602},
        Grad4.{0.08164729285680945 /N4, 0.4553054119602712 /N4, 0.08164729285680945 /N4, -0.8828161875373585},
        Grad4.{-0.08164729285680945 /N4, 0.8828161875373585 /N4, -0.08164729285680945 /N4, -0.4553054119602712},
        Grad4.{-0.15296486218853164 /N4, 0.7504883828755602 /N4, 0.4004672082940195 /N4, -0.5029860367700724},
        Grad4.{0.4004672082940195 /N4, 0.7504883828755602 /N4, -0.15296486218853164 /N4, -0.5029860367700724},
        Grad4.{0.3239847771997537 /N4, 0.6740059517812944 /N4, 0.3239847771997537 /N4, -0.5794684678643381},
        Grad4.{0.5794684678643381 /N4, -0.6740059517812944 /N4, -0.3239847771997537 /N4, -0.3239847771997537},
        Grad4.{0.5029860367700724 /N4, -0.7504883828755602 /N4, -0.4004672082940195 /N4, 0.15296486218853164},
        Grad4.{0.5029860367700724 /N4, -0.7504883828755602 /N4, 0.15296486218853164 /N4, -0.4004672082940195},
        Grad4.{0.4553054119602712 /N4, -0.8828161875373585 /N4, 0.08164729285680945 /N4, 0.08164729285680945},
        Grad4.{0.8828161875373585 /N4, -0.4553054119602712 /N4, -0.08164729285680945 /N4, -0.08164729285680945},
        Grad4.{0.7504883828755602 /N4, -0.5029860367700724 /N4, -0.15296486218853164 /N4, 0.4004672082940195},
        Grad4.{0.7504883828755602 /N4, -0.5029860367700724 /N4, 0.4004672082940195 /N4, -0.15296486218853164},
        Grad4.{0.6740059517812944 /N4, -0.5794684678643381 /N4, 0.3239847771997537 /N4, 0.3239847771997537},
        Grad4.{0.5794684678643381 /N4, -0.3239847771997537 /N4, -0.6740059517812944 /N4, -0.3239847771997537},
        Grad4.{0.5029860367700724 /N4, -0.4004672082940195 /N4, -0.7504883828755602 /N4, 0.15296486218853164},
        Grad4.{0.5029860367700724 /N4, 0.15296486218853164 /N4, -0.7504883828755602 /N4, -0.4004672082940195},
        Grad4.{0.4553054119602712 /N4, 0.08164729285680945 /N4, -0.8828161875373585 /N4, 0.08164729285680945},
        Grad4.{0.8828161875373585 /N4, -0.08164729285680945 /N4, -0.4553054119602712 /N4, -0.08164729285680945},
        Grad4.{0.7504883828755602 /N4, -0.15296486218853164 /N4, -0.5029860367700724 /N4, 0.4004672082940195},
        Grad4.{0.7504883828755602 /N4, 0.4004672082940195 /N4, -0.5029860367700724 /N4, -0.15296486218853164},
        Grad4.{0.6740059517812944 /N4, 0.3239847771997537 /N4, -0.5794684678643381 /N4, 0.3239847771997537},
        Grad4.{0.5794684678643381 /N4, -0.3239847771997537 /N4, -0.3239847771997537 /N4, -0.6740059517812944},
        Grad4.{0.5029860367700724 /N4, -0.4004672082940195 /N4, 0.15296486218853164 /N4, -0.7504883828755602},
        Grad4.{0.5029860367700724 /N4, 0.15296486218853164 /N4, -0.4004672082940195 /N4, -0.7504883828755602},
        Grad4.{0.4553054119602712 /N4, 0.08164729285680945 /N4, 0.08164729285680945 /N4, -0.8828161875373585},
        Grad4.{0.8828161875373585 /N4, -0.08164729285680945 /N4, -0.08164729285680945 /N4, -0.4553054119602712},
        Grad4.{0.7504883828755602 /N4, -0.15296486218853164 /N4, 0.4004672082940195 /N4, -0.5029860367700724},
        Grad4.{0.7504883828755602 /N4, 0.4004672082940195 /N4, -0.15296486218853164 /N4, -0.5029860367700724},
        Grad4.{0.6740059517812944 /N4, 0.3239847771997537 /N4, 0.3239847771997537 /N4, -0.5794684678643381},
        Grad4.{0.03381941603233842 /N4, 0.03381941603233842 /N4, 0.03381941603233842 /N4, 0.9982828964265062},
        Grad4.{-0.044802370851755174 /N4, -0.044802370851755174 /N4, 0.508629699630796 /N4, 0.8586508742123365},
        Grad4.{-0.044802370851755174 /N4, 0.508629699630796 /N4, -0.044802370851755174 /N4, 0.8586508742123365},
        Grad4.{-0.12128480194602098 /N4, 0.4321472685365301 /N4, 0.4321472685365301 /N4, 0.7821684431180708},
        Grad4.{0.508629699630796 /N4, -0.044802370851755174 /N4, -0.044802370851755174 /N4, 0.8586508742123365},
        Grad4.{0.4321472685365301 /N4, -0.12128480194602098 /N4, 0.4321472685365301 /N4, 0.7821684431180708},
        Grad4.{0.4321472685365301 /N4, 0.4321472685365301 /N4, -0.12128480194602098 /N4, 0.7821684431180708},
        Grad4.{0.37968289875261624 /N4, 0.37968289875261624 /N4, 0.37968289875261624 /N4, 0.753341017856078},
        Grad4.{0.03381941603233842 /N4, 0.03381941603233842 /N4, 0.9982828964265062 /N4, 0.03381941603233842},
        Grad4.{-0.044802370851755174 /N4, 0.044802370851755174 /N4, 0.8586508742123365 /N4, 0.508629699630796},
        Grad4.{-0.044802370851755174 /N4, 0.508629699630796 /N4, 0.8586508742123365 /N4, -0.044802370851755174},
        Grad4.{-0.12128480194602098 /N4, 0.4321472685365301 /N4, 0.7821684431180708 /N4, 0.4321472685365301},
        Grad4.{0.508629699630796 /N4, -0.044802370851755174 /N4, 0.8586508742123365 /N4, -0.044802370851755174},
        Grad4.{0.4321472685365301 /N4, -0.12128480194602098 /N4, 0.7821684431180708 /N4, 0.4321472685365301},
        Grad4.{0.4321472685365301 /N4, 0.4321472685365301 /N4, 0.7821684431180708 /N4, -0.12128480194602098},
        Grad4.{0.37968289875261624 /N4, 0.37968289875261624 /N4, 0.753341017856078 /N4, 0.37968289875261624},
        Grad4.{0.03381941603233842 /N4, 0.9982828964265062 /N4, 0.03381941603233842 /N4, 0.03381941603233842},
        Grad4.{-0.044802370851755174 /N4, 0.8586508742123365 /N4, -0.044802370851755174 /N4, 0.508629699630796},
        Grad4.{-0.044802370851755174 /N4, 0.8586508742123365 /N4, 0.508629699630796 /N4, -0.044802370851755174},
        Grad4.{-0.12128480194602098 /N4, 0.7821684431180708 /N4, 0.4321472685365301 /N4, 0.4321472685365301},
        Grad4.{0.508629699630796 /N4, 0.8586508742123365 /N4, -0.044802370851755174 /N4, -0.044802370851755174},
        Grad4.{0.4321472685365301 /N4, 0.7821684431180708 /N4, -0.12128480194602098 /N4, 0.4321472685365301},
        Grad4.{0.4321472685365301 /N4, 0.7821684431180708 /N4, 0.4321472685365301 /N4, -0.12128480194602098},
        Grad4.{0.37968289875261624 /N4, 0.753341017856078 /N4, 0.37968289875261624 /N4, 0.37968289875261624},
        Grad4.{0.9982828964265062 /N4, 0.03381941603233842 /N4, 0.03381941603233842 /N4, 0.03381941603233842},
        Grad4.{0.8586508742123365 /N4, -0.044802370851755174 /N4, -0.044802370851755174 /N4, 0.508629699630796},
        Grad4.{0.8586508742123365 /N4, -0.044802370851755174 /N4, 0.508629699630796 /N4, -0.044802370851755174},
        Grad4.{0.7821684431180708 /N4, -0.12128480194602098 /N4, 0.4321472685365301 /N4, 0.4321472685365301},
        Grad4.{0.8586508742123365 /N4, 0.508629699630796 /N4, -0.044802370851755174 /N4, -0.044802370851755174},
        Grad4.{0.7821684431180708 /N4, 0.4321472685365301 /N4, -0.12128480194602098 /N4, 0.4321472685365301},
        Grad4.{0.7821684431180708 /N4, 0.4321472685365301 /N4, 0.4321472685365301 /N4, -0.12128480194602098},
        Grad4.{0.753341017856078 /N4, 0.37968289875261624 /N4, 0.37968289875261624 /N4, 0.37968289875261624},
    ];

    for i: 0..PSIZE-1 res[i] = arr[i%160];

    return res;    
}();
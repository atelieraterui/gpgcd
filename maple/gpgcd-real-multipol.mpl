# gpgcd-real-multipol
# GPGCD library for multiple polynomials with the real coefficients
# Maple (computer algebra system) source code, tested on Maple 13
# Copyright (c) 2009--2010, Akira Terui

with (LinearAlgebra);
with (PolynomialTools);
with (ArrayTools);
#with (ListTools); if ListTools is loaded, ListTools:-Transpose is
# called where LinearAlgebra:-Transpose should be called

gpgcd_real_multipol := proc (pol, var, deggcd, stopcriterion, numite_bound)

# gpgcd_real_multipol
# the main driver of the GPGCD method for multiple polynomials with
# real coefficients

# Inputs:
# pol: a list of input polynomials
# var: main variable
# deggcd: the degree of approximate gcd
# stopcriterion: stop criterion for itarations
# (the 2-norm of the update vector in iterations)
# numite_bound: upper bound for the number of iterations
# if the number of iterations > numite_bound, then the program stops
# itaratins

# Outputs:
# perturbation, perturbation2, gcd, newpols, cofactors, numite
# perturbation: the 2-norm of the perturbation terms added to f and g
# perturbation2: the square of 'perturbation'
# gcd: the approximate gcd calculated
# newpols: a list of perturbed polynomials
# cofactors: co-factors of polynomials in newpols satisfying
# pnew_1 * anew_1 + ... + pnew_m * anew_m = 0

local
    i, j,
    numpol, # the number of input polynomials
    d, # array of degrees of input polynomials
    pv, # array of coefficient vectors of input polynomials
    cv, # array of coefficient vectors of cofactors
    degv, # deggcd - 1, degree of polynomials vanish
    smat, # generalized subresultant matrix
    v0, # the initial vector for the modified Newton method
    v0dim, # dimension of v0
    vv, # new iterate in the iterations
    numite, # the number of iterations
    newpv, # array of calculated perturbed polynomials
    gv, # array of coefficient vectors of calculated gcds
    tmpv, # temporary vector for calculating residues
    tmp, # tmporary value for norms of residues
    res, # array of norm of residues
    resindex, # the index of gcd which gives the smallest residue
    gcd; # the final answer of gcd

    # Initial setup

    numpol := nops(pol);

    pv := Array(1 .. numpol);
    for i to numpol do
        pv[i] := CoefficientVector['reverse'](pol[i], var);
    end do;

#    print(pv);

    d := Array(1 .. numpol);
    for i to numpol do
        d[i] := degree(pol[i], var);
    end do;

#    print(degp);

    degv := deggcd - 1;
    smat := gensubresmat_vect(pv, numpol, d, degv);

#    print(smat);

    # Initialize the coefficient vector for itaration

    v0 := gpgcd_init_real_multipol(pv, smat, deggcd);
    v0dim := Dimension(v0);

#    print(v0);

# DEBUG();

#    pv, cv := vect2pols_real_multipol(v0, numpol, d, degv);

#    return v0, pv, cv, d;

    # Call the iteration routine

    numite, vv := modnewtoniterate_real_multipol(v0, numpol, d, degv,
                      stopcriterion, numite_bound);

    # numite is the number of itarations taken in the optimization
    # vv is the result of the optimization

    if numite = numite_bound then
        userinfo(1, gpgcd_real_multipol, `does not converge within iteration`, numite);
    else
        userinfo(1, gpgcd_real_multipol, `converges at iteration`, numite);
    end if;

#    return numite, vv;

    # Construct polynomials from vv satisfying P_1 * U_1 + ... P_k * U_k = 0

    v0 := vv[1 .. v0dim];
    newpv, cv := vect2pols_real_multipol(v0, numpol, d, degv);

    for i from 2 to numpol do
        cv[i] := -1.0 * cv[i];
    end do;

#    return newpv, cv;

    # Construct 'numpol' candidates of GCDs
    #
    # 1) Construct GCD gv[i] from newpv[i] (perturbed polynomials) and cv[i]
    # (calculated cofactors)

    gv := Array(1 .. numpol);

    for i to numpol do
        gv[i] := LeastSquares(coefvect2convmat(cv[i], deggcd), newpv[i]);
    end do;

    # 2) Compare which of gv[i] is the most appropriate as the common
    # divisor by calculating the norm of residue from the original
    # polynomials

    res := Array(1 .. numpol);

    for i to numpol do
        for j to numpol do
            tmpv := (coefvect2convmat(cv[j], deggcd) . gv[i]) - pv[j];
            res[i] := res[i] + (tmpv . tmpv);
        end do;
        res[i] := sqrt(res[i]);
    end do;

    userinfo(2, gpgcd_real_multipol, `residue: `, res);

    # calculate gv[i] which gives the minimum residue

    tmp := res[1];
    resindex := 1;
    for i from 2 to numpol do
        if (res[i] < tmp) then
            tmp := res[i];
            resindex := i;
        end if;
    end do;

    userinfo(2, gpgcd_real_multipol, `the index of the calculated GCDs which gives the minimum perturbation from the given polynomials is`, resindex);

    # Recalculate perturbed polynomials and the norm of relative
    # perturbations using gv[resindex]
    # relative perturbations of newpv[i] from pv[i] is put into res[i]

    for i to numpol do
        newpv[i] := coefvect2convmat(cv[i], deggcd) . gv[resindex];
        res[i] := Norm(newpv[i] - pv[i], 2) / Norm(pv[i], 2);
    end do;

    # Convert vectors (gv[resindex], newpv[i] and cv[i]) to polynomials

    userinfo(2, gpgcd_real_multipol, `newpv =`, newpv);

    newpv := map(coefvect2pol, newpv, var);
    cv := map(coefvect2pol, cv, var);
    gcd := coefvect2pol(gv[resindex], var);

    # Return the result

    return res, gcd, newpv, cv, numite;

#    return perturbation, perturbation2, gcd, newpols, cofactors, numite;
end proc;

gensubresmat := proc (p, numpol, d, var, degv)

# gensubresmat: create a generalized subresultant matrix of
# polynomials

# Inputs:
# p: list of polynomials
# numpol: the number of polynomials in pols
# d: list of degrees in pols
# var: main variable
# degv: the degree of the subresultant matrix

# Output:
# the (degv)-th generalized subresultant matrix of pols

local
    i, j, # indices
    c; # array of convolution matrices of pols
    # with (LinearAlgebra);

    c := Array(1 .. numpol - 1, 1 .. numpol);
    for i to numpol - 1 do
        c[i, 1] := convmat(p[i + 1], var, d[1] - 1 - degv);
        c[i, i + 1] := convmat(p[1], var, d[i + 1] - 1 - degv);
        for j from 2 to i do
            c[i, j] := Matrix(d[i + 1] + d[1] - degv, d[j] - degv);
        end do;
        for j from i + 2 to numpol do
            c[i, j] := Matrix(d[i + 1] + d[1] - degv, d[j] - degv);
        end do;
    end do;

    return convert(convert(c, listlist), Matrix);
end proc;

gensubresmat_vect := proc (pv, numpol, d, degv)

# gensubresmat_vect: create a subresultant matrix from the
# coefficient vectors of univariate polynomials

# Inputs:
# pv: Array of coefficient vectors of polynomials, s.t.
#        pv = [pv[1], ... , pv[n]],
#        pv[i] := [f[i]_m[i], ... , f[i]_0] satisfying
#        p[i] = f_[i]_m[i] x^m[i] + ... + f[i]_0 x^0,
# numpol: the number of input polynomials
# d: Array of degrees of pv
# degv: the degree of the subresulant matrix

# Output:
# smat: the (degp)-th subresultant matrix of f and g

local
    i, j, # indices
    c; # array of convolution matrices of pols

    # with (LinearAlgebra);

    c := Array(1 .. numpol - 1, 1 .. numpol);
    for i to numpol - 1 do
        c[i, 1] := coefvect2convmat(pv[i + 1], d[1] - 1 - degv);
        c[i, i + 1] := coefvect2convmat(pv[1], d[i + 1] - 1 - degv);
        for j from 2 to i do
            c[i, j] := Matrix(d[i + 1] + d[1] - degv, d[j] - degv);
        end do;
        for j from i + 2 to numpol do
            c[i, j] := Matrix(d[i + 1] + d[1] - degv, d[j] - degv);
        end do;
    end do;

    return convert(convert(c, listlist), Matrix);
end proc;

gpgcd_init_real_multipol := proc (pv, smat, deggcd)

# gpgcd_init_real_multipol: calculate the initial value for iterations

# Inputs:
# pv: Array of coefficient vectors of polynomials, s.t.
#        pv = [pv[1], ... , pv[n]],
#        pv[i] := [f[i]_m[i], ... , f[i]_0] satisfying
#        p[i] = f_[i]_m[i] x^m[i] + ... + f[i]_0 x^0,
# smat: generalized subresultant matrix of polynomials in pv
# deggcd: degree of gcd

# Output:
# result: an initial vector for iterations

local
    pvf, # pv in double
    pvfvect, # Vector of coefficients in pvf
    smatf, # smat in double
    cofactorcoef, # vector of coefficients of cofactors
    i;
    # with(PolynomialTools);
    # with(ArrayTools);
    # with(LinearAlgebra);

    pvf := evalf(pv);
    pvfvect := Vector['column'](convert(map(Transpose, pvf), list));

#    print(pvfvect);

    smatf := evalf(smat);
    cofactorcoef := -1.0 *
    SingularValues(smatf,output = ('Vt'))[ColumnDimension(smatf)];
    cofactorcoef := Vector['column'](cofactorcoef);

#DEBUG();

#    print(cofactorcoef);

    return convert(<pvfvect, cofactorcoef>, Vector['column'], datatype=float[8]);
end proc;

modnewtoniterate_real_multipol := proc (inipoint, numpol, d, degv,
                     stopcriterion, numite)

# modnewtoniterate_real_multipol: iteration routine of modified Newton method

# Inputs:
# inipoint: the coefficient vectors of initial polynomials
# numpol: the number of polynomials
# d: array of degrees of polynomials
# degv: (the degree of approximate gcd) - 1
# stopcriterion: stop criterion for itarations
# (the 2-norm of the update vector in iterations)
# numite: an upper bound for the number of itrations
# if the number of iterations > numite, then iteration stops

# Outputs:
# numite, dv0
# numite: the number of iterations taken
# dv0: the coefficient vector of calculated polynomials

local
    i, j,
    vv,
    vvdim,
    dv, # update vector for vector of coefficients in polynomials
    dv0,
    # solution vector of updates. Note: dv0 contains dv along with the
    # Lagrange multiplirs
    dvnorm; # the norm of dv

    vv := inipoint;
    vvdim := Dimension(inipoint);
    dv := Vector(vvdim);

#DEBUG();

    for i to numite do
        userinfo(2, modnewtoniterate_real_multipol, `Iteration`, i, print());

        # calculate the next iterate

        dv0 := modnewtonsolve_real_multipol(inipoint, vv, numpol, d, degv);
        dv := dv0[1 .. vvdim];
        dvnorm := Norm(dv, 2);
        vv := vv + dv;

        userinfo(2, modnewtoniterate_real_multipol, `vv =`, vv, print());
        userinfo(2, modnewtoniterate_real_multipol, `dv =`, dv, print());
        userinfo(2, modnewtoniterate_real_multipol, `Norm(dv) =`, dvnorm, print());

        if dvnorm < stopcriterion then
            userinfo(1, modnewtoniterate_real_multipol, `converges at iteration`, i);
            # return the cofficient vector
            # along with the Lagrange multipliers
            for j to vvdim do
                dv0[j] := vv[j];
            end do;
            return i, dv0;
        end if;
    end do;

    # if the value does not converge within the threshold number of
    # iterations, then return the result (the cofficient vector
    # along with the Lagrange multipliers) at that time

    for j to vvdim do
        dv0[j] := vv[j];
    end do;
    return numite, dv0;
end proc;

modnewtonsolve_real_multipol := proc (inipoint, v0, numpol, d, degv)

# modnewtonsolve_real_multipol: calculating ONE iteration of the modified
# Newton method

# Inputs:
# inipoint: the vector of coefficients of initial polynomials
# v0: the vector of coefficients of current polynomials
# numpol: the number of polynomials
# d: array of the degree of polynomials
# degv: (the degree of approximate GCD) - 1

# Output:
# LinearSolve(jmat, df): the output of LinearSolve
# note that the output is a Maple 'Vector'

local i,
    pv, # array of coefficient vectors of polynomials
    cv, # array of coefficient vectors of cofactors
    cvlist, # cv converted to flat list
    pdim, # dimension of coefficients in pols (not cofactors) appears in v0
    smat, # generalized subresultant matrix
    jmat, # Jacobian matrix
    coefMat, # coefficient matrix of the linear system
    constraintMat, # constraint matrix
    constraintVect, # constraint vector
    constraintValue, # constraint value as constraintMat . constraintVector
    inipointdiff, # difference of v0 from inipoint
    df; # the right-hand-side of the linear system to be solved

    # with(LinearAlgebra);

    # get coefficient vectors in each polynomial from v0

#DEBUG();

    pv, cv := vect2pols_real_multipol(v0, numpol, d, degv);
    cvlist := ListTools:-Flatten(convert(map(convert, cv, list), list));
    pdim := Dimension(v0) - nops(cvlist);

    # calculating Jacobian matrix of the constraint

    smat := gensubresmat_vect(pv, numpol, d, degv);
    jmat := jacobianmat_real_multipol(smat, cv, cvlist, numpol, d, degv);

    # calculating coefficient matrix of the linear system

    coefMat := < IdentityMatrix(ColumnDimension(jmat)), jmat |
                 Transpose(jmat), Matrix(RowDimension(jmat)) >;

    # calculating df, the right-hand-side of the linear system

    constraintMat, constraintVect := constraintMatVect_real_multipol(smat, cvlist);
    constraintValue := -1.0 * (constraintMat . constraintVect);

#    inipointdiff := Vector['column'](inipoint)[1 .. pdim]
#    - Vector['column'](v0)[1 .. pdim];
    inipointdiff := Vector['column'](inipoint) - Vector['column'](v0);

#    print(inipointdiff);

#    df := <inipointdiff, Vector['column'](nops(cvlist)), constraintValue>;
    df := <inipointdiff, constraintValue>;

#DEBUG();

    # solving the linear system

    return convert(LinearSolve(coefMat, df), Vector, datatype=float[8]);
end proc;

vect2pols_real_multipol := proc (v, numpol, d, degv)

# vect2pols_real_multipol: convert vector to arrays of coefficients in
# polynomials

# Inputs:
# v: the vector of coefficients of current polynomials
# numpol: the number of polynomials
# d: array of the degree of polynomials in v
# degv: (the degree of approximate GCD) - 1

# Outputs:
# pv: array of coefficient vectors of polynomials
# cv: array of coefficient vectors of cofactors
# note that the output is a Maple 'vector'

local
    i, j,
    pv, # see Outputs in the above
    cv; # see Outputs in the above

    pv := Array(1 .. numpol);
    cv := Array(1 .. numpol);

#DEBUG();

    j := 1;
    for i to numpol do
        pv[i] := Vector['column'](v[j .. j + d[i]]);
        j := j + d[i] + 1;
    end do;
#DEBUG();
    for i to numpol do
        cv[i] := Vector['column'](v[j .. j + d[i] - degv - 1]);
        j := j + d[i] - degv;
    end do;

    return pv, cv;
end proc;

jacobianmat_real_multipol := proc (smat, cv, cvlist, numpol, d, degv)

# jacobianmat_real_multipol: calculate the Jacobian matrix for a iteration

# Inputs:
# smat: generalized subresultant matrix
# cv: array of coefficient vectors of cofactors
# cvlist: cv converted to flat list
# numpol: the number of polynomials
# d: array of the degree of polynomials in pv
# degv: (the degree of approximate GCD) - 1

# Output
# jmat: the Jacobian matrix

local
    i, j, #indices
    jmat, # Jacobian matrix, to be returned
    jl, # the left block of jmat
    r, # the row dimension of the i-th row block
    offset1,
    offset2;

    # with (LinearAlgebra);

#DEBUG();

    jl := Array(1 .. numpol - 1, 1 .. numpol);
    for i to numpol - 1 do
        jl[i, 1] := coefvect2convmat(cv[i + 1], d[1]);
        jl[i, i + 1] := coefvect2convmat(cv[1], d[i + 1]);
        r := Dimension(cv[i + 1]) + d[1];
        for j from 2 to i do
            jl[i, j] := Matrix(r, d[j] + 1);
        end do;
        for j from i + 2 to numpol do
            jl[i, j] := Matrix(r, d[j] + 1);
        end do;
    end do;
    jl := Matrix(convert(jl, listlist));

#DEBUG();

    jmat :=  < Matrix(1, ColumnDimension(jl)), jl |
               2 * Matrix['row'](cvlist), smat >;

    return jmat
end proc;

constraintMatVect_real_multipol := proc (smat, cvlist)

# constraintMatVect_real_multipol: calculates constraint matrix and vector
# for caluclating
# "constraint value" = constraintMat . constraintVect

# Inputs:
# smat: generalized subresultant matrix
# cvlist: cv (array vectors of coefficients in cofactors) converted to
# flat list

local
    i, j,
    constraintMat, # constraint matrix
    constraintVect; # constraint vector

    # with(LinearAlgebra);

    constraintMat :=
    < Matrix[row](cvlist), smat |
    Matrix(1, 1, evalf(-1)), Matrix(RowDimension(smat), 1) >;

    constraintVect :=
    < Vector['column'](cvlist), Vector['column'](1, evalf(1)) >;

    return constraintMat, constraintVect;
end proc;

polynorm := proc (pol, var, norm)

# polynorm: calculate the polynomial norm

# Inputs:
# pol: a univariate polynomial
# var: the main variable
# order: the order of norm to calculate the (order)-the norm

# Output:
# VectorNorm(CoefficientVector(pol, var), norm):
# the (order)-th norm of pol w.r.t. variable 'var'

    # with (LinearAlgebra);
    # with (PolynomialTools);
    return VectorNorm(CoefficientVector(pol, var), norm);
end proc;

convmat := proc (pol, var, deg)

# convmat: construct the convoluiton matrix of a polynomial

# Inputs:
# pol: an input polynomial
# var: the main variable
# deg: the degree of convolution

# Output:
# coefmat: the convolution matrix

local i, j, polcoef, polcoefdim, coefmatrowdim, coefmatcoldim, coefmat;
    polcoef := CoefficientVector['reverse'](pol ,var);
    polcoefdim := Dimension(polcoef);
    coefmatrowdim := polcoefdim + deg;
    coefmatcoldim := deg + 1;
    coefmat := Matrix(coefmatrowdim, coefmatcoldim);
    for j to coefmatcoldim do
        for i to polcoefdim do
            coefmat[i + j - 1, j] := polcoef[i];
        end;
    end;
    return coefmat
end proc;

coefvect2convmat := proc (polcoef, deg)

# coefvect2convmat: construct the convoluiton matrix of a polynomial
# from the coefficient vector

# Inputs:
# polcoef: the coefficien vector of an input polynomial
# deg: the degree of convolution

# Output:
# coefmat: the convolution matrix

local
    i, j,
    polcoefdim, coefmatrowdim, coefmatcoldim, coefmat;
    polcoefdim := Dimension(polcoef);
    coefmatrowdim := polcoefdim + deg;
    coefmatcoldim := deg + 1;
    coefmat := Matrix(coefmatrowdim, coefmatcoldim);
    for j to coefmatcoldim do
        for i to polcoefdim do
            coefmat[i + j - 1, j] := polcoef[i];
        end;
    end;
    return coefmat
end proc;

coefvect2pol := proc(vect, var)

# coefvect2pol: construct a polynomial from the coefficient vector

# Inputs:
# vect: the coefficient vector of an input polynomial
# var: the main variable

# Output:
# the polynomial whose coefficiet vector is equal to 'vect'

local vectdim;
    vectdim := Dimension(vect);
    return add (vect[i] * var ^ (vectdim - i), i = 1 .. vectdim)
end proc;

# below are utility routines only used in experiments

convmat2 := proc(f, g, var, deg)
local i, j, fcoef, gcoef, fcoefdim, gcoefdim, coefmatrowdim, coefmatcoldim,
    coefmat;
    fcoef := CoefficientVector['reverse'](f, var);
    gcoef := CoefficientVector['reverse'](g, var);
    fcoefdim := Dimension(fcoef);
    gcoefdim := Dimension(gcoef);
    coefmatrowdim := fcoefdim + gcoefdim + 2 * deg;
    coefmatcoldim := deg + 1;
    coefmat := Matrix(coefmatrowdim, coefmatcoldim);
    for j to coefmatcoldim do
        for i to fcoefdim do
            coefmat[i + j - 1, j] := fcoef[i];
        end;
        for i to gcoefdim do
            coefmat[fcoefdim + deg + i + j - 1, j] := gcoef[i];
        end;
    end;
    return coefmat
end proc;

doublepol2coefvect := proc(f, g, var)
    local i, fcoef, gcoef, fcoefdim, gcoefdim, vect;
    fcoef := CoefficientVector['reverse'](f, var);
    gcoef := CoefficientVector['reverse'](g, var);
    fcoefdim := Dimension(fcoef);
    gcoefdim := Dimension(gcoef);
    vect := Vector(fcoefdim + gcoefdim);
    for i to fcoefdim do
        vect[i] := fcoef[i];
    end;
    for i to gcoefdim do
        vect[fcoefdim + i] := gcoef[i];
    end;
    return vect
end proc;

# gpgcd-complex
# GPGCD library for polynomials with the complex coefficients
# Maple (computer algebra system) source code, tested on Maple 12
# Copyright (c) 2009, Akira Terui

with (LinearAlgebra);
with (PolynomialTools);
with (ArrayTools);

gpgcd_complex := proc (f, g, var, deggcd, stopcriterion, numite_bound)

# gpgcd_complex
# the main driver of the GPGCD method for polynomials with complex coefficients

# Inputs:
# f, g: input polynomials
# var: main variable
# deggcd: the degree of approximate gcd
# stopcriterion: stop criterion for itarations
# (the 2-norm of the update vector in iterations)
# numite_bound: upper bound for the number of iterations
# if the number of iterations > numite_bound, then the program stops itaratins

# Outputs:
# perturbation, perturbation2, gcd, fnew, gnew, anew, bnew, numite
# perturbation: the 2-norm of the perturbation terms added to f and g
# perturbation2: the square of 'perturbation'
# gcd: the approximate gcd calculated
# fnew: a perturbed polynomial from f
# gnew: a perturbed polynomial from g
# anew, bnew: co-factors of fnew and gnew satisfying
# fnew * anew + gnew * bnew = 0

local
    i,
    degf,
    degg,
    degp,
    dega,
    degb,
    smat,
    v0,
    const0,
    vv,
    numite,
    fnew,
    gnew,
    anew,
    bnew,
    fnewCoef_re,
    gnewCoef_re,
    fnewCoef_im,
    gnewCoef_im,
    anewCoef_re,
    bnewCoef_re,
    anewCoef_im,
    bnewCoef_im,
    fnewcoef,
    gnewcoef,
    vvdim,
    anewConvMat_re,
    anewConvMat_im,
    anewconvmat,
    bnewConvMat_re,
    bnewConvMat_im,
    bnewconvmat,
    gcd,
    d1,
    d2,
    d1coef,
    d1Coef_re,
    d1Coef_im,
    d2coef,
    d2Coef_re,
    d2Coef_im,
    residue1,
    residue2,
    perturbation,
    perturbation2;

    # Initial setup

    degf := degree(f, var);
    degg := degree(g, var);
    degp := deggcd - 1;
    dega := degg - degp - 1;
    degb := degf - degp - 1;

    # Calculate initial value for iterations

    v0 := gpgcd_init_complex(f, g, var, deggcd);
    const0 := v0; # keep the initial value in const0

    # Call the iteration routine

    numite, vv := modnewtoniterate_complex(const0, v0, degp,
                                           stopcriterion, numite_bound);

    # numite is the number of itarations taken in the optimization
    # vv is the result of the optimization

    if numite = numite_bound then
        userinfo(1, gpgcd_complex, `does not converge within iteration`,
                 numite, print());
    else
        userinfo(1, gpgcd_complex, `converges at iteration`, numite, print());
    end if;

    # extract coefficient vectors of new polynomials from vv

    fnewCoef_re := vv[1]; # the real part of fnew
    gnewCoef_re := vv[2]; # the real part of gnew
    fnewCoef_im := vv[3]; # the imaginary part of fnew
    gnewCoef_im := vv[4]; # the imaginary part of gnew
    anewCoef_re := vv[5]; # the real part of anew
    bnewCoef_re := vv[6]; # the real part of bnew
    anewCoef_im := vv[7]; # the imaginary part of anew
    bnewCoef_im := vv[8]; # the imaginary part of bnew

    # Construct polynomials fnew, gnew, anew and bnew

    fnew := coefvect2pol_complex(fnewCoef_re, fnewCoef_im, var);
    gnew := coefvect2pol_complex(gnewCoef_re, gnewCoef_im, var);
    anew := coefvect2pol_complex(anewCoef_re, anewCoef_im, var);
    bnew := coefvect2pol_complex(bnewCoef_re, bnewCoef_im, var);
    bnew := - bnew;

    userinfo(2, gpgcd_complex, `fnew * anew - gnew * bnew =`,
             print(expand(fnew * anew - gnew * bnew)), print());

    # Construct GCDs by two ways
    #
    # 1) Construct d1, a GCD from bnew and fnew

    # construct convolution matrices from coefficient vectors of b

    bnewConvMat_re := coefvect2convmat(bnewCoef_re, deggcd);
    bnewConvMat_im := coefvect2convmat(bnewCoef_im, deggcd);

    # construct coefficient matrix for the least squares

    bnewconvmat := <bnewConvMat_re, bnewConvMat_im |
                    - bnewConvMat_im, bnewConvMat_re>;

    # construct right-hand-side vector for the least sequares

    fnewcoef := Vector['column']([fnewCoef_re, fnewCoef_im]);

    # calculate the least squares solution

    d1coef := LeastSquares(bnewconvmat, fnewcoef);

    # construct the coefficient vectors of a GCD

    d1Coef_re := d1coef[1 .. deggcd + 1];
    d1Coef_im := d1coef[deggcd + 2 .. 2 * (deggcd + 1)];

    # construct a GCD itself

    d1 := coefvect2pol_complex(d1Coef_re, d1Coef_im, var);
    userinfo(2, gpgcd_complex, `d1 =`, print(d1), print());

    # 2) Construct d2, a GCD from anew and gnew

    # construct convolution matrices from coefficient vectors of a

    anewConvMat_re := coefvect2convmat(anewCoef_re, deggcd);
    anewConvMat_im := coefvect2convmat(anewCoef_im, deggcd);

    # construct coefficient matrix for the least squares

    anewconvmat := <anewConvMat_re, anewConvMat_im |
                    - anewConvMat_im, anewConvMat_re>;

    # construct right-hand-side vector for the least sequares

    gnewcoef := Vector['column']([gnewCoef_re, gnewCoef_im]);

    # calculate the least squares solution

    d2coef := LeastSquares(anewconvmat, gnewcoef);

    # construct the coefficient vectors of a GCD

    d2Coef_re := d2coef[1 .. deggcd + 1];
    d2Coef_im := d2coef[deggcd + 2 .. 2 * (deggcd + 1)];

    # construct a GCD itself

    d2 := coefvect2pol_complex(d2Coef_re, d2Coef_im, var);
    if (Re(lcoeff(d1)) * Re(lcoeff(d2))) < 0 then
        d2 := - 1.0 * d2;
    end if;
    userinfo(2, gpgcd_complex, `d2 =`, print(d2), print());

    # control sign of anew and bnew

    if Re(lcoeff(d1 * bnew)) * Re(lcoeff(fnew)) < 0 then
        bnew := - 1.0 * bnew;
        userinfo(2, gpgcd_complex, `the sign of bnew has been changed:`,
                print(bnew), print());
    end if;
    if Re(lcoeff(d2 * anew)) * Re(lcoeff(gnew)) < 0 then
        anew := - 1.0 * anew;
        userinfo(2, gpgcd_complex, `the sign of anew has been changed:`,
                print(anew), print());
    end if;

    # 3) Compare whether d1 or d2 is more appropriate as the common divisor

    residue1 := sqrt(
        polynorm(expand(d1 * bnew - fnew), var, 2)^2 +
        polynorm(expand(d1 * anew - gnew), var, 2)^2
                    );
    residue2 := sqrt(
        polynorm(expand(d2 * bnew - fnew), var, 2)^2 +
        polynorm(expand(d2 * anew - gnew), var, 2)^2
                    );
    userinfo(2, gpgcd_complex,
             `sqrt(|d1 * bnew - fnew|^2 + |d1 * anew - gnew|^2) =`,
             residue1, print());
    userinfo(2, gpgcd_complex,
             `sqrt(|d2 * bnew - fnew|^2 + |d2 * anew - gnew|^2) =`,
             residue2, print());
    if residue1 <= residue2 then
        gcd := d1;
        userinfo(2, gpgcd_complex,
                 `the selected common divisor is d1`, print());
    else
        gcd := d2;
        userinfo(2, gpgcd_complex,
                 `the selected common divisor is d2`, print());
    end;

    # Re-define fnew and gnew

    fnew := expand(gcd * bnew);
    gnew := expand(gcd * anew);

    # Calulate perturbations

    perturbation2 := polynorm(fnew - f, var, 2) ^ 2 +
    polynorm(gnew - g, var, 2) ^ 2;
    perturbation := sqrt(perturbation2);

    # Return the result

    return perturbation, perturbation2, gcd, fnew, gnew, anew, bnew, numite;
end proc;

subresmat := proc (f, g, var, degp)

# subresmat: create a subresultant matrix of polynomials
# this function uses SylvesterMatrix in LinearAlgebra package

# Inputs:
# f, g: polynomials
# var: main variable
# degp: the degree of the subresultant matrix

# Output:
# the (degp)-th subresultant matrix of f and g

local degf, degg, smat;
    # with (LinearAlgebra);
    degf := degree(f, var);
    degg := degree(g, var);
    smat := SylvesterMatrix(f, g, var);
    return Transpose(SubMatrix(smat,[1 .. degg-degp,
                                     degg+1 .. degg+(degf-degp)],
                                [1 .. ColumnDimension(smat)-degp]))
end proc;

subresmat_vect := proc (fvect, gvect, degp)

# subresmat_vect: create a subresultant matrix from the
# coefficient vectors of univariate polynomials

# Inputs:
# fvect: coefficient vector of f
#        for f = f_m x^m + ... + f_0 x^0,
#        fvect := [f_m, ... , f_0]
# gvect: coefficient vector of g
#        for g = g_n x^n + ... + g_0 x^0,
#        gvect := [g_n, ... , g_0]
# degp: the degree of the subresulant matrix

# Output:
# smat: the (degp)-th subresultant matrix of f and g

local
    i, j, # indices
    degf, # degree of f
    degg, # degree of g
    smat; # subresultant matrix to be returned

    # with (LinearAlgebra);

    degf := Dimension(fvect) - 1;
    degg := Dimension(gvect) - 1;

    # definition of the subresultant matrix

    smat := Matrix(degf + degg - degp, degf + degg - 2 * degp);

    # making columns consisting of the coefficients in f

    for i to degg - degp do
        for j to degf + 1 do
            smat[i - 1 + j, i] := fvect[j];
        end do;
    end do;

    # making columns consisting of the coefficients in g

    for i to degf - degp do
        for j to degg + 1 do
            smat[i - 1 + j, degg - degp + i] := gvect[j];
        end do;
    end do;

    return smat;
end proc;

sylvester_like_matrix := proc (fvect, fcol, gvect, gcol)

# sylvester_like_matrix: create a sylvester-like matrix from the
# two vectors by specifying the number of column blocks

# Inputs:
# fvect: a column vector: fvect := [f_1, ... , f_m]
# fcol:  the number of columns to put fvect
# gvect: a column vector: gvect := [g_1, ... , g_n]
# gcol:  the number of columns to put gvect

# Assertion:
# the number of rows created with fvect and gvect
# must be coincide such that (m + fcol - 1) = (n + gcol - 1);
# otherwise stop with error message.

# Output:
# mat: the Sylvester-like matrix

local
    i, j, # indices
    fdim, # dimension of fvect
    gdim, # dimension of gvect
    fRowDim, # the row dimension of submatrix consisting of fvect
    gRowDim, # the row dimension of submatrix consisting of gvect
    mat;  # the matrix to be returned


    # with (LinearAlgebra);

    fdim := Dimension(fvect);
    gdim := Dimension(gvect);

    # definition of the matrix

    fRowDim := fdim + fcol - 1;
    gRowDim := gdim + gcol - 1;

    if (fRowDim <> gRowDim) then
        error "The numbers of rows do not match: fvect (%1) and gvect (%2)",
        fRowDim, gRowDim;
    end if;

    mat := Matrix(fRowDim, fcol + gcol);

    # making columns consisting of fvect

    for i to fcol do
        for j to fdim do
            mat[i - 1 + j, i] := fvect[j];
        end do;
    end do;

    # making columns consisting of gvect

    for i to gcol do
        for j to gdim do
            mat[i - 1 + j, fcol + i] := gvect[j];
        end do;
    end do;

    return mat;
end proc;

gpgcd_init_complex := proc (f, g, var, deggcd)

# gpgcd_init_complex: calculate the initial value for iterations

# Inputs:
# f, g: the given polynomials
# var: the main variable
# deggcd: degree of gcd

# Output:
# [fcoef_re, gcoef_re, fcoef_im, gcoef_im, acoef_re, bcoef_re,
# acoef_im, bcoef_im]:
# the list of coefficient vectors of the following polynomials
# fcoef_re: the real part of f
# gcoef_re: the real part of g
# fcoef_im: the imaginary part of f
# gcoef_im: the imaginary part of g
# acoef_re: the real part of a
# bcoef_re: the real part of b
# acoef_im: the imaginary part of a
# bcoef_im: the imaginary part of b

local
    i,    # index
    degf, # degree of f
    degg, # degree of g
    degp, # degree of p
    dega, # degree of a
    degb, # degree of b
    fcoef, # coefficient vector of f
    gcoef, # coefficient vector of g
    fcoef_re, # coefficient vector of the real part of f
    fcoef_im, # coefficient vector of the imaginary part of f
    gcoef_re, # coefficient vector of the real part of g
    gcoef_im, # coefficient vector of the imaginary part of g
    smat_re,  # subresultant matrix of the real part of f
    smat_im,  # subresultant matrix of the imaginary part of g
    smat,     # matrix for calculating the initial values
    abcoef,   # coefficients extracted from the right singular vector of smat
    acoef_re, # coefficient vector of the real part of a
    acoef_im, # coefficient vector of the imaginary part of a
    bcoef_re, # coefficient vector of the real part of b
    bcoef_im; # coefficient vector of the imaginary part of b

    # with(PolynomialTools);
    # with(ArrayTools);
    # with(LinearAlgebra);

    degf := degree(f,var);
    degg := degree(g,var);
    degp := deggcd - 1;
    dega := degg - degp - 1;
    degb := degf - degp - 1;

    fcoef := Transpose(CoefficientVector['reverse'](f,var));
    gcoef := Transpose(CoefficientVector['reverse'](g,var));

    fcoef := evalf(fcoef);
    gcoef := evalf(gcoef);

    fcoef_re := map(Re, fcoef);
    fcoef_im := map(Im, fcoef);
    gcoef_re := map(Re, gcoef);
    gcoef_im := map(Im, gcoef);

    smat_re := subresmat_vect(fcoef_re, gcoef_re, degp);
    smat_im := subresmat_vect(fcoef_im, gcoef_im, degp);

    smat := Matrix(2 * RowDimension(smat_re), 2 * ColumnDimension(smat_re),
                   [[smat_re, - smat_im], [smat_im, smat_re]]);

    # calculate abcoef from the SVD of smat

    abcoef := -1.0 *
    SingularValues(smat,output = ('Vt'))[ColumnDimension(smat)];

    # extract acoef_re, bcoef_re, acoef_im, bcoef_im from abcoef
    # s.t. abcoef = [acoef_re, bcoef_re, acoef_im, bcoef_im]

    acoef_re := Vector(dega + 1, abcoef[1 .. dega + 1]);
    bcoef_re := Vector(degb + 1, abcoef[dega + 2 .. dega + degb + 2]);
    acoef_im := Vector(dega + 1,
                       abcoef[dega + degb + 3 ..
                              2 * (dega + 1) + degb + 1]);
    bcoef_im := Vector(degb + 1,
                       abcoef[2 * (dega + 1) + degb + 2 ..
                              2 * (dega + degb + 2)]);

    # return the list of coefficient vectors as the initial values

    return [fcoef_re, gcoef_re, fcoef_im, gcoef_im, acoef_re, bcoef_re,
            acoef_im, bcoef_im];

end proc;

modnewtoniterate_complex := proc (const0list, initlist, degp,
                                  stopcriterion, numite)

# modnewtoniterate_complex: iteration routine of modified Newton method

# Inputs:
# const0list, initlist: the list of coefficient vectors of initial polynomials
# degp: (the degree of approximate gcd) - 1
# stopcriterion: stop criterion for itarations
# (the 2-norm of the update vector in iterations)
# numite: an upper bound for the number of itrations
# if the number of iterations > numite, then iteration stops

# Outputs:
# numite, v0list
# numite: the number of iterations taken
# v0list: the list of coefficient vectors of calculated polynomials

local
    i, j,        # indices
    vv,          # current value of the iteration
    dv,          # update value of the iteration
    dv0,         # original solution of the linear system
    dvnorm,      # 2-norm of dv
    fcoef_re,    # coefficient vector of the real part of f
    fcoef_im,    # coefficient vector of the imaginary part of f
    gcoef_re,    # coefficient vector of the real part of g
    gcoef_im,    # coefficient vector of the imaginary part of g
    acoef_re,    # coefficient vector of the real part of a
    acoef_im,    # coefficient vector of the imaginary part of a
    bcoef_re,    # coefficient vector of the real part of b
    bcoef_im,    # coefficient vector of the imaginary part of b
    fcoef_reDim, # dimension of fcoef_re
    gcoef_reDim, # dimension of gcoef_re
    acoef_reDim, # dimension of acoef_re
    bcoef_reDim, # dimension of bcoef_re
    v0list,      # list of calculated coefficient vectors
    offset;      # offset of index used in vv

    v0list := initlist;

    fcoef_reDim  := Dimension(v0list[1]);
    gcoef_reDim  := Dimension(v0list[2]);
    acoef_reDim  := Dimension(v0list[5]);
    bcoef_reDim  := Dimension(v0list[6]);

    vv := Vector['row'](initlist);

    # repeat the iteration for numite times

    for i to numite do
        userinfo(2, modnewtoniterate_complex, `Iteration`, i, print());

        # calculate the next iterate

        dv0 := convert(modnewtonsolve_complex(const0list, v0list, degp),
                       Vector['row']);

        # extract vector of the coefficients from the solution
        # (the rest values are the Lagrange multipliers)

        dv := dv0[1 ..
                  2 * (fcoef_reDim + gcoef_reDim + acoef_reDim + bcoef_reDim)];

        # DEBUG();

        vv := vv + dv;
        dvnorm := Norm(dv, 2);
        userinfo(2, modnewtoniterate_complex, `vv =`, vv, print());
        userinfo(2, modnewtoniterate_complex, `dv =`, dv, print());
        userinfo(2, modnewtoniterate_complex, `Norm(dv) =`, dvnorm, print());

        # Make up list of coefficient vectors

        fcoef_re := vv[1 .. fcoef_reDim];
        offset := fcoef_reDim;
        gcoef_re := vv[offset + 1 .. offset + gcoef_reDim];
        offset := offset + gcoef_reDim;
        fcoef_im := vv[offset + 1 .. offset + fcoef_reDim];
        offset := offset + fcoef_reDim;
        gcoef_im := vv[offset + 1 .. offset + gcoef_reDim];
        offset := offset + gcoef_reDim;
        acoef_re := vv[offset + 1 .. offset + acoef_reDim];
        offset := offset + acoef_reDim;
        bcoef_re := vv[offset + 1 .. offset + bcoef_reDim];
        offset := offset + bcoef_reDim;
        acoef_im := vv[offset + 1 .. offset + acoef_reDim];
        offset := offset + acoef_reDim;
        bcoef_im := vv[offset + 1 .. offset + bcoef_reDim];
        offset := offset + bcoef_reDim;

        v0list := [fcoef_re, gcoef_re, fcoef_im, gcoef_im,
                   acoef_re, bcoef_re, acoef_im, bcoef_im];

        userinfo(2, modnewtoniterate_complex, `v0list =`, v0list, print());

        # check the stopcriterion
        # if calculated value satisfies the criterion, the break the
        # roop and return the result

        if dvnorm < stopcriterion then
            userinfo(1, modnewtoniterate_complex,
                     `converges at iteration`, i, print());
            return i, v0list;
        end if;
    end do;

    # if the value does not converge within the threshold number of
    # iterations, then return the result at that time

    return numite, v0list;

end proc;

modnewtonsolve_complex := proc (const0list, v0list, degp)

# modnewtonsolve_complex: calculating ONE iteration of the modified
# Newton method

# Inputs:
# const0list: the list of coefficient vectors of initial
# polynomials, s.t.
# const0list = [f0coef_re, g0coef_re, f0coef_im, g0coef_im,
#               a0coef_re, b0coef_re, a0coef_im, b0coef_im]
#
# v0list: the list of coefficient vectors of current
# polynomials, s.t.
# v0list     = [fcoef_re, gcoef_re, fcoef_im, gcoef_im,
#               acoef_re, bcoef_re, acoef_im, bcoef_im]
# degp: (the degree of approximate GCD) - 1

# Output:
# LinearSolve(coefMat, df): the output of LinearSolve
# note that the output is a Maple 'vector'

local
    f0coef_re, # coefficient vector of the real part of f0
    f0coef_im, # coefficient vector of the imaginary part of f0
    g0coef_re, # coefficient vector of the real part of g0
    g0coef_im, # coefficient vector of the imaginary part of g0
    a0coef_re, # coefficient vector of the real part of a0
    a0coef_im, # coefficient vector of the imaginary part of a0
    b0coef_re, # coefficient vector of the real part of b0
    b0coef_im, # coefficient vector of the imaginary part of b0
    fcoef_re, # coefficient vector of the real part of f
    fcoef_im, # coefficient vector of the imaginary part of f
    gcoef_re, # coefficient vector of the real part of g
    gcoef_im, # coefficient vector of the imaginary part of g
    acoef_re, # coefficient vector of the real part of a
    acoef_im, # coefficient vector of the imaginary part of a
    bcoef_re, # coefficient vector of the real part of b
    bcoef_im, # coefficient vector of the imaginary part of b
    smat_re,  # subresultant matrix of the real part of f
    smat_im,  # subresultant matrix of the imaginary part of g
    aColDim,
    bColDim,
    abmat_re,
    abmat_im,
    jmat,
    coefMat,
    constraintMat,
    constraintVect,
    constraintValue,
    const0diff,
    df;

    # with(LinearAlgebra);

    f0coef_re := const0list[1];
    g0coef_re := const0list[2];
    f0coef_im := const0list[3];
    g0coef_im := const0list[4];
    a0coef_re := const0list[5];
    b0coef_re := const0list[6];
    a0coef_im := const0list[7];
    b0coef_im := const0list[8];

    fcoef_re  := v0list[1];
    gcoef_re  := v0list[2];
    fcoef_im  := v0list[3];
    gcoef_im  := v0list[4];
    acoef_re  := v0list[5];
    bcoef_re  := v0list[6];
    acoef_im  := v0list[7];
    bcoef_im  := v0list[8];

    smat_re := subresmat_vect(fcoef_re, gcoef_re, degp);
    smat_im := subresmat_vect(fcoef_im, gcoef_im, degp);

    aColDim := Dimension(fcoef_re);
    bColDim := Dimension(gcoef_re);

    abmat_re := sylvester_like_matrix(acoef_re, aColDim, bcoef_re, bColDim);
    abmat_im := sylvester_like_matrix(acoef_im, aColDim, bcoef_im, bColDim);

    # calculating Jacobian matrix of the constraint

    jmat := jacobianmat_complex(smat_re, smat_im, abmat_re, abmat_im,
                                acoef_re, bcoef_re, acoef_im, bcoef_im);

    # calculating coefficient matrix of the linear system

    coefMat := < IdentityMatrix(ColumnDimension(jmat)), jmat |
                 Transpose(jmat), Matrix(RowDimension(jmat)) >;

    # calculating df, the right-hand-side of the linear system

    constraintMat, constraintVect :=
    constraintMatVect_complex(smat_re, smat_im,
                              acoef_re, bcoef_re, acoef_im, bcoef_im);
    constraintValue := -1.0 * (constraintMat . constraintVect);

    const0diff := Vector['column'](const0list) - Vector['column'](v0list);
    df := <const0diff, constraintValue>;

    # solving the linear system

    return LinearSolve(coefMat, df)
end proc;

jacobianmat_complex := proc (smat_re, smat_im, abmat_re, abmat_im,
                           acoef_re, bcoef_re, acoef_im, bcoef_im)

# jacobianmat_complex: calculate the Jacobian matrix for a iteration

# Inputs:
# smat_re: the real part of the subresultant matrix of f and g
# smat_im: the imaginary part of the subresultant matrix of f and g
# abmat_re: the real part of the sylvester-like matrix of A and B
# abmat_im: the imaginary part of the sylvester-like matrix of A and B
# acoef_re: the real part of the coefficient vector of A
# bcoef_re: the real part of the coefficient vector of B
# acoef_im: the imaginary part of the coefficient vector of A
# bcoef_im: the imaginary part of the coefficient vector of B

# Output
# jmat: the Jacobian matrix s.t.
# jmat = [  0...0     0...0     2*t(acoef_re,bcoef_re,acoef_im,bcoef_im) ]
#        [ abmat_re -abmat_im     smat_re             -smat_im           ]
#        [ abmat_im  abmat_re     smat_im              smat_re           ]

local
    i, j, # indices
    jmat,
    smat_reRowDim,
    smat_reColDim,
    smat_imRowDim,
    smat_imColDim,
    abmat_reRowDim,
    abmat_reColDim,
    abmat_imRowDim,
    abmat_imColDim;

    # with (LinearAlgebra);

    smat_reRowDim, smat_reColDim := Dimension(smat_re);
    smat_imRowDim, smat_imColDim := Dimension(smat_im);
    abmat_reRowDim, abmat_reColDim := Dimension(abmat_re);
    abmat_imRowDim, abmat_imColDim := Dimension(abmat_im);

    jmat := <Matrix(1, ColumnDimension(abmat_re)), abmat_re, abmat_im |
             Matrix(1, ColumnDimension(abmat_re)), - abmat_im, abmat_re |
             Matrix([2*acoef_re, 2*bcoef_re]), smat_re, smat_im |
             Matrix([2*acoef_im, 2*bcoef_im]), - smat_im, smat_re>;

    return jmat
end proc;

constraintMatVect_complex := proc (smat_re, smat_im, acoef_re, bcoef_re,
                                   acoef_im, bcoef_im)
# constraintMatVect_complex: calculates constraint matrix and vector
# for caluclating
# "constraint value" = constraintMat . constraintVect

# Inputs:
# smat_re: the real part of the subresultant matrix of f and g
# smat_im: the imaginary part of the subresultant matrix of f and g
# acoef_re: the real part of the coefficient vector of A
# bcoef_re: the real part of the coefficient vector of B
# acoef_im: the imaginary part of the coefficient vector of A
# bcoef_im: the imaginary part of the coefficient vector of B

# Outputs:
# constraintMat, constraintVect: a matrix and vector, respectively, s.t.
# constraintMat = [ (acoef_re, bcoef_re, acoef_im, bcoef_im), -1 ]
#                 [        smat_re,           - smat_im,        0 ]
#                 [        smat_im,             smat_re,        0 ],
# constraintVect = t((acoef_re, bcoef_re, acoef_im, bcoef_im, 1))

local
    smatRow,
    constraintMat,
    constraintVect;

    # with(LinearAlgebra);

    smatRow := RowDimension(smat_re);

    constraintMat :=
    <Matrix([acoef_re, bcoef_re]), smat_re, smat_im |
     Matrix([acoef_im, bcoef_im]), - smat_im, smat_re |
     Matrix([evalf(-1)]), Matrix(smatRow,1), Matrix(smatRow,1)>;

    constraintVect :=
    <Vector['column'](acoef_re),
     Vector['column'](bcoef_re),
     Vector['column'](acoef_im),
     Vector['column'](bcoef_im),
     Vector['column'](1,evalf(1))>;

    return constraintMat, constraintVect;

end proc;

polynorm := proc (pol, var, order)

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
    return VectorNorm(CoefficientVector(pol, var), order);
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

coefvect2pol_complex := proc(vect_re, vect_im, var)

# coefvect2pol_complex: construct a polynomial from the coefficient
# vectors, with the real and the imaginary parts separately

# Inputs:
# vect_re: the coefficient vector of the real parts of an input polynomial
# vect_im: the coefficient vector of the imaginary parts of an input polynomial
# var: the main variable

# Output:
# the polynomial whose coefficiet vector is equal to
# 'vect_re' + 'vect_im' * sqrt(-1)

local vectdim;
    vectdim := Dimension(vect_re);
    return add ((vect_re[i] + vect_im[i] * I) * var ^ (vectdim - i),
                i = 1 .. vectdim)
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

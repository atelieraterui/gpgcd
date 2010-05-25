# gpgcd-real-multipol
# GPGCD library for multiple polynomials with the real coefficients
# Maple (computer algebra system) source code, tested on Maple 13
# Copyright (c) 2009--2010, Akira Terui

with (LinearAlgebra);
with (PolynomialTools);
with (ArrayTools);

gpgcd_real_multipol := proc (pols, var, deggcd, stopcriterion, numite_bound)

# gpgcd_real_multipol
# the main driver of the GPGCD method for multiple polynomials with
# real coefficients

# Inputs:
# pols: a list of input polynomials
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
    i,
    numpol, # the number of input polynomials
#    degf,
#    degg,
    degp, # array of degrees of input polynomials
    degv, # deggcd - 1, degree of polynomials vanish
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
    fnewcoef,
    gnewcoef,
    anewconvmat,
    bnewconvmat,
    gcd,
    d1,
    d2,
    d1coef,
    d2coef,
    residue1,
    residue2,
    perturbation,
    perturbation2;

    # Initial setup

    numpol := nops(pols);
    degp := Array(1 .. numpol);
    for i to numpol do
        degp[i] := degree(pols[i], var);
    end do;

    print(degp);

    degv := deggcd - 1;
#    dega := degg - degv - 1;
#    degb := degf - degv - 1;
    smat := subresmat(pols, numpol, degp, var, degv);

    # Initialize the coefficient vector for itaration

#    v0 := gpgcd_init_real(f, g, var, deggcd);
#    const0 := v0; # keep the initial value in const0

    # Call the iteration routine

#    numite, vv := modnewtoniterate_real(const0, v0, degf, degg, degp,
#                      stopcriterion, numite_bound);

    # numite is the number of itarations taken in the optimization
    # vv is the result of the optimization

#    if numite = numite_bound then
#        userinfo(1, gpgcd_real, `does not converge within iteration`, numite);
#    else
#        userinfo(1, gpgcd_real, `converges at iteration`, numite);
#    end if;

    # Construct polynomials from vv satisfying fnew * anew + gnew * bnew = 0

#    fnew := add(vv[i] * var ^ (degf + 1 - i), i = 1 .. degf + 1);
#    gnew := add(vv[degf + 1 + i] * var ^ (degg + 1 - i), i = 1 .. degg + 1);
#    anew := add(vv[degf + degg + 2 + i] * var ^ (dega + 1 - i),
#                i = 1 .. dega + 1);
#    bnew := add(vv[degf + degg + dega + 3 + i] * var ^ (degb + 1 - i),
#                i = 1 .. degb + 1);
#    bnew := - bnew;

    # Construct GCDs by two ways
    #
    # 1) Construct d1, a GCD from anew and fnew

#    bnewconvmat := convmat(bnew, var, deggcd);
#    fnewcoef := CoefficientVector['reverse'](fnew, var);
#    d1coef := LeastSquares(bnewconvmat, fnewcoef);
#    d1 := coefvect2pol(d1coef, var);
#    userinfo(2, gpgcd_real, `d1 = `, print(d1));

    # 2) Construct d2, a GCD from bnew and gnew

#    anewconvmat := convmat(anew, var, deggcd);
#    gnewcoef := CoefficientVector['reverse'](gnew, var);
#    d2coef := LeastSquares(anewconvmat, gnewcoef);
#    d2 := coefvect2pol(d2coef, var);
#    userinfo(2, gpgcd_real, `d2 = `, print(d2));

    # 3) Compare whether d1 or d2 is more appropriate as the common divisor

#    residue1 := sqrt(
#        polynorm(expand(d1 * bnew - fnew), var, 2)^2 +
#        polynorm(expand(d1 * anew - gnew), var, 2)^2
#                    );
#    residue2 := sqrt(
#        polynorm(expand(d2 * bnew - fnew), var, 2)^2 +
#        polynorm(expand(d2 * anew - gnew), var, 2)^2
#                    );
#    userinfo(2, gpgcd_real, `sqrt(|d1 * bnew - fnew|^2 + |d1 * anew - gnew|^2) =`,
#             residue1);
#    userinfo(2, gpgcd_real, `sqrt(|d2 * bnew - fnew|^2 + |d2 * anew - gnew|^2) =`,
#             residue2);
#    if residue1 <= residue2 then
#        gcd := d1;
#        userinfo(2, gpgcd_real, `the selected common divisor is d1`);
#    else
#        gcd := d2;
#        userinfo(2, gpgcd_real, `the selected common divisor is d2`);
#    end;

    # Re-define fnew and gnew

#    fnew := expand(gcd * bnew);
#    gnew := expand(gcd * anew);

    # Calulate perturbations

#    perturbation2 := polynorm(fnew - f, var, 2) ^ 2 +
#    polynorm(gnew - g, var, 2) ^ 2;
#    perturbation := sqrt(perturbation2);

    # Return the result

#    return perturbation, perturbation2, gcd, newpols, cofactors, numite;
end proc;

subresmat := proc (pols, numpol, degp, var, degv)

# subresmat: create a subresultant matrix of polynomials
# this function uses SylvesterMatrix in LinearAlgebra package

# Inputs:
# f, g: polynomials
# var: main variable
# degv: the degree of the subresultant matrix

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

gpgcd_init_real := proc (f, g, var, deggcd)

# gpgcd_init_complex: calculate the initial value for iterations

# Inputs:
# f, g: the given polynomials
# var: the main variable
# deggcd: degree of gcd

# Output:
# result: an initial vector for iterations

local degf, degg, degp, dega, degb, fcoef, gcoef, smat, abcoef, result, i;
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
    smat := subresmat(f, g, var, degp);
    abcoef := -1.0 *
    SingularValues(smat,output = ('Vt'))[ColumnDimension(smat)];
    result := Vector(Dimension(fcoef)+Dimension(gcoef)+Dimension(abcoef));
    for i to degf+1 do
        result[i] := fcoef[i];
    end do;
    for i to degg+1 do
        result[degf+1+i] := gcoef[i];
    end do;
    for i to Dimension(abcoef) do
        result[degf+degg+2+i] := abcoef[i];
    end do;
    return result
end proc;

modnewtoniterate_real := proc (const0, inipoint, degf, degg, degp,
                     stopcriterion, numite)

# modnewtoniterate_real: iteration routine of modified Newton method

# Inputs:
# const0, inipoint: the coefficient vectors of initial polynomials
# degf: the degree of f
# degg: the degree of g
# degp: (the degree of approximate gcd) - 1
# stopcriterion: stop criterion for itarations
# (the 2-norm of the update vector in iterations)
# numite: an upper bound for the number of itrations
# if the number of iterations > numite, then iteration stops

# Outputs:
# numite, dv0
# numite: the number of iterations taken
# dv0: the coefficient vector of calculated polynomials

local i, j, vv, vvdim, dv, dv0, dvnorm;
    vv := inipoint;
    vvdim := Dimension(inipoint);
    dv := Vector(vvdim);
    for i to numite do
        userinfo(2, modnewtoniterate_real, `Iteration`, i, print());

        # calculate the next iterate

        dv0 := modnewtonsolve_real(const0, vv, degf, degg, degp);

        for j to vvdim do
            dv[j] := dv0[j];
        end do;
        dvnorm := Norm(dv, 2);
        vv := vv + dv;

        userinfo(2, modnewtoniterate_real, `vv =`, vv, print());
        userinfo(2, modnewtoniterate_real, `dv =`, dv, print());
        userinfo(2, modnewtoniterate_real, `Norm(dv) =`, dvnorm, print());

        if dvnorm < stopcriterion then
            userinfo(1, modnewtoniterate_real, `converges at iteration`, i);
            for j to vvdim do
                dv0[j] := vv[j];
            end do;
            return i, dv0;
        end if;
    end do;

    # if the value does not converge within the threshold number of
    # iterations, then return the result at that time

    for j to vvdim do
        dv0[j] := vv[j];
    end do;
    return numite, dv0;
end proc;

modnewtonsolve_real := proc (const0, v0, degf, degg, degp)

# modnewtonsolve_real: calculating ONE iteration of the modified
# Newton method

# Inputs:
# const0: the vector of coefficients of initial polynomials
# v0list: the vector of coefficients of current polynomials
# degf: the degree of f
# degg: the degree of g
# degp: (the degree of approximate GCD) - 1

# Output:
# LinearSolve(jmat, df): the output of LinearSolve
# note that the output is a Maple 'vector'

local i, dega, degb, fvect, gvect, avect, bvect, jmat, constraintvalue,
    const0diff, dfdim, df;
    # with(LinearAlgebra);
    dega := degg - degp - 1;
    degb := degf - degp - 1;
    fvect := v0[1 .. degf + 1];
    gvect := v0[degf + 2 .. degf + degg + 2];
    avect := v0[degf + degg + 3 .. degf + degg + dega + 3];
    bvect := v0[degf + degg + dega + 4 .. degf + degg + dega + degb + 4];
    jmat := jacobianmat_real(fvect, gvect, avect, bvect, degp);
    constraintvalue := -1.0 * constrainteval_real(fvect, gvect, avect, bvect, degp);
    # const0diff := Vector(degf + degg + dega + degb + 4);
    const0diff := const0 - v0;
    dfdim := 3 * (degf + degg - degp + 1);
    df := Vector(dfdim);
#    for i to degf + degg + 2 do
    for i to degf + degg + dega + degb + 4 do
        df[i] := const0diff[i];
    end do;
    for i to degf + degg - degp + 1 do
        df[degf + degg + dega + degb + 4 + i] := constraintvalue[i];
    end do;
    userinfo(2, modnewtonsolve_real, `df =`, df, print());
    return LinearSolve(jmat, df)
end proc;

jacobianmat_real := proc (fvect, gvect, avect, bvect, degp)

# jacobianmat_real: calculate the Jacobian matrix for a iteration

# Inputs:
# fvect: the coefficient vector of f
# gvect: the coefficient vector of g
# avect: the coefficient vector of A
# bvect: the coefficient vector of B
# degp: (the degree of approximate GCD) - 1

# Output
# jmat: the Jacobian matrix

local
    i, j, #indices
    degf,
    degg,
    dega,
    degb,
    jmatdim,
    jmat,
    offset1,
    offset2;

    # with (LinearAlgebra);

    degf := Dimension(fvect) - 1;
    degg := Dimension(gvect) - 1;
    dega := Dimension(avect) - 1;
    degb := Dimension(bvect) - 1;
    jmatdim := 3*(degf + degg - degp + 1);
    jmat := Matrix(jmatdim);
    for i to degf+degg+dega+degb+4 do
        jmat[i,i] := 1.0;
    end do;
    offset2 := degf+degg+dega+degb+5;

    # Putting coefficients of a and b for the constraint of the norm
    for i to dega+1 do
        jmat[degf+degg+2+i, offset2] := -2.0*avect[i];
        jmat[offset2, degf+degg+2+i] := 2.0*avect[i];
    end do;
    for i to degb+1 do
        jmat[degf+degg+dega+3+i, offset2] := -2.0*bvect[i];
        jmat[offset2, degf+degg+dega+3+i] := 2.0*bvect[i];
    end do;

    # Putting the coefficients of a
    for i to degf+1 do
        for j to dega+1 do
            jmat[i, offset2+i+j-1] := -1.0*avect[j];
            jmat[offset2+i+j-1, i] := avect[j];
        end do;
    end do;

    # Putting the coefficients of b
    offset1 := degf + 1;
    for i to degg+1 do
        for j to degb+1 do
            jmat[offset1+i, offset2+i+j-1] := -1.0*bvect[j];
            jmat[offset2+i+j-1, offset1+i] := bvect[j];
        end do;
    end do;

    # Putting the coefficients of f
    offset1 := degf + degg + 2;
    for i to dega+1 do
        for j to degf+1 do
            jmat[offset1+i, offset2+i+j-1] := -1.0*fvect[j];
            jmat[offset2+i+j-1, offset1+i] := fvect[j];
        end do;
    end do;

    # Putting the coefficients of g
    offset1 := degf + degg + dega + 3;
    for i to degb+1 do
        for j to degg+1 do
            jmat[offset1+i, offset2+i+j-1] := -1.0*gvect[j];
            jmat[offset2+i+j-1, offset1+i] := gvect[j];
        end do;
    end do;
    return jmat
end proc;

constrainteval_real := proc (fvect, gvect, avect, bvect, degp)

# constrainteval_real: evaluate the constraint for current iterate

# Inputs:
# fvect: the coefficient vector of f
# gvect: the coefficient vector of g
# avect: the coefficient vector of A
# bvect: the coefficient vector of B
# degp: (the degree of approximate GCD) - 1

# Output
# valuevect: the value of constraint (in vector)

local i, degf, degg, dega, degb, smat, abvect, smatvect, valuevect;
    # with(LinearAlgebra);
    degf := Dimension(fvect) - 1;
    degg := Dimension(gvect) - 1;
    dega := Dimension(avect) - 1;
    degb := Dimension(bvect) - 1;
#    smat := constraintsmat_real(fvect, gvect, degp);
    smat := subresmat_vect(fvect, gvect, degp);
    abvect := Vector(dega + degb + 2);
    for i to dega + 1 do
        abvect[i] := avect[i];
    end do;
    for i to degb + 1 do
        abvect[dega + 1 + i] := bvect[i];
    end do;
    smatvect := smat . abvect;
    valuevect := Vector(degf + degg - degp + 1);
    for i to degf + degg - degp do
        valuevect[i + 1] := smatvect[i];
    end do;
    valuevect[1] := DotProduct(abvect, abvect) - 1.0;
    return valuevect
end proc;

constraintsmat_real := proc (fvect, gvect, degp)
local i, j, degf, degg, dega, degb, smat;
    # with(LinearAlgebra);
    degf := Dimension(fvect) - 1;
    degg := Dimension(gvect) - 1;
    dega := degg - degp - 1;
    degb := degf - degp - 1;
    smat := Matrix(degf + degg - degp, dega + degb + 2);
    for j to dega + 1 do
        for i to degf + 1 do
            smat[j + i - 1, j] := fvect[i];
        end do;
    end do;
    for j to degb + 1 do
        for i to degg + 1 do
            smat[j + i - 1, dega + 1 + j] := gvect[i];
        end do;
    end do;
    return smat
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

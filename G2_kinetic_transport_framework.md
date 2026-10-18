# G2-S gyron kinetic theory: geometry, collision invariants, and transport framework

**Status:** Research working document, 18 October 2026. Exact geometric and collision identities are distinguished from formal kinetic definitions and provisional hydrodynamic closures. No soliton existence or stability result is claimed.

## 1. Scope and assumptions

We consider identical axisymmetric rigid gyrons with fixed shape parameters \(\alpha,\beta,\chi\), mass \(m\), and length scale \(R\). Binary contacts are frictionless, smooth, and perfectly elastic. Particle orientation is \(Q\in SO(3)\), laboratory velocity \(v\), and intrinsic angular momentum \(L\). The gas is dilute when a local Boltzmann collision operator is used. For nonconvex shapes, contacts must be restricted to the actual boundary of the excluded relative-center region; locally normal-matched surface points are not necessarily accessible.

The principal objectives are to derive collision-controlled transport coefficients, preserve spin and spin-squared balance laws, and eventually test vortex-soliton ansatzes. Phenomenological constants must not be substituted for the shape-dependent collision calculation.

## 2. G2 shape and inertia

In body coordinates,

\[
X(t,\phi)=q(t)e_r(\phi)+\chi Rt e_3,\quad
q(t)=R\sqrt{1-t^2}\,p(t^2),\quad p(x)=1+\alpha x+\beta x^2,
\quad -1\le t\le1.
\]

Assume \(\chi>0\), \(p(x)>0\) for \(0\le x\le1\), and a regular non-self-intersecting surface. Uniform body density is assumed for the following inertia formulas. Define

\[
A_1=\int_{-1}^1(1-t^2)p(t^2)^2dt,\quad
A_2=\int_{-1}^1(1-t^2)^2p(t^2)^4dt,\quad
B_1=\int_{-1}^1t^2(1-t^2)p(t^2)^2dt.
\]

Then

\[
I_\parallel=\frac{mR^2}{2}\frac{A_2}{A_1},\qquad
I_\perp=mR^2\left(\chi^2\frac{B_1}{A_1}+\frac14\frac{A_2}{A_1}\right).
\]

For axis \(a=Qe_3\), \(I(Q)=I_\perp(1-aa^\top)+I_\parallel aa^\top\), and \(\omega=I(Q)^{-1}L\). The axial spin \(\ell=a\cdot L\) is conserved by torque-free rotation and by the smooth normal collision impulse.

## 3. Curvature and convexity

With \(x=t^2\),

\[
q'(t)=-\frac{Rt}{\sqrt{1-t^2}}h(x),\quad
h(x)=1-2\alpha+(3\alpha-4\beta)x+5\beta x^2,
\]
\[
q''(t)=-\frac{R}{(1-t^2)^{3/2}}W(x),
\]
\[
W(x)=1-2\alpha+(9\alpha-12\beta)x+(-6\alpha+33\beta)x^2-20\beta x^3.
\]

The outward body normal and surface area element are

\[
n_b=\frac{\chi R e_r-q'e_3}{\sqrt{\chi^2R^2+q'^2}},\qquad
dS=q\sqrt{\chi^2R^2+q'^2}\,dt\,d\phi.
\]

The principal curvatures (outward-positive convention for convex surfaces) are

\[
k_\phi=\frac{\chi R}{q\sqrt{\chi^2R^2+q'^2}},\qquad
k_m=-\frac{\chi Rq''}{(\chi^2R^2+q'^2)^{3/2}}.
\]

For regular profiles with \(p>0\), meridional convexity is equivalent to \(W(x)\ge0\) over \([0,1]\); strict convexity away from poles corresponds to \(W>0\). Since \(W(0)=1-2\alpha\), \(\alpha>1/2\) guarantees equatorial concavity, but \(\beta\) can cause concavity away from the equator even if \(\alpha\le1/2\). The convex/concave division changes contact geometry, not the collision conservation laws.

The squared normal-torque lever arm at a body surface point is

\[
H_G(t)=|X\times n_b|^2
=\frac{(q q'+\chi^2R^2t)^2}{\chi^2R^2+q'^2}.
\]

## 4. Collision law and invariants

At a physically admissible first contact of two particles, choose the unit normal \(n\) outward from particle 1 toward particle 2, with contact offsets \(r_1,r_2\). Define

\[
g=v_1+\omega_1\times r_1-v_2-\omega_2\times r_2,\quad
g_n=g\cdot n<0,
\]
\[
h_i=r_i\times n,\quad H=|h_1|^2+|h_2|^2,\quad
K=\frac2m+\frac{H}{I_\perp}.
\]

The elastic normal impulse on particle 1 is \(J=-2g_n n/K\), with updates

\[
\Delta v_1=J/m,\quad\Delta v_2=-J/m,\quad
\Delta L_1=r_1\times J,\quad\Delta L_2=-r_2\times J.
\]

The pair conserves total linear momentum, total orbital-plus-intrinsic angular momentum, and total kinetic energy. Since \(a_i\cdot(r_i\times n)=0\) for a surface of revolution, each \(\ell_i=a_i\cdot L_i\) is separately collision-invariant. Therefore all functions \(F(\ell)\) are collision invariants. This produces an infinite family of conserved axial-spin moments. Smooth spherical particles have \(H=0\), so their collisions do not exchange translation and rotation energy; their translational viscosity need not vanish.

For the pair's spin-squared increment,

\[
\Delta (|L_1|^2+|L_2|^2)
=-\frac{4g_n}{K}(L_1\cdot h_1-L_2\cdot h_2)
+\frac{4g_n^2}{K^2}H.
\]

Because \(I_\perp\) is fixed and each axial-spin square is invariant, the rotational-energy change obeys

\[
\Delta E_{\rm rot}=\frac{\Delta(|L_1|^2+|L_2|^2)}{2I_\perp}.
\]

## 5. Contact geometry for both convex and concave bodies

Let \(G\) be the solid body and \(E(Q_1,Q_2)=Q_1G\oplus(-Q_2G)\) the excluded relative-center set. This Minkowski sum is valid for concave as well as convex solids. The **physical collision surface** is the regular, externally accessible part of \(\partial E\), with the boundary normal and relative-velocity convention chosen consistently.

For a regular local contact branch \(b\), the surface points have opposed outward normals. Write \(B_{ib}\) for the signed inverse shape operators in a common tangent plane. On charts where the Gauss maps are invertible, the normal-parameterized area Jacobian is

\[
J_b=|\det(B_{1b}+B_{2b})|.
\]

The regular-branch measure can be expressed schematically as

\[
d\Sigma_G=\sum_b A_b J_b\,d\Omega_n,
\]

where \(A_b\) enforces nonoverlap, external accessibility, and unique counting of the same configuration-space boundary point. In a strictly convex regime, there is one branch and \(A_b=1\). Concave geometry can both introduce and occlude candidate branches; **it is not necessarily a strictly positive additive correction**.

For principal radii \(\rho_{m i},\rho_{\phi i}\) and relative tangent-frame angle \(\theta\),

\[
\begin{aligned}
\det(B_1+B_2)={}&\rho_{m1}\rho_{\phi1}+\rho_{m2}\rho_{\phi2}\\
&+\cos^2\theta(\rho_{m1}\rho_{\phi2}+\rho_{\phi1}\rho_{m2})\\
&+\sin^2\theta(\rho_{m1}\rho_{m2}+\rho_{\phi1}\rho_{\phi2}).
\end{aligned}
\]

At vanishing principal curvature, multiple simultaneous contacts, or singular normal charts, use the configuration-space boundary measure directly or a surface-coordinate/coarea construction. The branch formula above is **not** claimed globally valid without this qualification.

## 6. Kinetic equation and moment definitions

The one-particle distribution is \(f(x,v,Q,L,t)\), with rigid-body streaming and binary collisions:

\[
\partial_t f+v\cdot\nabla_x f+\mathcal R_G f=C_G[f,f].
\]

Write \(d\Gamma=dv\,dQ\,dL\), with normalized Haar measure for \(Q\). Define

\[
n=\int f\,d\Gamma,\quad \rho=mn,\quad
u=\frac1n\int vf\,d\Gamma,\quad c=v-u,
\]
\[
s=\int Lf\,d\Gamma,\quad
\mathsf H=\int I(Q)f\,d\Gamma,\quad
\Omega=\mathsf H^{-1}s,
\]
\[
N_2=\int|L|^2f\,d\Gamma,\quad
N_a=\int(a\cdot L)^2f\,d\Gamma,\quad
\mathsf A=\int(aa^\top-\mathsf1/3)f\,d\Gamma.
\]

The total rotational energy density is exactly

\[
E_{\rm rot}=\frac{N_2-N_a}{2I_\perp}+\frac{N_a}{2I_\parallel}.
\]

If \(e_s\) denotes rotational peculiar energy relative to \(\Omega\), then

\[
E_{\rm rot}=e_s+\frac12s\cdot\mathsf H^{-1}s.
\]

A scalar \(T_s\) defined by \(e_s=3nk_BT_s/2\) is only a parametrization of the total peculiar rotational energy. **Because axial spin is separately conserved, its temperature is not generally the temperature that equilibrates with translation.** For an unpolarized Maxwellian transverse sector define \(T_\perp\) by \(E_\perp=nk_BT_\perp\); the axial population may be arbitrary and remains separately conserved by collisions.

## 7. Eight schematic coupled moment balances

The following are balance-law structures, not a proven closed eight-field model. Fluxes include kinetic and, where finite particle size is retained, collisional-transfer parts. The tensor orientation and coherent-spin work terms must be derived consistently from angular-momentum and total-energy conservation.

\[
\begin{aligned}
&\partial_t\rho+\nabla\cdot(\rho u)=0,\\
&\rho D_tu=\nabla\cdot\sigma,\\
&\partial_ts+\nabla\cdot(u\otimes s+J_s)=\tau,\\
&\partial_te_t+\nabla\cdot(e_tu+q_t)=\mathcal P_t-Q_{t\perp},\\
&\partial_te_s+\nabla\cdot(e_su+q_s)=\mathcal P_s+Q_{t\perp},\\
&\partial_tN_2+\nabla\cdot(N_2u+j_2)=\mathcal C_2,\\
&\partial_tN_a+\nabla\cdot(N_au+j_a)=0,\\
&\partial_t\mathsf A+\nabla\cdot(u\mathsf A+J_A)=R_A.
\end{aligned}
\]

The exact relation between \(N_2,N_a,s\) and \(e_s\) means these eight lines contain redundant variables. In a final formulation one must choose independent moments and impose the energy identity, including corresponding flux/source identities. Conservation of every \(F(\ell)\) means that the full equilibrium manifold is generally infinite-dimensional; a finite closure requires an explicit restricted family of axial-spin distributions.

## 8. Provisional constitutive structure and coefficient inventory

For a near-isotropic, local, first-gradient closure, one may write

\[
\sigma=-p\mathsf1+2\mu_G D^\circ+\zeta_G(\nabla\cdot u)\mathsf1+\sigma^a+\sigma^N,
\]

where \(\sigma^a\) may involve rotational viscosity \(\mu_{r,G}\) coupling vorticity to intrinsic spin. Heat fluxes may be represented by a \(2\times2\) conductivity matrix \(\kappa_{tt},\kappa_{ts},\kappa_{st},\kappa_{ss}\). Spin, full spin-squared, axial spin-squared and nematic fluxes may contain self- and cross-diffusion terms. A linearized translation/transverse-rotation exchange coefficient \(\Lambda_G\) and orientational relaxation parameters may appear.

The previously listed **18 nominal scalar entries** are a bookkeeping inventory of a simplified closure, not 18 demonstrated independent physical parameters: three viscosities; nine self/cross diffusion entries; four thermal conductivities; two nominal relaxation entries. Tensorial anisotropy, Onsager relations, conservation identities and the infinite axial-spin invariants can alter this count. Orientational dynamics includes reversible free rotation, so a single nematic relaxation time is generally inadequate.

In general, every local coefficient depends on position and time **through the evolving state**:

\[
C_G(x,t)=C_G[\mathcal U(x,t);\alpha,\beta,\chi].
\]

Fixed \(\alpha,\beta,\chi\) do not imply constant viscosity or diffusivity. Near a homogeneous isotropic equilibrium, constant background values provide a useful first soliton-screening approximation. Farther from equilibrium the response may be anisotropic, nonlinear, nonlocal or memory-dependent.

## 9. Exact Gaussian reduction of energy exchange

For a spatially homogeneous unpolarized trial distribution with independent Gaussian translation and **transverse** rotation, let \(T_t\) and \(T_\perp\) be their temperatures. For a fixed contact geometry set

\[
w=(v_1-v_2)\cdot n,\quad
b=\omega_1\cdot h_1-\omega_2\cdot h_2,\quad g_n=w+b,
\]
\[
\kappa=H/I_\perp,\quad K=2/m+\kappa,\quad
V=\operatorname{Var}(g_n)=k_B(2T_t/m+\kappa T_\perp).
\]

The pair rotational-energy increment is

\[
\Delta E_{\rm rot}=-2g_nb/K+2\kappa g_n^2/K^2.
\]

Gaussian conditioning gives

\[
\boxed{\mathbb E[(-g_n)_+\Delta E_{\rm rot}]
=\frac{4k_B}{m}\sqrt{\frac{2V}{\pi}}\frac{\kappa}{K^2}(T_t-T_\perp).}
\]

This result is exact **within the specified Gaussian trial family** for each admissible regular contact, and shows exchange vanishes for \(T_t=T_\perp\) or \(H=0\).

Let \(dQ_1,dQ_2\) be normalized uniform orientation measures. The rotational collision-energy source per unit volume is

\[
Q_{\rm rot}=\frac{n^2}{2}\int dQ_1dQ_2d\Omega_n
\sum_b A_bJ_b\,\frac{4k_B}{m}\sqrt{\frac{2V_b}{\pi}}
\frac{\kappa_b}{K_b^2}(T_t-T_\perp).
\]

At common temperature \(T\), \(Q_{\rm rot}=\Lambda_G(T_t-T_\perp)+O((T_t-T_\perp)^2)\), with

\[
\boxed{\Lambda_G=\frac{2n^2k_B}{m}\sqrt{\frac{2k_BT}{\pi}}
\int dQ_1dQ_2d\Omega_n\sum_b A_bJ_b\frac{\kappa_b}{K_b^{3/2}}.}
\]

**Correct heat capacities:** The axial-spin population is collision-frozen. At fixed unpolarized axial distribution, translational energy is \(3nk_BT_t/2\) and transverse rotational energy is \(nk_BT_\perp\) (two rotational degrees of freedom). Thus, within this two-temperature homogeneous trial closure,

\[
\frac{d}{dt}(T_t-T_\perp)
=-\frac{5\Lambda_G}{3nk_B}(T_t-T_\perp),\qquad
\boxed{\tau_{t\perp}^{-1}=\frac{5\Lambda_G}{3nk_B}.}
\]

The earlier factor \(4\Lambda_G/(3nk_B)\), obtained by assigning three exchanging rotational degrees of freedom, is **not valid** when the axial-spin population is conserved. This distinction matters for the G2-S model.

Introduce \(\iota_\perp=I_\perp/(mR^2)\), \(\widehat H=H/R^2\), and \(\widehat J=J/R^2\). Define

\[
\mathcal F_{t\perp}(\alpha,\beta,\chi)
=\int dQ_1dQ_2d\Omega_n\sum_b A_b\widehat J_b
\frac{\widehat H_b/\iota_\perp}
{(2+\widehat H_b/\iota_\perp)^{3/2}}.
\]

Then

\[
\boxed{\Lambda_G=2n^2k_BR^2\sqrt{\frac{2k_BT}{\pi m}}\mathcal F_{t\perp},\qquad
\tau_{t\perp}^{-1}=\frac{10}{3}nR^2\sqrt{\frac{2k_BT}{\pi m}}\mathcal F_{t\perp}.}
\]

The integral has not yet been numerically evaluated, and the relaxation time is a **Gaussian-moment projected rate**, not necessarily an eigenvalue of the full linearized kinetic operator.

## 10. Collision matrices and remaining coefficients

For a stationary equilibrium \(f_0\), the symmetric binary collision bracket is formally

\[
[A,B]_G=\frac14\int f_{0,1}f_{0,2}(-g_n)_+
\Delta A\Delta B\,d\Sigma_G\,d\Gamma_1d\Gamma_2.
\]

The nullspace contains number, linear momentum, total energy, and all \(F(\ell)\). In a Galerkin approximation, the collision matrix \(K_{AB}=[\Phi_A,\Phi_B]_G\) and mass matrix \(M_{AB}=\langle\Phi_A\Phi_B\rangle_{f_0}\) yield projected relaxation eigenvalues from \(Kv=\lambda Mv\). Viscosities, diffusivities, and thermal conductivities require solving the corresponding projected linear-response equations with free rotational streaming and finite-size collisional transfer included where appropriate. Positivity, reciprocal relations in correct thermodynamic variables, and exact conservation constraints must be checked.

## 11. Soliton-search strategy

1. Derive the full collision and transport functions for fixed \(\alpha,\beta,\chi\), including accessible concave contact branches.
2. Compute near-isotropic reference coefficients and their dimensionless temperature/density scalings.
3. Screen steady axisymmetric and traveling/helical vortex ansatzes with constant background coefficients as a **controlled approximation**, not a fundamental assumption.
4. Continue candidate solutions with state-dependent coefficients \(C_G[\mathcal U(x,t)]\), retaining terms such as \(\nabla\mu_G\).
5. Evaluate anisotropic transport, conservation identities, stability spectra, and far-field behavior. A dissipative localized structure is not automatically a stable soliton.

## 12. Open mathematical tasks and verification checks

- Derive a globally valid nonconvex collision measure including first-contact accessibility, branch multiplicity, zero-curvature charts and any multiple-contact singularities.
- Verify the orientation/contact normalization against the hard-sphere Boltzmann cross section.
- Evaluate \(\mathcal F_{t\perp}(\alpha,\beta,\chi)\) by numerical quadrature and verify continuity through convexity boundaries where applicable.
- Determine how general axial-spin distributions affect equilibrium and projected transport.
- Derive full angular-momentum-consistent spin stresses, energy fluxes, and spin-squared flux identities before closing the moment system.
- Calculate shear viscosity, rotational viscosity, spin transport and spin-squared diffusion from the projected collision operator; test positivity and reciprocity.
- Avoid treating the 18 nominal scalar coefficients as independent fit parameters.

**Interpretation:** This is a shape-explicit kinetic/transport framework with one analytically reduced energy-exchange integral. It is not yet a completed calculation of all transport coefficients or a demonstrated soliton theory.

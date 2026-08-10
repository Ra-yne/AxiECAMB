# <a name="top"></a>AxiECAMB

Axion Effective-method in CAMB: a CAMB-based Boltzmann code with effective methods to accurately compute ultralight axion (ULA) observables in linear theory.

*[AxiECAMB, AxionCAMB, CAMB](#intro)
*[Getting Started](#basics)
*[Method](#physics)
*[Isocurvature](#isocurvature)
*[Warnings](#warnings)

----------------------------------------------------------------------
#### <a name="intro"></a>AxiECAMB, AxionCAMB, CAMB

AxiECAMB (Axion Effective-method in CAMB) executes the effective method for ultralight axions in https://arxiv.org/abs/2412.15192. The method is an extension of (Passaglia and Hu 2022, https://arxiv.org/abs/2201.10238) to treat corrections for metric perturbations and the CMB, so that linear observables such as the CMB and matter power spectra can achieve subpercent accuracy for ULA parameter estimation.

AxiECAMB heavily modifies the existing ultralight axion Boltzmann code axionCAMB (https://github.com/dgrin1/axionCAMB), which itself is based on CAMB (http://camb.info/, https://github.com/cmbant/CAMB) version Nov 2013. 

When using AxiECAMB, please cite https://arxiv.org/abs/2412.15192. For its isocurvature features, please cite https://arxiv.org/abs/XXXX.XXXXX.

----------------------------------------------------------------------
#### <a name="basics"></a>Getting Started

If cloned to local for the first time, run make to compile the code. You can modify Makefile according to your compiler. For basic usage, modify params.ini according to desired cosmology and run ./camb params.ini. 

The main addition to params.ini due to ULAs, as compared to LCDM, are 

 - m_ax: the ULA mass in units of eV;
 - omaxh2: the ULA abundance;
 - axfrac: the ULA fraction -- the code determines internally whether m_ax/H_0 >= 10 which assigns this fraction to that of the dark matter or dark energy (see "Method" section below)

The parameter use_axfrac specifies whether to use the ULA fraction parametrization.

Parameters new in v1.1:

| parameter | default | meaning |
|---|---|---|
| movH_switch | 10 | The user-input m/H_* at which the exact Klein-Gordon evolution switches to the effective fluid approximation. The baseline/minimum value is 10, and smaller values are refused. See "Method". |
| iso_input_mode | 2 | selects how the isocurvature amplitude is specified: 1 = A_iso, 2 = r_a_iso, 3 = Hinf. See "Isocurvature". |
| A_iso | *required in mode 1* | isocurvature amplitude at k0. |
| r_a_iso | *required in mode 2* | isocurvature amplitude relative to A_s. |
| Hinf | *required in mode 3* | log10(H_inflation/GeV). |
| n_a_iso | 0 | isocurvature tilt, tensor convention (0 = scale invariant). |
| iso_consistency | F | impose n_a_iso = -r_a_iso/8. |

The three amplitude parameters are read only in their own mode, so only one of them need be present.

----------------------------------------------------------------------
#### <a name="physics"></a>Method

AxiECAMB adopts ULAs with a quadratic potential and operates in synchronous gauge. It solves the exact Klein-Gordon (KG) equation from the start time to m/H = movH_switch, and switches to an effective fluid approximation that accurately tracks the time-average of the fast oscillations in the ULA KG system. For light axion masses where m/H~10 is around equality, m/H is raised above 10 and the switch is prevented from occuring around recombination.  For the lightest ULAs m/H0<10, which is around 1.4e-32 eV, KG is solved to the present day.

AxiECAMB re-solves movH_switch so that the Klein-Gordon-to-fluid transition does not fall during recombination, and for masses whose switch is near recombination, AxiECAMB adjusts the switch for the correct phase, so the realized switch can differ from the requested value (see Appendices B and E of https://arxiv.org/abs/2412.15192). Both the requested and the realized value are reported at feedback_level > 0.

Larger movH_switch values track more oscillation cycles exactly and only improve accuracy when paired with a higher accuracy_boost. The verified setting is movH_switch = 10 with accuracy_boost = 1, which is sufficient for data analysis; for other combinations, check convergence of your target observables.

In addition to these methods, AxiECAMB fixes numerous bugs and inaccuracies in AxionCAMB such that it is typically ~50% faster. 

With default settings for accuracy and k<1e3 Mpc-1, AxiECAMB produces accurate results for m<~1e-18eV. For heavier masses, ULA effects only appear at higher k, and fitting functions for matter power spectra specified in https://arxiv.org/abs/2201.10238 is more efficient. 

Our threshold for dark matter vs dark energy like ULAs is whether the mass exceeds m_ETA = 10 H0(fiducial) = 1.44e-32 eV. From v1.1 this boundary is fixed in eV rather than recomputed from m/H0 = 10 for each sampled cosmology, so that the parameter definition does not change when H0 is varied in parameter estimation. Above it, the axfrac (fraction of ULA-abundance) parameter represents the ULA fraction in dark matter, and ULA contributes to the total matter transfer function. Below it, axfrac refers to the ULA fraction in dark energy instead, and the matter transfer function excludes the ULA. In this case KG is solved to the present day and present day cosmological parameters are specified as their instantaneous z=0 values not their time average.

----------------------------------------------------------------------
#### <a name="isocurvature"></a>Isocurvature

New in v1.1. Set axion_isocurvature = T to compute axion isocurvature spectra; it is F by default, and every parameter in this section is read only when it is T.

The primordial isocurvature spectrum is

    Delta^2_iso(k) = A_iso (k/k0)^n_a_iso,     k0 = pivot_scalar

specified in one of three ways, selected by the integer iso_input_mode:

 1. A_iso -- the amplitude at k0 directly, independent of the inflationary origin;
 2. r_a_iso (default) -- A_iso is derived from r_a_iso, A_s and the initial field value phi_ini;
 3. Hinf -- log10(H_inflation/GeV). r_a_iso is derived for a massless spectator field, then treated as in mode 2.

Each mode requires only its own amplitude and ignores the other two. Setting iso_consistency = T imposes the single-field consistency relation n_a_iso = -r_a_iso/8 in every mode; if F, n_a_iso is taken from params.ini. n_a_iso follows the tensor convention: n_a_iso = 0 is scale invariant, corresponding to n = 1 in the CAMB scalar convention.

At feedback_level > 0, AxiECAMB reports the derived initial field value phi_ini/Mpl together with the amplitude expressed in all three parameterizations.

**Adiabatic and isocurvature spectra are computed in separate runs and are NOT summed by AxiECAMB.** Run each separately and sum them externally. 

The fiducial-template Cl interpolation is forced off (use_spline_template = F) for isocurvature runs, which must never be interpolated against an adiabatic template.

Isocurvature transfer functions are written to transfer_out.dat in the same CAMB column convention as the adiabatic ones, i.e. delta_i/[k^2 delta_ax(tau_ini)]. To obtain T_i = delta_i/delta_ax(tau_ini) as defined in the accompanying paper, multiply by k^2 with k in Mpc^-1, i.e. k = (column 1) * h, since column 1 is k/h.

----------------------------------------------------------------------
#### <a name="warnings"></a>Warnings

Tensor spectra are not supported in this version and are deferred to AxiECAMB 2.0; setting get_tensor_cls = T is refused with an error.

The ULA GrowthRate subroutine pertaining to SZ analysis (and other measures of the growth rate, e.g. redshift-space distortions) is currently disabled.

If the matter power spectra or transfer function is needed for z>0, one needs to be careful regarding whether z is before or after the switch point to fluid approximation, given that they are defined differently before and after the switch. The code will display a warning sign if a z>0 is requested that is before the switch.

The non-linear mode of the code is inherited from axionCAMB and not extensively tested. Turning on the non-linear mode will trigger a warning message. Halofit has not been recalibrated for axions, thus where the effective slope leaves the range Halofit was fitted over, the code falls back to linear theory for that scale.

accurateBB in params.ini should be assigned as T if accurate BB polarization power spectra are needed.

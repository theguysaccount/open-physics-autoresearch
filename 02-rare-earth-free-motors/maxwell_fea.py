"""
maxwell_fea.py — Exact 2D Nonlinear Magnetostatic FEA (Newton-Raphson + Exact
Sliding-Band Rotation + Arkkio Torque + Virtual Work + Bertotti Iron Loss +
AC Proximity Copper Loss + Cold/Hot Demagnetization) AND Exact 2D Plane-Stress
Centrifugal Elasticity FEA with General Rotational Periodicity across 4/6/8-pole topologies.
"""

from __future__ import annotations
from typing import Dict, Tuple, Optional
import math
import numpy as np
from numba import njit
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from materials import (
    MU0, NU0, eval_nu_dnu_steel,
    STEEL_RHO, STEEL_E_MOD, STEEL_POISSON, STEEL_SIGMA_Y,
    KH_VIRGIN, KH_CUT_MULT, KE_STEEL, KA_STEEL,
    CU_RHO_120C, CU_WIRE_DIAM, MAGNET_DB,
)
from motor_geometry import (
    L_STK, K_FILL_MAX, N_THETA_GAP,
    REG_AIR_BARRIER,
    REG_ROTOR_STEEL, REG_ROTOR_STEEL_CUT, REG_STATOR_STEEL, REG_STATOR_STEEL_CUT,
    REG_MAGNET_BASE, REG_SLOT_BASE,
    MotorDesign, MotorMeshBundle,
)


# Electrical offset of rotor d-axis (theta_e = 90 deg) relative to Phase A magnetic axis (-60 deg)
# Universal for ALL integer-slot full-pitch windings [+A]*q + [-C]*q + [+B]*q across 4p, 6p, 8p!
GAMMA_E0_RAD = math.radians(150.0)


@njit(cache=True)
def _precompute_tri_geometry(p: np.ndarray, t: np.ndarray):
    E = t.shape[1]
    areas = np.empty(E, dtype=np.float64)
    gx = np.empty((E, 3), dtype=np.float64)
    gy = np.empty((E, 3), dtype=np.float64)
    xc = np.empty(E, dtype=np.float64)
    yc = np.empty(E, dtype=np.float64)
    for e in range(E):
        i0, i1, i2 = t[0, e], t[1, e], t[2, e]
        x0, y0 = p[0, i0], p[1, i0]
        x1, y1 = p[0, i1], p[1, i1]
        x2, y2 = p[0, i2], p[1, i2]
        det = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
        areas[e] = 0.5 * det
        inv_det = 1.0 / det
        gx[e, 0] = (y1 - y2) * inv_det
        gx[e, 1] = (y2 - y0) * inv_det
        gx[e, 2] = (y0 - y1) * inv_det
        gy[e, 0] = (x2 - x1) * inv_det
        gy[e, 1] = (x0 - x2) * inv_det
        gy[e, 2] = (x1 - x0) * inv_det
        xc[e] = (x0 + x1 + x2) / 3.0
        yc[e] = (y0 + y1 + y2) / 3.0
    return areas, gx, gy, xc, yc


@njit(cache=True)
def _assemble_nonlinear_maxwell(
    A_free: np.ndarray,
    t: np.ndarray,
    areas: np.ndarray,
    gx: np.ndarray,
    gy: np.ndarray,
    mat_code: np.ndarray,
    nu_lin: np.ndarray,
    Brx: np.ndarray,
    Bry: np.ndarray,
    Jz: np.ndarray,
    dof_map: np.ndarray,
    dof_sign: np.ndarray,
    n_free: int,
):
    E = t.shape[1]
    rows = np.empty(9 * E, dtype=np.int64)
    cols = np.empty(9 * E, dtype=np.int64)
    vals = np.empty(9 * E, dtype=np.float64)
    res = np.zeros(n_free, dtype=np.float64)
    Bx_arr = np.empty(E, dtype=np.float64)
    By_arr = np.empty(E, dtype=np.float64)
    nnz = 0

    for e in range(E):
        area = areas[e]
        n0, n1, n2 = t[0, e], t[1, e], t[2, e]
        d0, d1, d2 = dof_map[n0], dof_map[n1], dof_map[n2]
        s0, s1, s2 = dof_sign[n0], dof_sign[n1], dof_sign[n2]

        a0 = s0 * A_free[d0] if d0 >= 0 else 0.0
        a1 = s1 * A_free[d1] if d1 >= 0 else 0.0
        a2 = s2 * A_free[d2] if d2 >= 0 else 0.0

        gx0, gx1, gx2 = gx[e, 0], gx[e, 1], gx[e, 2]
        gy0, gy1, gy2 = gy[e, 0], gy[e, 1], gy[e, 2]

        dAx = a0 * gx0 + a1 * gx1 + a2 * gx2
        dAy = a0 * gy0 + a1 * gy1 + a2 * gy2
        Bx_arr[e] = dAy
        By_arr[e] = -dAx
        B2 = dAx * dAx + dAy * dAy

        mc = mat_code[e]
        if mc == 1:
            nu, dnu = eval_nu_dnu_steel(B2, False)
        elif mc == 2:
            nu, dnu = eval_nu_dnu_steel(B2, True)
        else:
            nu = nu_lin[e]
            dnu = 0.0

        gd0 = gx0 * dAx + gy0 * dAy
        gd1 = gx1 * dAx + gy1 * dAy
        gd2 = gx2 * dAx + gy2 * dAy

        f_j = (area / 3.0) * Jz[e]
        if mc == 3:
            brx = Brx[e]
            bry = Bry[e]
            fm0 = area * nu * (brx * gy0 - bry * gx0)
            fm1 = area * nu * (brx * gy1 - bry * gx1)
            fm2 = area * nu * (brx * gy2 - bry * gx2)
        else:
            fm0 = 0.0
            fm1 = 0.0
            fm2 = 0.0

        r0 = area * nu * gd0 - (f_j + fm0)
        r1 = area * nu * gd1 - (f_j + fm1)
        r2 = area * nu * gd2 - (f_j + fm2)

        if d0 >= 0:
            res[d0] += s0 * r0
        if d1 >= 0:
            res[d1] += s1 * r1
        if d2 >= 0:
            res[d2] += s2 * r2

        ds = (d0, d1, d2)
        ss = (s0, s1, s2)
        gxs = (gx0, gx1, gx2)
        gys = (gy0, gy1, gy2)
        gds = (gd0, gd1, gd2)
        c_dnu = 2.0 * area * dnu
        c_nu = area * nu

        for i in range(3):
            di = ds[i]
            if di < 0:
                continue
            si = ss[i]
            gxi = gxs[i]
            gyi = gys[i]
            gdi = gds[i]
            for j in range(3):
                dj = ds[j]
                if dj < 0:
                    continue
                kij = c_nu * (gxi * gxs[j] + gyi * gys[j]) + c_dnu * gdi * gds[j]
                rows[nnz] = di
                cols[nnz] = dj
                vals[nnz] = si * ss[j] * kij
                nnz += 1

    return rows[:nnz], cols[:nnz], vals[:nnz], res, Bx_arr, By_arr


@njit(cache=True)
def _compute_coenergy_1pole(
    Bx: np.ndarray, By: np.ndarray, areas: np.ndarray,
    mat_code: np.ndarray, nu_lin: np.ndarray, Brx: np.ndarray, Bry: np.ndarray
) -> float:
    gl_x = np.array([0.06943184420297371, 0.33000947820757187,
                     0.6699905217924281,  0.9305681557970262])
    gl_w = np.array([0.1739274225687269,  0.3260725774312731,
                     0.3260725774312731,  0.1739274225687269])
    W_mag = 0.0
    for e in range(len(areas)):
        bx, by = Bx[e], By[e]
        B2 = bx * bx + by * by
        mc = mat_code[e]
        if mc == 0:
            wm = 0.5 * nu_lin[e] * B2
        elif mc == 3:
            dbx = bx - Brx[e]
            dby = by - Bry[e]
            wm = 0.5 * nu_lin[e] * (dbx * dbx + dby * dby)
        else:
            is_cut = (mc == 2)
            acc = 0.0
            for q in range(4):
                u = gl_x[q] * B2
                nu_q, _ = eval_nu_dnu_steel(u, is_cut)
                acc += gl_w[q] * nu_q
            wm = 0.5 * B2 * acc
        W_mag += wm * areas[e]
    return W_mag


class MotorFEASolver:
    """Reusable nonlinear 2D Maxwell FEA solver supporting 4p, 6p, and 8p topologies."""

    def __init__(self, design: MotorDesign, mesh: MotorMeshBundle, mag_temp_mode: str = "HOT_140C"):
        self.design = design
        self.mesh = mesh
        self.mag_temp_mode = mag_temp_mode
        self.p_pairs = mesh.p_pairs
        self.n_poles = 2 * mesh.p_pairs
        self.q_spp = mesh.q_spp
        self.n_spp = mesh.slots_per_pole
        self.th_pole = mesh.theta_pole

        # Build 1-pole winding phase & polarity arrays: [+A]*q + [-C]*q + [+B]*q
        self.slot_phase = np.array([0] * self.q_spp + [2] * self.q_spp + [1] * self.q_spp, dtype=np.int32)
        self.slot_pol = np.array([+1.0] * self.q_spp + [-1.0] * self.q_spp + [+1.0] * self.q_spp, dtype=np.float64)

        self.n_rot_nodes = mesh.rot_p.shape[1]
        self.n_sta_nodes = mesh.sta_p.shape[1]
        self.n_tot_nodes = self.n_rot_nodes + self.n_sta_nodes
        self.n_rot_elems = mesh.rot_t.shape[1]
        self.n_sta_elems = mesh.sta_t.shape[1]
        self.n_tot_elems = self.n_rot_elems + self.n_sta_elems

        self.p = np.hstack([mesh.rot_p, mesh.sta_p])
        self.t = np.hstack([mesh.rot_t, mesh.sta_t + self.n_rot_nodes])
        self.reg = np.concatenate([mesh.rot_reg, mesh.sta_reg])
        self.areas, self.gx, self.gy, self.xc, self.yc = _precompute_tri_geometry(self.p, self.t)

        self.arkkio_elems = np.concatenate([
            mesh.rot_arkkio_elems,
            mesh.sta_arkkio_elems + self.n_rot_elems
        ])
        self.arkkio_r = np.hypot(self.xc[self.arkkio_elems], self.yc[self.arkkio_elems])
        self.arkkio_cos = self.xc[self.arkkio_elems] / self.arkkio_r
        self.arkkio_sin = self.yc[self.arkkio_elems] / self.arkkio_r
        self.arkkio_dg = mesh.r_g3 - mesh.r_g1

        self.mat_code = np.zeros(self.n_tot_elems, dtype=np.int32)
        self.nu_lin = np.full(self.n_tot_elems, NU0, dtype=np.float64)
        self.Brx = np.zeros(self.n_tot_elems, dtype=np.float64)
        self.Bry = np.zeros(self.n_tot_elems, dtype=np.float64)
        self.mag_dir_x = np.zeros(self.n_tot_elems, dtype=np.float64)
        self.mag_dir_y = np.zeros(self.n_tot_elems, dtype=np.float64)

        self.mat_code[(self.reg == REG_ROTOR_STEEL) | (self.reg == REG_STATOR_STEEL)] = 1
        self.mat_code[(self.reg == REG_ROTOR_STEEL_CUT) | (self.reg == REG_STATOR_STEEL_CUT)] = 2

        self.mag_props = MAGNET_DB[design.magnet_grade]
        self.set_magnet_temperature(mag_temp_mode)

        elem_idx_all = np.arange(self.n_tot_elems)
        self.slot_elem_masks = [
            (self.reg == (REG_SLOT_BASE + s)) & (elem_idx_all >= self.n_rot_elems)
            for s in range(self.n_spp)
        ]
        self.all_slot_elems = np.where(
            (self.reg >= REG_SLOT_BASE) & (self.reg < REG_SLOT_BASE + self.n_spp)
            & (elem_idx_all >= self.n_rot_elems)
        )[0]
        self.rot_steel_elems = np.where(
            (self.reg == REG_ROTOR_STEEL) | (self.reg == REG_ROTOR_STEEL_CUT)
        )[0]
        self.sta_steel_elems = np.where(
            (self.reg == REG_STATOR_STEEL) | (self.reg == REG_STATOR_STEEL_CUT)
        )[0]
        self.mag_elems = np.where(self.reg >= REG_MAGNET_BASE)[0]

        # Exact 3D end-winding length (scales inversely with p_pairs -> realistic 4p vs 6p vs 8p penalty!)
        r_slot_mid = 0.5 * (mesh.r_si + (0.100 - design.h_yoke_mm * 1e-3))
        tau_p = math.pi * r_slot_mid / self.p_pairs
        self.l_ew = 1.55 * tau_p + 0.025
        self.slot_area_mean = float(np.mean(mesh.slot_areas))
        Nc = design.N_turns
        self.R_ph = (
            CU_RHO_120C * (self.n_poles * self.q_spp * (Nc ** 2) * (L_STK + self.l_ew))
            / (K_FILL_MAX * self.slot_area_mean)
        )
        self.L_ew = MU0 * self.n_poles * self.q_spp * (Nc ** 2) * self.l_ew * 0.30
        self._dof_cache: Dict[int, Tuple[np.ndarray, np.ndarray, int]] = {}

    def set_magnet_temperature(self, mode: str):
        self.mag_temp_mode = mode
        mp = self.mag_props
        Br = mp["Br_140C"] if mode == "HOT_140C" else (mp["Br_m30C"] if mode == "COLD_M30C" else mp["Br_20C"])
        nu_m = NU0 / mp["mu_r"]
        self.Brx[:] = 0.0
        self.Bry[:] = 0.0
        for m_tag, (mx, my) in self.mesh.rot_mag_dir.items():
            mask = (self.reg == m_tag)
            self.mat_code[mask] = 3
            self.nu_lin[mask] = nu_m
            self.Brx[mask] = Br * mx
            self.Bry[mask] = Br * my
            self.mag_dir_x[mask] = mx
            self.mag_dir_y[mask] = my

    def get_dof_map(self, step_k: int) -> Tuple[np.ndarray, np.ndarray, int]:
        step_k = step_k % (2 * N_THETA_GAP)
        if step_k in self._dof_cache:
            return self._dof_cache[step_k]

        m = self.mesh
        nr = self.n_rot_nodes
        nt = self.n_tot_nodes
        master = np.arange(nt, dtype=np.int64)
        sign = np.ones(nt, dtype=np.float64)
        is_dirichlet = np.zeros(nt, dtype=np.bool_)

        is_dirichlet[m.rot_shaft_nodes] = True
        is_dirichlet[m.sta_outer_nodes + nr] = True

        for n0, np_ in zip(m.sta_per_0 + nr, m.sta_per_p + nr):
            if not is_dirichlet[np_]:
                master[np_] = n0
                sign[np_] = -1.0

        for n0, np_ in zip(m.rot_per_0, m.rot_per_p):
            if not is_dirichlet[np_]:
                master[np_] = n0
                sign[np_] = -1.0

        sta_mid_g = m.sta_gap_mid + nr
        for j in range(N_THETA_GAP + 1):
            r_node = m.rot_gap_mid[j]
            tot = j + step_k
            n_wraps = tot // N_THETA_GAP
            j_s = tot % N_THETA_GAP
            s_node = sta_mid_g[j_s]
            s_sgn = -1.0 if (n_wraps % 2 == 1) else 1.0
            master[r_node] = master[s_node]
            sign[r_node] = s_sgn * sign[s_node]

        for _ in range(3):
            for i in range(nt):
                if is_dirichlet[i]:
                    continue
                p_i = master[i]
                if p_i != i:
                    if is_dirichlet[p_i]:
                        is_dirichlet[i] = True
                    else:
                        sign[i] *= sign[p_i]
                        master[i] = master[p_i]

        dof_map = np.full(nt, -1, dtype=np.int64)
        dof_sign = np.zeros(nt, dtype=np.float64)
        free_masters = np.where((~is_dirichlet) & (master == np.arange(nt)))[0]
        master_to_free = np.full(nt, -1, dtype=np.int64)
        master_to_free[free_masters] = np.arange(len(free_masters), dtype=np.int64)

        for i in range(nt):
            if not is_dirichlet[i]:
                dof_map[i] = master_to_free[master[i]]
                dof_sign[i] = sign[i]

        n_free = len(free_masters)
        self._dof_cache[step_k] = (dof_map, dof_sign, n_free)
        return dof_map, dof_sign, n_free

    def solve_step(
        self, step_k: int, i_abc: Tuple[float, float, float],
        A_init: Optional[np.ndarray] = None, tol: float = 1e-6, max_iter: int = 18
    ) -> dict:
        dof_map, dof_sign, n_free = self.get_dof_map(step_k)
        Jz = np.zeros(self.n_tot_elems, dtype=np.float64)
        Nc = self.design.N_turns
        for s in range(self.n_spp):
            ph = self.slot_phase[s]
            pol = self.slot_pol[s]
            Jz[self.slot_elem_masks[s]] = pol * Nc * i_abc[ph] / self.mesh.slot_areas[s]

        A_free = A_init.copy() if (A_init is not None and len(A_init) == n_free) else np.zeros(n_free, dtype=np.float64)

        converged = False
        res_norm_0 = None
        for it in range(max_iter):
            rows, cols, vals, res, Bx, By = _assemble_nonlinear_maxwell(
                A_free, self.t, self.areas, self.gx, self.gy,
                self.mat_code, self.nu_lin, self.Brx, self.Bry, Jz,
                dof_map, dof_sign, n_free,
            )
            rnorm = float(np.linalg.norm(res))
            if res_norm_0 is None:
                res_norm_0 = max(rnorm, 1e-9)
            if it > 0 and rnorm <= max(tol * res_norm_0, 1e-11):
                converged = True
                break

            K = sp.coo_matrix((vals, (rows, cols)), shape=(n_free, n_free)).tocsc()
            dA = spla.spsolve(K, -res)
            alpha = 1.0
            A_cand = A_free + dA
            if it >= 1:
                step_ratio = float(np.linalg.norm(dA)) / max(float(np.linalg.norm(A_free)), 1e-9)
                if step_ratio > 1.5:
                    for _ in range(4):
                        _, _, _, res_c, _, _ = _assemble_nonlinear_maxwell(
                            A_free + alpha * dA, self.t, self.areas, self.gx, self.gy,
                            self.mat_code, self.nu_lin, self.Brx, self.Bry, Jz,
                            dof_map, dof_sign, n_free,
                        )
                        if float(np.linalg.norm(res_c)) < rnorm * 1.05:
                            break
                        alpha *= 0.5
                    A_cand = A_free + alpha * dA
            A_free = A_cand

        _, _, _, _, Bx, By = _assemble_nonlinear_maxwell(
            A_free, self.t, self.areas, self.gx, self.gy,
            self.mat_code, self.nu_lin, self.Brx, self.Bry, Jz,
            dof_map, dof_sign, n_free,
        )
        A_full = np.where(dof_map >= 0, dof_sign * A_free[np.maximum(dof_map, 0)], 0.0)

        ae = self.arkkio_elems
        Bxa, Bya = Bx[ae], By[ae]
        Br = Bxa * self.arkkio_cos + Bya * self.arkkio_sin
        Bth = -Bxa * self.arkkio_sin + Bya * self.arkkio_cos
        torque_arkkio = float(
            (self.n_poles * L_STK / (MU0 * self.arkkio_dg))
            * np.sum(self.arkkio_r * Br * Bth * self.areas[ae])
        )

        A_elem_mean = (A_full[self.t[0]] + A_full[self.t[1]] + A_full[self.t[2]]) / 3.0
        lam_abc = np.zeros(3, dtype=np.float64)
        for s in range(self.n_spp):
            mask = self.slot_elem_masks[s]
            A_bar = float(np.sum(A_elem_mean[mask] * self.areas[mask]) / self.mesh.slot_areas[s])
            lam_abc[self.slot_phase[s]] += self.n_poles * L_STK * Nc * self.slot_pol[s] * A_bar

        th_mech = step_k * (self.th_pole / N_THETA_GAP)
        gamma_e = GAMMA_E0_RAD + self.p_pairs * th_mech
        c0 = math.cos(gamma_e)
        c1 = math.cos(gamma_e - 2.0 * math.pi / 3.0)
        c2 = math.cos(gamma_e + 2.0 * math.pi / 3.0)
        s0 = -math.sin(gamma_e)
        s1 = -math.sin(gamma_e - 2.0 * math.pi / 3.0)
        s2 = -math.sin(gamma_e + 2.0 * math.pi / 3.0)
        lam_d = (2.0 / 3.0) * (lam_abc[0] * c0 + lam_abc[1] * c1 + lam_abc[2] * c2)
        lam_q = (2.0 / 3.0) * (lam_abc[0] * s0 + lam_abc[1] * s1 + lam_abc[2] * s2)

        min_B_mag = 99.0
        p02_B_mag = 99.0
        if len(self.mag_elems) > 0:
            B_par = Bx[self.mag_elems] * self.mag_dir_x[self.mag_elems] + By[self.mag_elems] * self.mag_dir_y[self.mag_elems]
            min_B_mag = float(np.min(B_par))
            p02_B_mag = float(np.percentile(B_par, 2.0))

        W_mag_full = self.n_poles * L_STK * _compute_coenergy_1pole(
            Bx, By, self.areas, self.mat_code, self.nu_lin, self.Brx, self.Bry
        )
        W_coenergy = float(lam_abc[0] * i_abc[0] + lam_abc[1] * i_abc[1] + lam_abc[2] * i_abc[2] - W_mag_full)

        return {
            "converged": converged,
            "iters": it + 1,
            "A_free": A_free,
            "A_full": A_full,
            "Bx": Bx,
            "By": By,
            "torque": torque_arkkio,
            "lam_abc": lam_abc,
            "lam_d": lam_d,
            "lam_q": lam_q,
            "min_B_mag": min_B_mag,
            "p02_B_mag": p02_B_mag,
            "W_coenergy": W_coenergy,
            "gamma_e": gamma_e,
        }

    def solve_dq_sweep(
        self, Id: float, Iq: float, rpm: float,
        n_steps: int = 6,
        k_h5: float = 0.0, phi_h5: float = 0.0,
        k_h7: float = 0.0, phi_h7: float = 0.0,
        include_end_winding: bool = True,
    ) -> dict:
        assert 24 % n_steps == 0
        stride = 24 // n_steps
        I1_pk = math.hypot(Id, Iq)
        I5_pk = k_h5 * I1_pk
        I7_pk = k_h7 * I1_pk

        torques = np.empty(n_steps, dtype=np.float64)
        lam_ds = np.empty(n_steps, dtype=np.float64)
        lam_qs = np.empty(n_steps, dtype=np.float64)
        Bx_sta = np.empty((n_steps, len(self.sta_steel_elems)), dtype=np.float64)
        By_sta = np.empty((n_steps, len(self.sta_steel_elems)), dtype=np.float64)
        Bx_rot = np.empty((n_steps, len(self.rot_steel_elems)), dtype=np.float64)
        By_rot = np.empty((n_steps, len(self.rot_steel_elems)), dtype=np.float64)
        Bx_slot = np.empty((n_steps, len(self.all_slot_elems)), dtype=np.float64)
        By_slot = np.empty((n_steps, len(self.all_slot_elems)), dtype=np.float64)

        min_B_mag = 99.0
        p02_B_mag = 99.0
        A_prev = None
        all_conv = True

        for idx in range(n_steps):
            step_k = idx * stride
            th_mech = step_k * (self.th_pole / N_THETA_GAP)
            ge = GAMMA_E0_RAD + self.p_pairs * th_mech
            ia = (Id * math.cos(ge) - Iq * math.sin(ge)
                  + I5_pk * math.cos(5.0 * ge + phi_h5) + I7_pk * math.cos(7.0 * ge + phi_h7))
            ib = (Id * math.cos(ge - 2.0 * math.pi / 3.0) - Iq * math.sin(ge - 2.0 * math.pi / 3.0)
                  + I5_pk * math.cos(5.0 * (ge - 2.0 * math.pi / 3.0) + phi_h5)
                  + I7_pk * math.cos(7.0 * (ge - 2.0 * math.pi / 3.0) + phi_h7))
            ic = -(ia + ib)

            res = self.solve_step(step_k, (ia, ib, ic), A_init=A_prev)
            all_conv = all_conv and res["converged"]
            A_prev = res["A_free"]

            torques[idx] = res["torque"]
            lam_ds[idx] = res["lam_d"]
            lam_qs[idx] = res["lam_q"]
            Bx_sta[idx] = res["Bx"][self.sta_steel_elems]
            By_sta[idx] = res["By"][self.sta_steel_elems]
            Bx_rot[idx] = res["Bx"][self.rot_steel_elems]
            By_rot[idx] = res["By"][self.rot_steel_elems]
            Bx_slot[idx] = res["Bx"][self.all_slot_elems]
            By_slot[idx] = res["By"][self.all_slot_elems]
            min_B_mag = min(min_B_mag, res["min_B_mag"])
            p02_B_mag = min(p02_B_mag, res["p02_B_mag"])

        T_avg = float(np.mean(torques))
        T_rip_pct = float((np.max(torques) - np.min(torques)) / max(abs(T_avg), 1.0) * 100.0)

        lam_d_avg = float(np.mean(lam_ds)) + (self.L_ew * Id if include_end_winding else 0.0)
        lam_q_avg = float(np.mean(lam_qs)) + (self.L_ew * Iq if include_end_winding else 0.0)

        omega_m = rpm * 2.0 * math.pi / 60.0
        omega_e = self.p_pairs * omega_m
        f_e = omega_e / (2.0 * math.pi)

        Vd = self.R_ph * Id - omega_e * lam_q_avg
        Vq = self.R_ph * Iq + omega_e * lam_d_avg
        V_pk_fund = math.hypot(Vd, Vq)
        L_q_eff = max(abs(lam_q_avg / Iq) if abs(Iq) > 1.0 else 0.001, 0.0004)
        V_pk_harm = omega_e * L_q_eff * (5.0 * I5_pk + 7.0 * I7_pk) * 0.45
        V_ph_pk = math.hypot(V_pk_fund, V_pk_harm)

        I_rms = math.sqrt(0.5 * (Id * Id + Iq * Iq + I5_pk * I5_pk + I7_pk * I7_pk))
        I_pk_eff = I_rms * math.sqrt(2.0)

        # 1. DC Copper Loss (120 C, incl. 3D end-winding)
        P_cu_dc = 3.0 * (I_rms ** 2) * self.R_ph

        # 2. AC Proximity Copper Loss from exact FEA slot leakage B(t)
        B_slot_sq_mean = np.mean(Bx_slot ** 2 + By_slot ** 2, axis=0)
        B_slot_rip_sq = np.var(Bx_slot, axis=0) + np.var(By_slot, axis=0)
        V_slot_cu = self.n_poles * L_STK * K_FILL_MAX * self.areas[self.all_slot_elems]
        P_cu_ac = float(np.sum(
            V_slot_cu * (CU_WIRE_DIAM ** 2) / (16.0 * CU_RHO_120C)
            * ((omega_e ** 2) * B_slot_sq_mean + ((6.0 * omega_e) ** 2) * B_slot_rip_sq)
        ))

        # 3. Bertotti Core Loss (Stator + Rotor, with 1.75x hysteresis in 0.5 mm cut-edge zone)
        sta_vol = self.n_poles * L_STK * self.areas[self.sta_steel_elems]
        sta_kh = np.where(self.reg[self.sta_steel_elems] == REG_STATOR_STEEL_CUT,
                          KH_VIRGIN * KH_CUT_MULT, KH_VIRGIN)
        B2_sta_mean = np.mean(Bx_sta ** 2 + By_sta ** 2, axis=0)
        B_sta_pk1 = np.sqrt(np.maximum(2.0 * B2_sta_mean, 1e-14))
        B2_sta_harm = np.var(Bx_sta, axis=0) + np.var(By_sta, axis=0)
        B_sta_pkh = np.sqrt(np.maximum(2.0 * B2_sta_harm, 1e-14))

        P_fe_sta = float(np.sum(sta_vol * (
            sta_kh * f_e * (B_sta_pk1 ** 2)
            + KE_STEEL * (f_e ** 2) * (B_sta_pk1 ** 2)
            + KA_STEEL * (f_e ** 1.5) * (B_sta_pk1 ** 1.5)
            + sta_kh * (6.0 * f_e) * (B_sta_pkh ** 2)
            + KE_STEEL * ((6.0 * f_e) ** 2) * (B_sta_pkh ** 2)
            + KA_STEEL * ((6.0 * f_e) ** 1.5) * (B_sta_pkh ** 1.5)
        )))

        rot_vol = self.n_poles * L_STK * self.areas[self.rot_steel_elems]
        rot_kh = np.where(self.reg[self.rot_steel_elems] == REG_ROTOR_STEEL_CUT,
                          KH_VIRGIN * KH_CUT_MULT, KH_VIRGIN)
        fx_r = np.fft.rfft(Bx_rot, axis=0) / n_steps
        fy_r = np.fft.rfft(By_rot, axis=0) / n_steps
        P_fe_rot = 0.0
        for h_idx in range(1, fx_r.shape[0]):
            scale = 2.0 if (h_idx < fx_r.shape[0] - 1 or n_steps % 2 == 1) else 1.0
            Bh = scale * np.sqrt(np.abs(fx_r[h_idx]) ** 2 + np.abs(fy_r[h_idx]) ** 2)
            fh = 6.0 * h_idx * f_e
            P_fe_rot += float(np.sum(rot_vol * (
                rot_kh * fh * (Bh ** 2)
                + KE_STEEL * (fh ** 2) * (Bh ** 2)
                + KA_STEEL * (fh ** 1.5) * (Bh ** 1.5)
            )))

        P_mech_loss = 12.0 * (rpm / 1000.0) + 0.18 * ((rpm / 1000.0) ** 2.6)
        P_shaft = T_avg * omega_m - P_mech_loss
        P_loss_motor = P_cu_dc + P_cu_ac + P_fe_sta + P_fe_rot + P_mech_loss
        P_elec_in = P_shaft + P_loss_motor

        S_va = 1.5 * V_ph_pk * I_pk_eff
        cos_phi = float(np.clip(P_elec_in / max(S_va, 1.0), 0.0, 0.999))
        J_rms_Amm2 = (self.design.N_turns * I_rms) / (K_FILL_MAX * self.slot_area_mean * 1e6)

        return {
            "converged": all_conv,
            "T_avg": T_avg,
            "T_ripple_pct": T_rip_pct,
            "torques": torques,
            "lam_d": lam_d_avg,
            "lam_q": lam_q_avg,
            "Vd": Vd,
            "Vq": Vq,
            "V_ph_pk": V_ph_pk,
            "I_rms": I_rms,
            "I_pk": I_pk_eff,
            "J_rms_Amm2": J_rms_Amm2,
            "S_kVA": S_va / 1000.0,
            "cos_phi": cos_phi,
            "P_shaft_W": P_shaft,
            "P_elec_W": P_elec_in,
            "P_cu_dc_W": P_cu_dc,
            "P_cu_ac_W": P_cu_ac,
            "P_fe_sta_W": P_fe_sta,
            "P_fe_rot_W": P_fe_rot,
            "P_mech_W": P_mech_loss,
            "eta_motor_pct": float(100.0 * P_shaft / max(P_elec_in, 1.0)) if P_shaft > 0 else 0.0,
            "min_B_mag_T": min_B_mag,
            "p02_B_mag_T": p02_B_mag,
        }


# ---------------------------------------------------------------------------
# Exact 2D Plane-Stress Centrifugal Elasticity FEA on Rotor Laminations
# ---------------------------------------------------------------------------
@njit(cache=True)
def _assemble_elasticity_2d(
    p: np.ndarray, t: np.ndarray, areas: np.ndarray,
    gx: np.ndarray, gy: np.ndarray, xc: np.ndarray, yc: np.ndarray,
    E_elem: np.ndarray, nu_elem: np.ndarray, rho_elem: np.ndarray,
    omega_rad_s: float,
):
    E_cnt = t.shape[1]
    N = p.shape[1]
    rows = np.empty(36 * E_cnt, dtype=np.int64)
    cols = np.empty(36 * E_cnt, dtype=np.int64)
    vals = np.empty(36 * E_cnt, dtype=np.float64)
    F = np.zeros(2 * N, dtype=np.float64)
    om2 = omega_rad_s * omega_rad_s
    nnz = 0

    for e in range(E_cnt):
        Em = E_elem[e]
        if Em <= 0.0:
            continue
        nu = nu_elem[e]
        area = areas[e]
        c_mat = Em / (1.0 - nu * nu)
        D11 = c_mat
        D12 = c_mat * nu
        D33 = c_mat * 0.5 * (1.0 - nu)

        fx = rho_elem[e] * om2 * xc[e] * (area / 3.0)
        fy = rho_elem[e] * om2 * yc[e] * (area / 3.0)

        dofs = np.empty(6, dtype=np.int64)
        B0 = np.empty(6, dtype=np.float64)
        B1 = np.empty(6, dtype=np.float64)
        B2 = np.empty(6, dtype=np.float64)
        for k in range(3):
            nk = t[k, e]
            dofs[2 * k] = 2 * nk
            dofs[2 * k + 1] = 2 * nk + 1
            F[2 * nk] += fx
            F[2 * nk + 1] += fy
            B0[2 * k] = gx[e, k]
            B0[2 * k + 1] = 0.0
            B1[2 * k] = 0.0
            B1[2 * k + 1] = gy[e, k]
            B2[2 * k] = gy[e, k]
            B2[2 * k + 1] = gx[e, k]

        for i in range(6):
            db0 = D11 * B0[i] + D12 * B1[i]
            db1 = D12 * B0[i] + D11 * B1[i]
            db2 = D33 * B2[i]
            ri = dofs[i]
            for j in range(6):
                kij = area * (B0[j] * db0 + B1[j] * db1 + B2[j] * db2)
                rows[nnz] = ri
                cols[nnz] = dofs[j]
                vals[nnz] = kij
                nnz += 1

    return rows[:nnz], cols[:nnz], vals[:nnz], F


def solve_centrifugal_stress(design: MotorDesign, mesh: MotorMeshBundle,
                             overspeed_rpm: float = 14400.0) -> dict:
    """Solve exact 2D plane-stress centrifugal elasticity FEA on the 1-pole rotor with
    exact rotational periodicity R(theta_pole) across 4p, 6p, or 8p topologies."""
    rot_p = mesh.rot_p
    rot_t = mesh.rot_t
    rot_reg = mesh.rot_reg

    is_steel = (rot_reg == REG_ROTOR_STEEL) | (rot_reg == REG_ROTOR_STEEL_CUT)
    is_mag = (rot_reg >= REG_MAGNET_BASE)
    is_bar_air = (rot_reg == REG_AIR_BARRIER)
    solid_mask = is_steel | is_mag | is_bar_air
    solid_t_raw = rot_t[:, solid_mask]
    solid_steel_flag = is_steel[solid_mask]
    solid_mag_flag = is_mag[solid_mask]

    used_nodes = np.unique(solid_t_raw)
    old_to_new = np.full(rot_p.shape[1], -1, dtype=np.int64)
    old_to_new[used_nodes] = np.arange(len(used_nodes), dtype=np.int64)
    sp_pts = rot_p[:, used_nodes]
    sp_t = old_to_new[solid_t_raw]
    N = sp_pts.shape[1]

    areas, gx, gy, xc, yc = _precompute_tri_geometry(sp_pts, sp_t)
    mag_rho = MAGNET_DB[design.magnet_grade]["rho"]
    E_elem = np.where(solid_steel_flag, STEEL_E_MOD, np.where(solid_mag_flag, 1.5e9, 15.0e6))
    nu_elem = np.where(solid_steel_flag, STEEL_POISSON, 0.25)
    rho_elem = np.where(solid_steel_flag, STEEL_RHO, np.where(solid_mag_flag, mag_rho, 0.0))

    omega = overspeed_rpm * 2.0 * math.pi / 60.0
    rows, cols, vals, F_full = _assemble_elasticity_2d(
        sp_pts, sp_t, areas, gx, gy, xc, yc, E_elem, nu_elem, rho_elem, omega
    )
    K_full = sp.coo_matrix((vals, (rows, cols)), shape=(2 * N, 2 * N)).tocsc()

    shaft_new = old_to_new[mesh.rot_shaft_nodes]
    shaft_new = shaft_new[shaft_new >= 0]
    fixed_nodes = set(int(x) for x in shaft_new)

    slave_to_master = {}
    for n0_old, np_old in zip(mesh.rot_per_0, mesh.rot_per_p):
        n0 = int(old_to_new[n0_old])
        np_ = int(old_to_new[np_old])
        if n0 >= 0 and np_ >= 0 and n0 not in fixed_nodes and np_ not in fixed_nodes:
            slave_to_master[np_] = n0

    free_nodes = [n for n in range(N) if n not in fixed_nodes and n not in slave_to_master]
    node_to_free = {n: idx for idx, n in enumerate(free_nodes)}
    n_free_dofs = 2 * len(free_nodes)

    cp = math.cos(mesh.theta_pole)
    sp_ = math.sin(mesh.theta_pole)
    p_rows, p_cols, p_vals = [], [], []
    for n in free_nodes:
        fi = node_to_free[n]
        p_rows.extend([2 * n, 2 * n + 1])
        p_cols.extend([2 * fi, 2 * fi + 1])
        p_vals.extend([1.0, 1.0])
    for np_, n0 in slave_to_master.items():
        fi = node_to_free[n0]
        p_rows.extend([2 * np_, 2 * np_, 2 * np_ + 1, 2 * np_ + 1])
        p_cols.extend([2 * fi, 2 * fi + 1, 2 * fi, 2 * fi + 1])
        p_vals.extend([cp, -sp_, sp_, cp])

    P = sp.coo_matrix((p_vals, (p_rows, p_cols)), shape=(2 * N, n_free_dofs)).tocsc()
    K_red = (P.T @ K_full @ P).tocsc()
    F_red = P.T @ F_full
    u_free = spla.spsolve(K_red, F_red)
    u_full = P @ u_free
    ux = u_full[0::2]
    uy = u_full[1::2]

    steel_idx = np.where(solid_steel_flag)[0]
    exx = ux[sp_t[0, steel_idx]] * gx[steel_idx, 0] + ux[sp_t[1, steel_idx]] * gx[steel_idx, 1] + ux[sp_t[2, steel_idx]] * gx[steel_idx, 2]
    eyy = uy[sp_t[0, steel_idx]] * gy[steel_idx, 0] + uy[sp_t[1, steel_idx]] * gy[steel_idx, 1] + uy[sp_t[2, steel_idx]] * gy[steel_idx, 2]
    gxy = (ux[sp_t[0, steel_idx]] * gy[steel_idx, 0] + ux[sp_t[1, steel_idx]] * gy[steel_idx, 1] + ux[sp_t[2, steel_idx]] * gy[steel_idx, 2]
           + uy[sp_t[0, steel_idx]] * gx[steel_idx, 0] + uy[sp_t[1, steel_idx]] * gx[steel_idx, 1] + uy[sp_t[2, steel_idx]] * gx[steel_idx, 2])

    c_st = STEEL_E_MOD / (1.0 - STEEL_POISSON ** 2)
    sxx = c_st * (exx + STEEL_POISSON * eyy)
    syy = c_st * (STEEL_POISSON * exx + eyy)
    txy = c_st * 0.5 * (1.0 - STEEL_POISSON) * gxy
    vm_elem = np.sqrt(np.maximum(sxx * sxx - sxx * syy + syy * syy + 3.0 * txy * txy, 0.0))

    vm_nodal_num = np.zeros(N, dtype=np.float64)
    vm_nodal_den = np.zeros(N, dtype=np.float64)
    for loc_i, e in enumerate(steel_idx):
        a = areas[e]
        v = vm_elem[loc_i]
        for k in range(3):
            nk = sp_t[k, e]
            vm_nodal_num[nk] += a * v
            vm_nodal_den[nk] += a
    valid_n = vm_nodal_den > 0
    vm_nodal = vm_nodal_num[valid_n] / vm_nodal_den[valid_n]

    vm_max_MPa = float(np.max(vm_elem) / 1e6)
    vm_nodal_max_MPa = float(np.max(vm_nodal) / 1e6)
    vm_p99_MPa = float(np.percentile(vm_elem, 99.5) / 1e6)
    max_disp_um = float(np.max(np.hypot(ux, uy)) * 1e6)
    sf_nodal = float((STEEL_SIGMA_Y / 1e6) / max(vm_nodal_max_MPa, 1e-6))

    return {
        "vm_nodal_max_MPa": vm_nodal_max_MPa,
        "vm_elem_max_MPa": vm_max_MPa,
        "vm_p995_MPa": vm_p99_MPa,
        "safety_factor": sf_nodal,
        "max_disp_um": max_disp_um,
        "vm_elem_steel_MPa": vm_elem / 1e6,
        "steel_t_local": sp_t[:, steel_idx],
        "steel_p_local": sp_pts,
    }


if __name__ == "__main__":
    import time
    from motor_geometry import (
        build_motor_mesh, preset_ndfeb_vipm, preset_pure_synrm_4p,
        preset_ferrite_pmasynrm_6p, preset_ferrite_pmasynrm_4p,
    )

    print("=== MULTI-TOPOLOGY (8p / 6p / 4p) SALIENCY, TORQUE CROSS-CHECK & STRESS BENCHMARK ===")
    for fn in (preset_ndfeb_vipm, preset_pure_synrm_4p, preset_ferrite_pmasynrm_6p, preset_ferrite_pmasynrm_4p):
        d = fn()
        t0 = time.perf_counter()
        m = build_motor_mesh(d, mesh_scale=1.15)
        t_mesh = (time.perf_counter() - t0) * 1000.0
        s = MotorFEASolver(d, m, "HOT_140C")
        r0 = s.solve_step(0, (0, 0, 0))
        lam_m = r0["lam_d"]
        # Pick gamma = 55 deg at Ipk = 180 A
        Ipk = 180.0
        Id = -Ipk * math.sin(math.radians(55.0))
        Iq =  Ipk * math.cos(math.radians(55.0))
        t1 = time.perf_counter()
        sw = s.solve_dq_sweep(Id, Iq, rpm=4000.0, n_steps=6)
        t_sw = (time.perf_counter() - t1) * 1000.0
        Ld = (sw["lam_d"] - lam_m) / Id
        Lq = sw["lam_q"] / Iq
        xi = Lq / max(Ld, 1e-9)
        st = solve_centrifugal_stress(d, m, 14400.0)
        # Virtual work cross-check at step 2
        th_m = 2 * (m.theta_pole / N_THETA_GAP)
        ge = GAMMA_E0_RAD + m.p_pairs * th_m
        ia = Id * math.cos(ge) - Iq * math.sin(ge)
        ib = Id * math.cos(ge - 2 * math.pi / 3) - Iq * math.sin(ge - 2 * math.pi / 3)
        ic = -(ia + ib)
        rp = s.solve_step(3, (ia, ib, ic))
        rc = s.solve_step(2, (ia, ib, ic), A_init=rp["A_free"])
        rm = s.solve_step(1, (ia, ib, ic), A_init=rc["A_free"])
        T_vw = (rp["W_coenergy"] - rm["W_coenergy"]) / (2.0 * m.theta_pole / N_THETA_GAP)
        vw_err = abs(rc["torque"] - T_vw) / max(abs(T_vw), 1.0) * 100.0

        print(f"\n[{d.name}] ({d.n_poles}p/{d.q_slots}s, mesh={t_mesh:.0f}ms, 6-step={t_sw:.0f}ms, tri={s.n_tot_elems})")
        print(f"  No-Load: lam_m={lam_m*1000:5.2f} mWb, lam_q={r0['lam_q']*1000:+5.2f} mWb | R_ph={s.R_ph*1e3:.2f} mOhm, L_ew={s.L_ew*1e6:.1f} uH")
        print(f"  @4000rpm g=55deg: T={sw['T_avg']:5.1f} Nm (VW err={vw_err:.2f}%), Rip={sw['T_ripple_pct']:4.1f}%, cos(phi)={sw['cos_phi']:.3f}, xi={xi:.2f}, eta={sw['eta_motor_pct']:.2f}%")
        print(f"  Stress @14.4krpm: vm_nodal={st['vm_nodal_max_MPa']:.1f} MPa, SF={st['safety_factor']:.2f}")

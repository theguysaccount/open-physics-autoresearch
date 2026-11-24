"""
motor_geometry.py — Watertight Parametric 2D Gmsh Geometry & Sliding-Band Mesh
Generator for Multi-Topology Traction Motors (4p/36s, 4p/48s, 6p/36s, 6p/54s, 8p/48s).

1-Pole Anti-Periodic Domain (theta in [0, pi / p_pairs]):
  - Exact node-matching sliding ring at r_mid = (R_ro + R_si)/2 with N_theta = 72 uniform
    segments per pole (divisible by 6, 9, and 12 slots/pole, and 24 steps = 60 deg elec!).
  - Pure structured 2-layer uniform air-gap annulus [r_g1, r_g3] generated analytically
    for zero-noise Arkkio Maxwell Stress torque integration.
  - Exact rotational periodicity on theta=0 and theta=theta_pole via Gmsh setPeriodic.
  - Supports 1 to 4 unconstrained asymmetric flux barriers, central + flux-concentrating
    wing magnets, structural center ribs, outer bridges >= 0.80 mm, and surface notches.
  - Automatic 0.5 mm cut-edge lamination tagging on all steel elements bordering air/slots/barriers.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Tuple
import math
import numpy as np
import gmsh


# Locked Motor Envelope Constants (80 kW peak / 40 kW cont EV Traction Class, 400 V DC Bus)
R_SO = 0.100                 # 100.0 mm outer stator radius (200 mm OD)
R_SH = 0.026                 # 26.0 mm shaft radius (52 mm shaft OD)
L_STK = 0.140                # 140.0 mm active stack length (4.398 L active vol)
K_FILL_MAX = 0.45            # Max realistic copper slot fill factor
CUT_EDGE_DEPTH = 0.0005      # 0.5 mm stamped cut-edge degraded steel zone
N_THETA_GAP = 72             # 72 uniform segments per pole (divisible by 6, 9, 12, and 3)

# Allowed Discrete Pole-Slot Topologies: (p_pairs, q_spp) -> (2p poles, Qs = 6*p*q slots)
VALID_TOPOLOGIES = [
    (4, 2),  # 8p / 48s (6 slots/pole, 45 deg pole)
    (3, 2),  # 6p / 36s (6 slots/pole, 60 deg pole)
    (3, 3),  # 6p / 54s (9 slots/pole, 60 deg pole)
    (2, 3),  # 4p / 36s (9 slots/pole, 90 deg pole)
    (2, 4),  # 4p / 48s (12 slots/pole, 90 deg pole)
]

# Physical Region Integer Tags
REG_AIR_GAP_ROT = 1
REG_AIR_GAP_STA = 2
REG_AIR_BARRIER = 3
REG_ROTOR_STEEL = 4
REG_ROTOR_STEEL_CUT = 5
REG_STATOR_STEEL = 6
REG_STATOR_STEEL_CUT = 7
REG_MAGNET_BASE = 100
REG_SLOT_BASE = 10           # Stator slots 0..(3*q - 1) -> tags 10..21


@dataclass
class BarrierSpec:
    """One flux barrier in pole-local coordinates (ordered inner=0 -> outer=n-1).
    Angles beta_L_frac and beta_R_frac are expressed as fractions of the pole half-angle
    theta_center (0.15 .. 0.92), so geometries transfer seamlessly across 4p, 6p, and 8p!"""
    r_frac: float            # Relative radial center within available rotor space [0.08..0.90]
    thick_mm: float          # Barrier thickness (mm) [1.8..9.0]
    half_width_frac: float   # Central flat half-width relative to sector width at u_b [0.06..0.72]
    beta_L_frac: float       # Left wing tip angle / theta_center [0.15..0.92]
    beta_R_frac: float       # Right wing tip angle / theta_center [0.15..0.92]
    bridge_mm: float         # Outer tangential steel bridge thickness (mm, >= 0.80)
    rib_mm: float            # Center radial structural steel rib width (mm, 0.0 or >= 0.80)
    mag_center_fill: float   # Magnet fill fraction in central flat bar [0.0..0.92]
    mag_wing_fill: float     # Magnet fill fraction in left/right wings [0.0..0.88]


@dataclass
class MotorDesign:
    """Complete parametric 4-layer EV/industrial traction motor + harmonic control genome."""
    name: str
    magnet_grade: str = "FERRITE_FB12B"   # "NONE", "FERRITE_FB12B", "NDFEB_N42SH"
    # Layer 1: Discrete Pole-Slot Topology
    p_pairs: int = 3                      # 2 (4-pole), 3 (6-pole), 4 (8-pole)
    q_spp: int = 2                        # slots/pole/phase -> Qs = 6 * p_pairs * q_spp
    # Layer 3: Stator & Air-Gap Micro-Geometry
    R_ro_mm: float = 62.0                 # Rotor outer radius [52.0, 70.0] mm
    gap_mm: float = 0.65                  # Mechanical air gap [0.65, 1.20] mm
    k_tooth: float = 0.48                 # Tooth width fraction of slot pitch [0.36, 0.64]
    h_yoke_mm: float = 17.0               # Stator back-iron thickness [10.0, 29.0] mm
    b_so_mm: float = 2.4                  # Stator slot opening width [1.6, 5.5] mm
    h_tt_mm: float = 1.1                  # Tooth-tip lip height [0.8, 2.2] mm
    N_turns: float = 3.5                  # Effective series conductors per slot
    # Rotor Barriers (1 to 4 barriers, ordered inner -> outer)
    barriers: List[BarrierSpec] = field(default_factory=list)
    # Optional Rotor Surface Flux-Modulation Notches
    notch_depth_mm: float = 0.0           # [0.0, 1.5] mm
    notch_angle_frac: float = 0.84        # Angle from pole center / theta_center [0.55, 0.93]
    notch_width_deg: float = 1.6          # Angular width [0.8, 3.5] deg
    # Layer 4: Co-Evolved Harmonic Current Injection (5th & 7th as fraction of fundamental peak I1)
    k_h5: float = 0.0                     # [0.0, 0.12]
    phi_h5_rad: float = 0.0               # [-pi, pi]
    k_h7: float = 0.0                     # [0.0, 0.08]
    phi_h7_rad: float = 0.0               # [-pi, pi]

    @property
    def n_poles(self) -> int:
        return 2 * self.p_pairs

    @property
    def q_slots(self) -> int:
        return 6 * self.p_pairs * self.q_spp

    @property
    def slots_per_pole(self) -> int:
        return 3 * self.q_spp

    @property
    def theta_pole(self) -> float:
        return math.pi / self.p_pairs

    @property
    def theta_center(self) -> float:
        return 0.5 * math.pi / self.p_pairs


@dataclass
class MotorMeshBundle:
    """Combined stator + rotor + sliding air-gap finite-element mesh bundle."""
    p_pairs: int
    q_spp: int
    slots_per_pole: int
    theta_pole: float
    theta_center: float
    # Rotor mesh
    rot_p: np.ndarray
    rot_t: np.ndarray
    rot_reg: np.ndarray
    rot_mag_dir: Dict[int, Tuple[float, float]]
    rot_per_0: np.ndarray
    rot_per_p: np.ndarray
    rot_gap_mid: np.ndarray
    rot_shaft_nodes: np.ndarray
    rot_arkkio_elems: np.ndarray
    # Stator mesh
    sta_p: np.ndarray
    sta_t: np.ndarray
    sta_reg: np.ndarray
    sta_per_0: np.ndarray
    sta_per_p: np.ndarray
    sta_gap_mid: np.ndarray
    sta_outer_nodes: np.ndarray
    sta_arkkio_elems: np.ndarray
    # Derived geometric scalars
    r_ro: float
    r_si: float
    r_g1: float
    r_mid: float
    r_g3: float
    slot_areas: np.ndarray
    area_rot_steel: float
    area_sta_steel: float
    area_magnets: float
    area_copper_1pole: float
    min_bridge_mm: float
    min_rib_mm: float


def _uv_to_xy(u: float, v: float, th_c: float) -> Tuple[float, float]:
    c = math.cos(th_c)
    s = math.sin(th_c)
    return (u * c - v * s, u * s + v * c)


def _rect_poly_xy(p1: Tuple[float, float], p2: Tuple[float, float],
                  thickness: float) -> List[Tuple[float, float]]:
    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]
    L = math.hypot(dx, dy)
    if L < 1e-7:
        return []
    nx = -dy / L * (0.5 * thickness)
    ny = dx / L * (0.5 * thickness)
    return [
        (p1[0] - nx, p1[1] - ny),
        (p2[0] - nx, p2[1] - ny),
        (p2[0] + nx, p2[1] + ny),
        (p1[0] + nx, p1[1] + ny),
    ]


def _make_occ_polygon(pts: List[Tuple[float, float]]) -> int:
    p_tags = [gmsh.model.occ.addPoint(x, y, 0.0) for (x, y) in pts]
    l_tags = [gmsh.model.occ.addLine(p_tags[i], p_tags[(i + 1) % len(p_tags)])
              for i in range(len(p_tags))]
    loop = gmsh.model.occ.addCurveLoop(l_tags)
    return gmsh.model.occ.addPlaneSurface([loop])


def _make_occ_sector(r_in: float, r_out: float, th0: float, th1: float) -> int:
    c = gmsh.model.occ.addPoint(0.0, 0.0, 0.0)
    p1 = gmsh.model.occ.addPoint(r_in * math.cos(th0), r_in * math.sin(th0), 0.0)
    p2 = gmsh.model.occ.addPoint(r_out * math.cos(th0), r_out * math.sin(th0), 0.0)
    p3 = gmsh.model.occ.addPoint(r_out * math.cos(th1), r_out * math.sin(th1), 0.0)
    p4 = gmsh.model.occ.addPoint(r_in * math.cos(th1), r_in * math.sin(th1), 0.0)
    l1 = gmsh.model.occ.addLine(p1, p2)
    a_out = gmsh.model.occ.addCircleArc(p2, c, p3)
    l2 = gmsh.model.occ.addLine(p3, p4)
    a_in = gmsh.model.occ.addCircleArc(p4, c, p1)
    loop = gmsh.model.occ.addCurveLoop([l1, a_out, l2, a_in])
    return gmsh.model.occ.addPlaneSurface([loop])


def _rot_matrix(angle: float) -> List[float]:
    c = math.cos(angle)
    s = math.sin(angle)
    return [
        c, -s, 0.0, 0.0,
        s,  c, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0,
    ]


def _extract_2d_mesh_from_gmsh() -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    node_coords = node_coords.reshape(-1, 3)
    tag_to_idx = {int(tag): idx for idx, tag in enumerate(node_tags)}

    tri_list = []
    reg_list = []
    for dim, surf_tag in gmsh.model.getEntities(2):
        etypes, _, enodes = gmsh.model.mesh.getElements(2, surf_tag)
        for et, en in zip(etypes, enodes):
            if et == 2:
                tris = en.reshape(-1, 3)
                mapped = np.vectorize(tag_to_idx.__getitem__)(tris)
                tri_list.append(mapped)
                reg_list.append(np.full(len(mapped), surf_tag, dtype=np.int32))

    if not tri_list:
        raise RuntimeError("Gmsh generated zero 2D triangles")
    t = np.vstack(tri_list).T.astype(np.int64)
    regs = np.concatenate(reg_list)
    p = node_coords[:, :2].T.copy()

    used = np.unique(t)
    remap = np.full(p.shape[1], -1, dtype=np.int64)
    remap[used] = np.arange(len(used), dtype=np.int64)
    p = p[:, used]
    t = remap[t]

    x0, y0 = p[0, t[0]], p[1, t[0]]
    x1, y1 = p[0, t[1]], p[1, t[1]]
    x2, y2 = p[0, t[2]], p[1, t[2]]
    area2 = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
    flip = area2 < 0
    if np.any(flip):
        t[1, flip], t[2, flip] = t[2, flip].copy(), t[1, flip].copy()
    return p, t, regs


def _apply_cut_edge_tagging(p: np.ndarray, t: np.ndarray, regs: np.ndarray,
                            steel_tag: int, cut_tag: int, depth: float = CUT_EDGE_DEPTH):
    is_steel = (regs == steel_tag)
    if not np.any(is_steel):
        return regs
    e0 = np.sort(t[[0, 1], :], axis=0)
    e1 = np.sort(t[[1, 2], :], axis=0)
    e2 = np.sort(t[[2, 0], :], axis=0)
    edges = np.hstack([e0, e1, e2])
    elem_ids = np.tile(np.arange(t.shape[1]), 3)
    steel_flag = is_steel[elem_ids]

    n_nodes = p.shape[1]
    keys = edges[0] * n_nodes + edges[1]
    order = np.argsort(keys)
    sk = keys[order]
    sf = steel_flag[order]

    same = sk[:-1] == sk[1:]
    diff_mat = sf[:-1] != sf[1:]
    iface_mask = same & diff_mat
    iface_edges = edges[:, order[:-1][iface_mask]]
    if iface_edges.shape[1] == 0:
        return regs

    emid = 0.5 * (p[:, iface_edges[0]] + p[:, iface_edges[1]])
    steel_idx = np.where(is_steel)[0]
    scent = (p[:, t[0, steel_idx]] + p[:, t[1, steel_idx]] + p[:, t[2, steel_idx]]) / 3.0

    from scipy.spatial import cKDTree
    tree = cKDTree(emid.T)
    dists, _ = tree.query(scent.T, k=1)
    cut_elems = steel_idx[dists <= depth]
    out = regs.copy()
    out[cut_elems] = cut_tag
    return out


def _append_structured_ring(
    p: np.ndarray, t: np.ndarray, regs: np.ndarray,
    ring_inner_nodes: np.ndarray, r_new: float, th_pole: float, reg_tag: int, outward: bool
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    n_pts = len(ring_inner_nodes)
    assert n_pts == N_THETA_GAP + 1
    thetas = np.linspace(0.0, th_pole, n_pts)
    new_pts = np.vstack([r_new * np.cos(thetas), r_new * np.sin(thetas)])
    n_old = p.shape[1]
    new_ids = np.arange(n_old, n_old + n_pts, dtype=np.int64)
    p_all = np.hstack([p, new_pts])

    tris = []
    for j in range(N_THETA_GAP):
        a0, a1 = ring_inner_nodes[j], ring_inner_nodes[j + 1]
        b0, b1 = new_ids[j], new_ids[j + 1]
        if outward:
            tris.append((a0, b0, b1))
            tris.append((a0, b1, a1))
        else:
            tris.append((b0, a0, a1))
            tris.append((b0, a1, b1))
    t_new = np.array(tris, dtype=np.int64).T
    e_old = t.shape[1]
    t_all = np.hstack([t, t_new])
    regs_all = np.concatenate([regs, np.full(t_new.shape[1], reg_tag, dtype=np.int32)])
    new_elem_idx = np.arange(e_old, e_old + t_new.shape[1], dtype=np.int64)
    return p_all, t_all, regs_all, new_ids, new_elem_idx


def _find_arc_nodes(p: np.ndarray, r_target: float, tol: float = 1.5e-5) -> np.ndarray:
    r = np.hypot(p[0], p[1])
    idx = np.where(np.abs(r - r_target) < tol)[0]
    th = np.arctan2(p[1, idx], p[0, idx])
    return idx[np.argsort(th)]


def _find_radial_nodes(p: np.ndarray, th_target: float, tol: float = 1.5e-5) -> np.ndarray:
    r = np.hypot(p[0], p[1])
    th = np.arctan2(p[1], p[0])
    dist = np.abs(p[0] * (-math.sin(th_target)) + p[1] * math.cos(th_target))
    idx = np.where((dist < tol) & (r > 1e-4) & (np.abs(th - th_target) < 0.05))[0]
    return idx[np.argsort(r[idx])]


def _verify_periodic_pair(p: np.ndarray, n0: np.ndarray, np_: np.ndarray, name: str):
    if len(n0) != len(np_):
        raise RuntimeError(f"{name} periodic node count mismatch: {len(n0)} vs {len(np_)}")
    r0 = np.hypot(p[0, n0], p[1, n0])
    rp = np.hypot(p[0, np_], p[1, np_])
    if np.max(np.abs(r0 - rp)) > 1e-4:
        raise RuntimeError(f"{name} periodic radial mismatch: max dr = {np.max(np.abs(r0 - rp)):.6f}")


def validate_design_geometry(design: MotorDesign) -> Tuple[bool, str, float, float]:
    """Fast pre-mesh geometric feasibility check enforcing strict manufacturability."""
    if (design.p_pairs, design.q_spp) not in VALID_TOPOLOGIES:
        return False, f"Invalid topology ({design.p_pairs}, {design.q_spp})", 0.0, 0.0
    if not (50.0 <= design.R_ro_mm <= 70.5):
        return False, "R_ro out of [50, 70.5] mm", 0.0, 0.0
    if not (0.65 <= design.gap_mm <= 1.25):
        return False, "gap out of [0.65, 1.25] mm", 0.0, 0.0
    if not (0.35 <= design.k_tooth <= 0.65):
        return False, "k_tooth out of [0.35, 0.65]", 0.0, 0.0
    if not (10.0 <= design.h_yoke_mm <= 30.0):
        return False, "h_yoke out of [10, 30] mm", 0.0, 0.0
    r_si = design.R_ro_mm + design.gap_mm
    h_slot = R_SO * 1000.0 - r_si - design.h_yoke_mm
    if h_slot < 9.5:
        return False, f"slot depth too shallow ({h_slot:.1f} mm)", 0.0, 0.0
    slot_pitch_mm = 2.0 * math.pi * r_si / design.q_slots
    tooth_w = slot_pitch_mm * design.k_tooth
    if design.b_so_mm < 1.6 or design.b_so_mm > (slot_pitch_mm - tooth_w) * 0.92:
        return False, "slot opening invalid", 0.0, 0.0
    if not (1 <= len(design.barriers) <= 4):
        return False, "must have 1..4 barriers", 0.0, 0.0

    min_br = 99.0
    min_rib = 99.0
    for b in design.barriers:
        if b.bridge_mm < 0.80:
            return False, f"bridge {b.bridge_mm:.2f} < 0.80 mm", 0.0, 0.0
        if b.rib_mm > 0.0 and b.rib_mm < 0.80:
            return False, f"rib {b.rib_mm:.2f} < 0.80 mm", 0.0, 0.0
        if not (1.80 <= b.thick_mm <= 9.50):
            return False, f"barrier thickness {b.thick_mm:.2f} out of [1.8, 9.5] mm", 0.0, 0.0
        if not (0.14 <= b.beta_L_frac <= 0.93 and 0.14 <= b.beta_R_frac <= 0.93):
            return False, "barrier wing angle fraction out of [0.14, 0.93]", 0.0, 0.0
        min_br = min(min_br, b.bridge_mm)
        if b.rib_mm > 0.0:
            min_rib = min(min_rib, b.rib_mm)
    if min_rib > 50.0:
        min_rib = 0.0

    # Strict radial and angular separation between consecutive barriers (guarantees >= 1.6 mm steel channels!)
    r_avail_mm = (design.R_ro_mm - 2.6) - (R_SH * 1000.0 + 3.5)
    for i in range(len(design.barriers) - 1):
        b_in = design.barriers[i]
        b_out = design.barriers[i + 1]
        dr_centers = (b_out.r_frac - b_in.r_frac) * r_avail_mm
        steel_gap = dr_centers - 0.5 * (b_in.thick_mm + b_out.thick_mm)
        if steel_gap < 1.50:
            return False, f"radial steel channel {steel_gap:.2f} < 1.50 mm", 0.0, 0.0
        if b_in.beta_L_frac - b_out.beta_L_frac < 0.11 or b_in.beta_R_frac - b_out.beta_R_frac < 0.11:
            return False, "barrier wings too close angularly", 0.0, 0.0

    return True, "ok", min_br, min_rib


def build_motor_mesh(design: MotorDesign, mesh_scale: float = 1.0) -> MotorMeshBundle:
    """Build complete watertight 1-pole Stator + Rotor + Sliding-Ring mesh via Gmsh OCC."""
    ok, msg, min_br, min_rib = validate_design_geometry(design)
    if not ok:
        raise ValueError(f"Invalid motor geometry: {msg}")

    th_pole = design.theta_pole
    th_center = design.theta_center
    n_spp = design.slots_per_pole
    r_ro = design.R_ro_mm * 1e-3
    gap = design.gap_mm * 1e-3
    r_si = r_ro + gap
    r_g1 = r_ro + 0.25 * gap
    r_mid = r_ro + 0.50 * gap
    r_g3 = r_ro + 0.75 * gap

    if not gmsh.isInitialized():
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.Verbosity", 0)

    def on_ray(xyz, th):
        return abs(xyz[0] * (-math.sin(th)) + xyz[1] * math.cos(th)) < 2.5e-5

    try:
        # =====================================================================
        # PART 1: BUILD STATOR MESH (r_g3 to R_SO) + APPEND [r_mid, r_g3] RING
        # =====================================================================
        gmsh.model.add("stator")
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.00025 * mesh_scale)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.0050 * mesh_scale)

        sta_air_outer = _make_occ_sector(r_g3, r_si, 0.0, th_pole)
        sta_iron_sector = _make_occ_sector(r_si, R_SO, 0.0, th_pole)

        dth_s = th_pole / n_spp
        h_tt = design.h_tt_mm * 1e-3
        h_wedge = 0.0008
        r_slot_bot = R_SO - design.h_yoke_mm * 1e-3
        r_tt_top = r_si + h_tt
        r_wedge_top = r_tt_top + h_wedge

        tooth_w = (2.0 * math.pi * r_si / design.q_slots) * design.k_tooth
        half_so = 0.5 * design.b_so_mm * 1e-3

        slot_cu_tags_raw = []
        slot_open_tags_raw = []
        for s_idx in range(n_spp):
            th_c = (s_idx + 0.5) * dth_s
            w_inner = max(0.0014, 2.0 * r_wedge_top * math.sin(0.5 * dth_s) - tooth_w)
            w_outer = max(w_inner, 2.0 * r_slot_bot * math.sin(0.5 * dth_s) - tooth_w)

            def to_xy(rr, vv):
                c, sn = math.cos(th_c), math.sin(th_c)
                return (rr * c - vv * sn, rr * sn + vv * c)

            r_open_in = r_si - 0.04 * gap
            pts_open = [
                to_xy(r_open_in, -half_so),
                to_xy(r_tt_top, -half_so),
                to_xy(r_wedge_top, -0.5 * w_inner),
                to_xy(r_wedge_top, +0.5 * w_inner),
                to_xy(r_tt_top, +half_so),
                to_xy(r_open_in, +half_so),
            ]
            pts_cu = [
                to_xy(r_wedge_top, -0.5 * w_inner),
                to_xy(r_slot_bot, -0.5 * w_outer),
                to_xy(r_slot_bot, +0.5 * w_outer),
                to_xy(r_wedge_top, +0.5 * w_inner),
            ]
            slot_open_tags_raw.append(_make_occ_polygon(pts_open))
            slot_cu_tags_raw.append(_make_occ_polygon(pts_cu))

        tools = [(2, sta_air_outer)] + [(2, s) for s in slot_open_tags_raw] + [(2, s) for s in slot_cu_tags_raw]
        _, out_map = gmsh.model.occ.fragment([(2, sta_iron_sector)], tools)
        gmsh.model.occ.synchronize()

        cu_surf_to_idx = {}
        for s_idx in range(n_spp):
            for (_, st) in out_map[2 + n_spp + s_idx]:
                cu_surf_to_idx[st] = s_idx
        air_surfs = set(st for (_, st) in out_map[1])
        for s_idx in range(n_spp):
            for (_, st) in out_map[2 + s_idx]:
                if st not in cu_surf_to_idx:
                    air_surfs.add(st)

        c0_list, cp_list = [], []
        for dim, ctag in gmsh.model.getEntities(1):
            bnds = gmsh.model.getBoundary([(1, ctag)], oriented=False)
            if len(bnds) == 2:
                pA = gmsh.model.getValue(0, bnds[0][1], [])
                pB = gmsh.model.getValue(0, bnds[1][1], [])
                rA, rB = math.hypot(pA[0], pA[1]), math.hypot(pB[0], pB[1])
                thA, thB = math.atan2(pA[1], pA[0]), math.atan2(pB[1], pB[0])
                if abs(rA - r_g3) < 2.5e-5 and abs(rB - r_g3) < 2.5e-5 and abs(thB - thA) > 0.5 * th_pole:
                    gmsh.model.mesh.setTransfiniteCurve(ctag, N_THETA_GAP + 1)
                elif on_ray(pA, 0.0) and on_ray(pB, 0.0) and abs(rB - rA) > 1e-5:
                    c0_list.append((0.5 * (rA + rB), ctag))
                elif on_ray(pA, th_pole) and on_ray(pB, th_pole) and abs(rB - rA) > 1e-5:
                    cp_list.append((0.5 * (rA + rB), ctag))

        c0_list.sort()
        cp_list.sort()
        rot_m = _rot_matrix(th_pole)
        for (_, t0), (_, tp) in zip(c0_list, cp_list):
            gmsh.model.mesh.setPeriodic(1, [tp], [t0], rot_m)

        for dim, ptag in gmsh.model.getEntities(0):
            xyz = gmsh.model.getValue(0, ptag, [])
            rr = math.hypot(xyz[0], xyz[1])
            if rr < r_wedge_top + 0.002:
                lc = 0.0012 * mesh_scale
            elif rr < r_slot_bot + 0.002:
                lc = 0.0029 * mesh_scale
            else:
                lc = 0.0044 * mesh_scale
            gmsh.model.mesh.setSize([(0, ptag)], lc)

        gmsh.model.mesh.generate(2)
        sta_p, sta_t, sta_surf_tags = _extract_2d_mesh_from_gmsh()

        sta_reg = np.full(sta_t.shape[1], REG_STATOR_STEEL, dtype=np.int32)
        for st_id in np.unique(sta_surf_tags):
            mask = (sta_surf_tags == st_id)
            if st_id in cu_surf_to_idx:
                sta_reg[mask] = REG_SLOT_BASE + cu_surf_to_idx[st_id]
            elif st_id in air_surfs:
                sta_reg[mask] = REG_AIR_GAP_STA
            else:
                sta_reg[mask] = REG_STATOR_STEEL

        sta_reg = _apply_cut_edge_tagging(sta_p, sta_t, sta_reg, REG_STATOR_STEEL, REG_STATOR_STEEL_CUT)

        sta_g3_nodes = _find_arc_nodes(sta_p, r_g3)
        if len(sta_g3_nodes) != N_THETA_GAP + 1:
            raise RuntimeError(f"Stator r_g3 node count {len(sta_g3_nodes)} != {N_THETA_GAP + 1}")
        sta_p, sta_t, sta_reg, sta_gap_mid, sta_arkkio_elems = _append_structured_ring(
            sta_p, sta_t, sta_reg, sta_g3_nodes, r_mid, th_pole, REG_AIR_GAP_STA, outward=False
        )
        sta_per_0 = _find_radial_nodes(sta_p, 0.0)
        sta_per_p = _find_radial_nodes(sta_p, th_pole)
        _verify_periodic_pair(sta_p, sta_per_0, sta_per_p, "Stator")
        sta_outer_nodes = _find_arc_nodes(sta_p, R_SO)

        # =====================================================================
        # PART 2: BUILD ROTOR MESH (R_SH to r_g1) + APPEND [r_g1, r_mid] RING
        # =====================================================================
        gmsh.model.add("rotor")
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.00022 * mesh_scale)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.0044 * mesh_scale)

        rot_iron_sector = _make_occ_sector(R_SH, r_ro, 0.0, th_pole)
        rot_air_inner = _make_occ_sector(r_ro, r_g1, 0.0, th_pole)

        # Optional rotor surface flux-modulation notches
        notch_surfs = []
        if design.notch_depth_mm >= 0.15:
            nd = design.notch_depth_mm * 1e-3
            nw = math.radians(design.notch_width_deg)
            for sgn in (-1.0, +1.0):
                th_n = th_center + sgn * (design.notch_angle_frac * th_center)
                if 0.012 < th_n - 0.5 * nw and th_n + 0.5 * nw < th_pole - 0.012:
                    pts_n = [
                        ((r_ro - nd) * math.cos(th_n - 0.5 * nw), (r_ro - nd) * math.sin(th_n - 0.5 * nw)),
                        ((r_ro + 0.04 * gap) * math.cos(th_n - 0.5 * nw), (r_ro + 0.04 * gap) * math.sin(th_n - 0.5 * nw)),
                        ((r_ro + 0.04 * gap) * math.cos(th_n + 0.5 * nw), (r_ro + 0.04 * gap) * math.sin(th_n + 0.5 * nw)),
                        ((r_ro - nd) * math.cos(th_n + 0.5 * nw), (r_ro - nd) * math.sin(th_n + 0.5 * nw)),
                    ]
                    notch_surfs.append(_make_occ_polygon(pts_n))

        r_avail_min = R_SH + 0.0035
        r_avail_max = r_ro - 0.0026
        span_r = r_avail_max - r_avail_min

        barrier_air_polys = []
        magnet_items: List[Tuple[int, Tuple[float, float]]] = []
        use_mag = (design.magnet_grade != "NONE")

        for b in design.barriers:
            u_b = r_avail_min + b.r_frac * span_r
            t_b = b.thick_mm * 1e-3
            br_m = b.bridge_mm * 1e-3
            r_clip = r_ro - br_m

            # Cap w_b safely inside pole sector (for 4-pole tan(45)=1.0, cap by chord inside r_clip!)
            w_sec_max = u_b * math.tan(th_center) - 0.0022 - 0.5 * t_b
            w_rad_max = math.sqrt(max(0.0, (r_clip - 0.0025) ** 2 - u_b ** 2))
            w_max = max(0.0020, min(w_sec_max, w_rad_max))
            w_b = float(np.clip(b.half_width_frac, 0.06, 0.72)) * w_max

            r_tip_ext = r_clip + 0.0020
            th_L = th_center - b.beta_L_frac * th_center
            th_R = th_center + b.beta_R_frac * th_center

            p_knee_L = _uv_to_xy(u_b, -w_b, th_center)
            p_knee_R = _uv_to_xy(u_b, +w_b, th_center)
            p_tip_L = (r_tip_ext * math.cos(th_L), r_tip_ext * math.sin(th_L))
            p_tip_R = (r_tip_ext * math.cos(th_R), r_tip_ext * math.sin(th_R))

            seg_c = _make_occ_polygon(_rect_poly_xy(p_knee_L, p_knee_R, t_b))
            seg_L = _make_occ_polygon(_rect_poly_xy(p_knee_L, p_tip_L, t_b))
            seg_R = _make_occ_polygon(_rect_poly_xy(p_knee_R, p_tip_R, t_b))
            fused, _ = gmsh.model.occ.fuse([(2, seg_c)], [(2, seg_L), (2, seg_R)])

            dth_edge = max(0.015, 0.0012 / r_clip)
            clip_sec = _make_occ_sector(R_SH + 0.0025, r_clip, dth_edge, th_pole - dth_edge)
            bar_clipped, _ = gmsh.model.occ.intersect(fused, [(2, clip_sec)])
            if not bar_clipped:
                continue

            if b.rib_mm >= 0.80:
                rib_w = b.rib_mm * 1e-3
                p_r1 = _uv_to_xy(u_b - 1.5 * t_b, 0.0, th_center)
                p_r2 = _uv_to_xy(u_b + 1.5 * t_b, 0.0, th_center)
                rib_rect = _make_occ_polygon(_rect_poly_xy(p_r1, p_r2, rib_w))
                bar_clipped, _ = gmsh.model.occ.cut(bar_clipped, [(2, rib_rect)])

            if not bar_clipped:
                continue
            barrier_air_polys.extend([s for (_, s) in bar_clipped])

            if use_mag:
                t_m = 0.96 * t_b
                mx_d, my_d = math.cos(th_center), math.sin(th_center)
                if b.mag_center_fill >= 0.12 and w_b > 0.0022:
                    if b.rib_mm >= 0.80:
                        half_rib = 0.5 * b.rib_mm * 1e-3 + 0.00025
                        if w_b - half_rib > 0.0015:
                            w_span = (w_b - half_rib) * b.mag_center_fill
                            for sgn in (-1.0, +1.0):
                                v0 = sgn * half_rib
                                v1 = sgn * (half_rib + w_span)
                                pts_m = _rect_poly_xy(
                                    _uv_to_xy(u_b, min(v0, v1), th_center),
                                    _uv_to_xy(u_b, max(v0, v1), th_center), t_m
                                )
                                if pts_m:
                                    m_raw = _make_occ_polygon(pts_m)
                                    bc_copy = gmsh.model.occ.copy(bar_clipped)
                                    m_in, _ = gmsh.model.occ.intersect([(2, m_raw)], bc_copy)
                                    for (_, ms) in m_in:
                                        magnet_items.append((ms, (mx_d, my_d)))
                    else:
                        w_span = w_b * b.mag_center_fill
                        pts_m = _rect_poly_xy(
                            _uv_to_xy(u_b, -w_span, th_center),
                            _uv_to_xy(u_b, +w_span, th_center), t_m
                        )
                        if pts_m:
                            m_raw = _make_occ_polygon(pts_m)
                            bc_copy = gmsh.model.occ.copy(bar_clipped)
                            m_in, _ = gmsh.model.occ.intersect([(2, m_raw)], bc_copy)
                            for (_, ms) in m_in:
                                magnet_items.append((ms, (mx_d, my_d)))

                if b.mag_wing_fill >= 0.12:
                    for p_knee, p_tip in ((p_knee_L, p_tip_L), (p_knee_R, p_tip_R)):
                        dx = p_tip[0] - p_knee[0]
                        dy = p_tip[1] - p_knee[1]
                        L_w = math.hypot(dx, dy)
                        if L_w < 0.003:
                            continue
                        f0 = 0.08
                        f1 = min(0.88, f0 + 0.78 * b.mag_wing_fill)
                        q0 = (p_knee[0] + f0 * dx, p_knee[1] + f0 * dy)
                        q1 = (p_knee[0] + f1 * dx, p_knee[1] + f1 * dy)
                        nx, ny = -dy / L_w, dx / L_w
                        if nx * mx_d + ny * my_d < 0.0:
                            nx, ny = -nx, -ny
                        pts_m = _rect_poly_xy(q0, q1, t_m)
                        if pts_m:
                            m_raw = _make_occ_polygon(pts_m)
                            bc_copy = gmsh.model.occ.copy(bar_clipped)
                            m_in, _ = gmsh.model.occ.intersect([(2, m_raw)], bc_copy)
                            for (_, ms) in m_in:
                                magnet_items.append((ms, (nx, ny)))

        tools_r = [(2, rot_air_inner)]
        n_notch = len(notch_surfs)
        n_bar = len(barrier_air_polys)
        n_mag = len(magnet_items)
        tools_r += [(2, s) for s in notch_surfs]
        tools_r += [(2, s) for s in barrier_air_polys]
        tools_r += [(2, s) for (s, _) in magnet_items]

        _, r_map = gmsh.model.occ.fragment([(2, rot_iron_sector)], tools_r)
        gmsh.model.occ.synchronize()

        rot_air_gap_surfs = set(st for (_, st) in r_map[1])
        rot_notch_surfs = set()
        for k in range(n_notch):
            for (_, st) in r_map[2 + k]:
                rot_notch_surfs.add(st)
        rot_bar_surfs = set()
        for k in range(n_bar):
            for (_, st) in r_map[2 + n_notch + k]:
                rot_bar_surfs.add(st)

        rot_mag_surf_to_tag: Dict[int, int] = {}
        rot_mag_dir: Dict[int, Tuple[float, float]] = {}
        for k in range(n_mag):
            m_tag = REG_MAGNET_BASE + k
            _, m_vec = magnet_items[k]
            rot_mag_dir[m_tag] = m_vec
            for (_, st) in r_map[2 + n_notch + n_bar + k]:
                rot_mag_surf_to_tag[st] = m_tag

        rc0_list, rcp_list = [], []
        for dim, ctag in gmsh.model.getEntities(1):
            bnds = gmsh.model.getBoundary([(1, ctag)], oriented=False)
            if len(bnds) == 2:
                pA = gmsh.model.getValue(0, bnds[0][1], [])
                pB = gmsh.model.getValue(0, bnds[1][1], [])
                rA, rB = math.hypot(pA[0], pA[1]), math.hypot(pB[0], pB[1])
                thA, thB = math.atan2(pA[1], pA[0]), math.atan2(pB[1], pB[0])
                if abs(rA - r_g1) < 2.5e-5 and abs(rB - r_g1) < 2.5e-5 and abs(thB - thA) > 0.5 * th_pole:
                    gmsh.model.mesh.setTransfiniteCurve(ctag, N_THETA_GAP + 1)
                elif on_ray(pA, 0.0) and on_ray(pB, 0.0) and abs(rB - rA) > 1e-5:
                    rc0_list.append((0.5 * (rA + rB), ctag))
                elif on_ray(pA, th_pole) and on_ray(pB, th_pole) and abs(rB - rA) > 1e-5:
                    rcp_list.append((0.5 * (rA + rB), ctag))

        rc0_list.sort()
        rcp_list.sort()
        for (_, t0), (_, tp) in zip(rc0_list, rcp_list):
            gmsh.model.mesh.setPeriodic(1, [tp], [t0], rot_m)

        for dim, ptag in gmsh.model.getEntities(0):
            xyz = gmsh.model.getValue(0, ptag, [])
            rr = math.hypot(xyz[0], xyz[1])
            if rr > r_ro - 0.0055:
                lc = 0.00095 * mesh_scale
            elif rr > R_SH + 0.005:
                lc = 0.00195 * mesh_scale
            else:
                lc = 0.0034 * mesh_scale
            gmsh.model.mesh.setSize([(0, ptag)], lc)

        gmsh.model.mesh.generate(2)
        rot_p, rot_t, rot_surf_tags = _extract_2d_mesh_from_gmsh()

        rot_reg = np.full(rot_t.shape[1], REG_ROTOR_STEEL, dtype=np.int32)
        for st_id in np.unique(rot_surf_tags):
            mask = (rot_surf_tags == st_id)
            if st_id in rot_mag_surf_to_tag:
                rot_reg[mask] = rot_mag_surf_to_tag[st_id]
            elif st_id in rot_air_gap_surfs or st_id in rot_notch_surfs:
                rot_reg[mask] = REG_AIR_GAP_ROT
            elif st_id in rot_bar_surfs:
                rot_reg[mask] = REG_AIR_BARRIER
            else:
                rot_reg[mask] = REG_ROTOR_STEEL

        rot_reg = _apply_cut_edge_tagging(rot_p, rot_t, rot_reg, REG_ROTOR_STEEL, REG_ROTOR_STEEL_CUT)

        rot_g1_nodes = _find_arc_nodes(rot_p, r_g1)
        if len(rot_g1_nodes) != N_THETA_GAP + 1:
            raise RuntimeError(f"Rotor r_g1 node count {len(rot_g1_nodes)} != {N_THETA_GAP + 1}")
        rot_p, rot_t, rot_reg, rot_gap_mid, rot_arkkio_elems = _append_structured_ring(
            rot_p, rot_t, rot_reg, rot_g1_nodes, r_mid, th_pole, REG_AIR_GAP_ROT, outward=True
        )
        rot_per_0_all = _find_radial_nodes(rot_p, 0.0)
        rot_per_p_all = _find_radial_nodes(rot_p, th_pole)
        rot_per_0 = rot_per_0_all[np.hypot(rot_p[0, rot_per_0_all], rot_p[1, rot_per_0_all]) < r_mid - 1e-6]
        rot_per_p = rot_per_p_all[np.hypot(rot_p[0, rot_per_p_all], rot_p[1, rot_per_p_all]) < r_mid - 1e-6]
        _verify_periodic_pair(rot_p, rot_per_0, rot_per_p, "Rotor")
        rot_shaft_nodes = _find_arc_nodes(rot_p, R_SH)

        def tri_areas(pp, tt):
            return 0.5 * ((pp[0, tt[1]] - pp[0, tt[0]]) * (pp[1, tt[2]] - pp[1, tt[0]])
                          - (pp[0, tt[2]] - pp[0, tt[0]]) * (pp[1, tt[1]] - pp[1, tt[0]]))

        a_sta = tri_areas(sta_p, sta_t)
        a_rot = tri_areas(rot_p, rot_t)
        slot_areas = np.array([np.sum(a_sta[sta_reg == (REG_SLOT_BASE + k)]) for k in range(n_spp)])
        area_sta_steel = float(np.sum(a_sta[(sta_reg == REG_STATOR_STEEL) | (sta_reg == REG_STATOR_STEEL_CUT)]))
        area_rot_steel = float(np.sum(a_rot[(rot_reg == REG_ROTOR_STEEL) | (rot_reg == REG_ROTOR_STEEL_CUT)]))
        area_magnets = float(np.sum(a_rot[rot_reg >= REG_MAGNET_BASE]))
        area_copper_1pole = float(np.sum(slot_areas)) * K_FILL_MAX

        return MotorMeshBundle(
            p_pairs=design.p_pairs, q_spp=design.q_spp, slots_per_pole=n_spp,
            theta_pole=th_pole, theta_center=th_center,
            rot_p=rot_p, rot_t=rot_t, rot_reg=rot_reg, rot_mag_dir=rot_mag_dir,
            rot_per_0=rot_per_0, rot_per_p=rot_per_p, rot_gap_mid=rot_gap_mid,
            rot_shaft_nodes=rot_shaft_nodes, rot_arkkio_elems=rot_arkkio_elems,
            sta_p=sta_p, sta_t=sta_t, sta_reg=sta_reg,
            sta_per_0=sta_per_0, sta_per_p=sta_per_p, sta_gap_mid=sta_gap_mid,
            sta_outer_nodes=sta_outer_nodes, sta_arkkio_elems=sta_arkkio_elems,
            r_ro=r_ro, r_si=r_si, r_g1=r_g1, r_mid=r_mid, r_g3=r_g3,
            slot_areas=slot_areas, area_rot_steel=area_rot_steel,
            area_sta_steel=area_sta_steel, area_magnets=area_magnets,
            area_copper_1pole=area_copper_1pole, min_bridge_mm=min_br, min_rib_mm=min_rib,
        )
    finally:
        gmsh.clear()


# ---------------------------------------------------------------------------
# Industry & Textbook Baseline Presets (Across 8p, 6p, 4p)
# ---------------------------------------------------------------------------
def preset_ndfeb_vipm() -> MotorDesign:
    """B1: Commercial Tesla/Prius-class 8-pole 48-slot NdFeB V-shape IPM motor."""
    return MotorDesign(
        name="B1_NdFeB_VIPM_8p48s",
        magnet_grade="NDFEB_N42SH", p_pairs=4, q_spp=2,
        R_ro_mm=64.0, gap_mm=0.75, k_tooth=0.49, h_yoke_mm=14.0,
        b_so_mm=2.5, h_tt_mm=1.2, N_turns=3.5,
        barriers=[
            BarrierSpec(r_frac=0.34, thick_mm=5.2, half_width_frac=0.14,
                        beta_L_frac=0.76, beta_R_frac=0.76, bridge_mm=1.25,
                        rib_mm=1.8, mag_center_fill=0.0, mag_wing_fill=0.80),
            BarrierSpec(r_frac=0.76, thick_mm=2.6, half_width_frac=0.22,
                        beta_L_frac=0.39, beta_R_frac=0.39, bridge_mm=1.15,
                        rib_mm=0.0, mag_center_fill=0.0, mag_wing_fill=0.0),
        ],
    )


def preset_pure_synrm_4p() -> MotorDesign:
    """B2: Textbook 4-pole 36-slot 3-barrier Magnet-Free SynRM (ABB/Vagati industrial standard)."""
    return MotorDesign(
        name="B2_Pure_SynRM_4p36s",
        magnet_grade="NONE", p_pairs=2, q_spp=3,
        R_ro_mm=59.5, gap_mm=0.65, k_tooth=0.48, h_yoke_mm=22.0,
        b_so_mm=2.5, h_tt_mm=1.1, N_turns=4.0,
        barriers=[
            BarrierSpec(r_frac=0.20, thick_mm=4.8, half_width_frac=0.42,
                        beta_L_frac=0.86, beta_R_frac=0.86, bridge_mm=1.05,
                        rib_mm=1.2, mag_center_fill=0.0, mag_wing_fill=0.0),
            BarrierSpec(r_frac=0.50, thick_mm=4.2, half_width_frac=0.36,
                        beta_L_frac=0.62, beta_R_frac=0.62, bridge_mm=1.00,
                        rib_mm=1.0, mag_center_fill=0.0, mag_wing_fill=0.0),
            BarrierSpec(r_frac=0.79, thick_mm=3.2, half_width_frac=0.26,
                        beta_L_frac=0.36, beta_R_frac=0.36, bridge_mm=0.95,
                        rib_mm=0.0, mag_center_fill=0.0, mag_wing_fill=0.0),
        ],
    )


def preset_ferrite_pmasynrm_6p() -> MotorDesign:
    """B3: Textbook 6-pole 36-slot 3-barrier Symmetric Ferrite PMa-SynRM."""
    return MotorDesign(
        name="B3_Ferrite_PMaSynRM_6p36s",
        magnet_grade="FERRITE_FB12B", p_pairs=3, q_spp=2,
        R_ro_mm=62.5, gap_mm=0.65, k_tooth=0.48, h_yoke_mm=17.0,
        b_so_mm=2.5, h_tt_mm=1.1, N_turns=4.5,
        barriers=[
            BarrierSpec(r_frac=0.20, thick_mm=5.4, half_width_frac=0.44,
                        beta_L_frac=0.86, beta_R_frac=0.86, bridge_mm=1.15,
                        rib_mm=1.4, mag_center_fill=0.88, mag_wing_fill=0.78),
            BarrierSpec(r_frac=0.51, thick_mm=4.4, half_width_frac=0.38,
                        beta_L_frac=0.61, beta_R_frac=0.61, bridge_mm=1.10,
                        rib_mm=1.1, mag_center_fill=0.88, mag_wing_fill=0.72),
            BarrierSpec(r_frac=0.80, thick_mm=3.2, half_width_frac=0.28,
                        beta_L_frac=0.35, beta_R_frac=0.35, bridge_mm=1.00,
                        rib_mm=0.0, mag_center_fill=0.88, mag_wing_fill=0.0),
        ],
    )


def preset_ferrite_pmasynrm_4p() -> MotorDesign:
    """B4: Textbook 4-pole 36-slot 4-barrier Ferrite PMa-SynRM (high-saliency flux-concentrating)."""
    return MotorDesign(
        name="B4_Ferrite_PMaSynRM_4p36s",
        magnet_grade="FERRITE_FB12B", p_pairs=2, q_spp=3,
        R_ro_mm=60.0, gap_mm=0.65, k_tooth=0.48, h_yoke_mm=21.5,
        b_so_mm=2.5, h_tt_mm=1.1, N_turns=4.0,
        barriers=[
            BarrierSpec(r_frac=0.14, thick_mm=4.4, half_width_frac=0.42,
                        beta_L_frac=0.89, beta_R_frac=0.89, bridge_mm=1.15,
                        rib_mm=1.4, mag_center_fill=0.88, mag_wing_fill=0.80),
            BarrierSpec(r_frac=0.39, thick_mm=4.0, half_width_frac=0.38,
                        beta_L_frac=0.70, beta_R_frac=0.70, bridge_mm=1.10,
                        rib_mm=1.2, mag_center_fill=0.88, mag_wing_fill=0.76),
            BarrierSpec(r_frac=0.63, thick_mm=3.5, half_width_frac=0.32,
                        beta_L_frac=0.50, beta_R_frac=0.50, bridge_mm=1.05,
                        rib_mm=1.0, mag_center_fill=0.88, mag_wing_fill=0.68),
            BarrierSpec(r_frac=0.85, thick_mm=2.6, half_width_frac=0.24,
                        beta_L_frac=0.28, beta_R_frac=0.28, bridge_mm=0.95,
                        rib_mm=0.0, mag_center_fill=0.85, mag_wing_fill=0.0),
        ],
    )

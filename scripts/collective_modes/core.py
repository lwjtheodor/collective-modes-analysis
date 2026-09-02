"""Geometry, correlation, and mode construction shared by every command."""

from __future__ import annotations

import numpy as np

from .dump import Frame


def stable_order(frame: Frame) -> Frame:
    if "id" not in frame.fields:
        return frame
    order = np.argsort(frame.column("id").astype(np.int64))
    return Frame(frame.timestep, frame.bounds, frame.fields, frame.values[order])


def coordinates(frame: Frame, unwrapped: bool = False) -> np.ndarray:
    xyz = np.column_stack([frame.column(name) for name in ("x", "y", "z")]).astype(float)
    if unwrapped and {"ix", "iy", "iz"}.issubset(frame.fields):
        images = np.column_stack([frame.column(name) for name in ("ix", "iy", "iz")])
        xyz += images * frame.box_lengths
    return xyz


def velocities(frame: Frame) -> np.ndarray:
    return np.column_stack([frame.column(name) for name in ("vx", "vy", "vz")]).astype(float)


def cylindrical_basis(xyz: np.ndarray, axis_xy: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    dxy = xyz[:, :2] - axis_xy[None, :]
    radius = np.hypot(dxy[:, 0], dxy[:, 1])
    er = dxy / np.maximum(radius[:, None], 1e-12)
    etheta = np.column_stack((-er[:, 1], er[:, 0]))
    theta = np.arctan2(dxy[:, 1], dxy[:, 0])
    return radius, theta, np.stack((er, etheta), axis=-1)


def velocity_in_frame(velocity: np.ndarray, mode: str, wall_velocity: np.ndarray | None = None) -> np.ndarray:
    if mode == "lab":
        return velocity
    if mode == "selected_com":
        return velocity - velocity.mean(axis=0, keepdims=True)
    if mode == "wall_relative":
        if wall_velocity is None:
            raise ValueError("wall_relative frame requires CNT atom velocities")
        return velocity - wall_velocity[None, :]
    raise ValueError(f"unknown velocity frame {mode!r}")


def cylindrical_currents(xyz: np.ndarray, velocity: np.ndarray, lz_A: float, axis_xy: np.ndarray, rcnt_A: float, n_values: np.ndarray, m_values: np.ndarray) -> dict[str, np.ndarray]:
    """Return complex currents indexed [n,m] for Jz, Jr, Jtheta, L, Tinplane, Tr."""
    _, theta, basis = cylindrical_basis(xyz, axis_xy)
    er, etheta = basis[:, :, 0], basis[:, :, 1]
    vr = np.sum(velocity[:, :2] * er, axis=1)
    vtheta = np.sum(velocity[:, :2] * etheta, axis=1)
    kz = 2.0 * np.pi * n_values / lz_A
    ktheta = m_values / rcnt_A
    phase = np.exp(-1j * (kz[:, None, None] * xyz[None, None, :, 2] + ktheta[None, :, None] * theta[None, None, :]))
    jz = np.einsum("nmp,p->nm", phase, velocity[:, 2])
    jr = np.einsum("nmp,p->nm", phase, vr)
    jtheta = np.einsum("nmp,p->nm", phase, vtheta)
    q = np.hypot(kz[:, None], ktheta[None, :])
    long = np.divide(kz[:, None] * jz + ktheta[None, :] * jtheta, q, out=np.zeros_like(jz), where=q > 0)
    tinplane = np.divide(-ktheta[None, :] * jz + kz[:, None] * jtheta, q, out=np.zeros_like(jz), where=q > 0)
    return {"Jz": jz, "Jr": jr, "Jtheta": jtheta, "L": long, "Tinplane": tinplane, "Tr": jr, "kz_inv_A": kz, "ktheta_inv_A": ktheta}


def complex_acf(series: np.ndarray, max_lag: int, demean: bool = True) -> np.ndarray:
    """All-origin ordered Re<delta A(t+tau) delta A(t)*> for last axis channels."""
    values = np.asarray(series, dtype=np.complex128)
    if demean:
        values = values - values.mean(axis=0, keepdims=True)
    ntime = values.shape[0]
    if not 0 <= max_lag < ntime:
        raise ValueError(f"max_lag={max_lag} incompatible with {ntime} frames")
    size = 1 << (2 * ntime - 1).bit_length()
    transform = np.fft.fft(values, size, axis=0)
    raw = np.fft.ifft(transform * np.conjugate(transform), axis=0)[: max_lag + 1]
    return raw / (ntime - np.arange(max_lag + 1))[:, None]


def complex_cross_acf(a: np.ndarray, b: np.ndarray, max_lag: int) -> np.ndarray:
    """Ordered Re<delta A(t+tau) delta B(t)*>; C_AB and C_BA must be called separately."""
    a = np.asarray(a, dtype=np.complex128) - np.mean(a, axis=0, keepdims=True)
    b = np.asarray(b, dtype=np.complex128) - np.mean(b, axis=0, keepdims=True)
    ntime = len(a)
    if len(b) != ntime or not 0 <= max_lag < ntime:
        raise ValueError("cross-correlation shapes/max lag incompatible")
    size = 1 << (2 * ntime - 1).bit_length()
    raw = np.fft.ifft(np.fft.fft(a, size, axis=0) * np.conjugate(np.fft.fft(b, size, axis=0)), axis=0)[: max_lag + 1]
    return raw / (ntime - np.arange(max_lag + 1))[:, None]


def integrate_vacf(vacf: np.ndarray, dt_ps: float) -> tuple[np.ndarray, np.ndarray]:
    """Return consistent 1D MSD and alpha=t I/J from a dimensional VACF."""
    time = np.arange(len(vacf), dtype=float) * dt_ps
    integral = np.concatenate(([0.0], np.cumsum((vacf[1:] + vacf[:-1]) * 0.5 * dt_ps)))
    second = np.concatenate(([0.0], np.cumsum((integral[1:] + integral[:-1]) * 0.5 * dt_ps)))
    msd = 2.0 * second
    alpha = np.full_like(time, np.nan)
    valid = (time > 0) & (second != 0)
    alpha[valid] = time[valid] * integral[valid] / second[valid]
    return msd, alpha

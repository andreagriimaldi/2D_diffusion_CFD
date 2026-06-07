# CFD Examen Parcial 
## 2D Steady-State Heat Diffusion by the Finite Volume Method

FIE-UNDEF · Introducción a la Dinámica de Fluidos Computacional · 2026-S1

---

## 1. Problem statement

Metal plate `10 × 10 × 1 cm` containing two metallic sectors:

| Sector | Size | `k` [W/m/K] | Generation `q̇` [W/m³] |
|--------|------|-------------|-----------------------|
| A      | 4 × 4 × 1 cm | 100   | 5.0 × 10⁵ |
| B      | remainder    | 550   | 0 |

Governing equation (pure diffusion, no convection):

$$\frac{\partial}{\partial x}\left(k\,\frac{\partial T}{\partial x}\right) + \frac{\partial}{\partial y}\left(k\,\frac{\partial T}{\partial y}\right) + \dot{q} = 0$$

### Boundary conditions

| Edge | Type | Condition |
|------|------|-----------|
| North | Neumann | insulated, ∂T/∂n = 0 |
| South | Dirichlet | T = 30 °C |
| East / West | Robin (convective) | −k ∂T/∂n = h (T − T∞), h = 200 W/m²/K, T∞ = 15 °C |
| Front / Back (z) | Neumann | insulated → 2D problem |

### Configurations

- **Fig. 1 (base):** sector A in the top-left corner.
- **Fig. 2 (alt):** sector A moved to the bottom-right interior.

---

## 2. Tasks

- [ ] Solve `T(x,y)` by FVM (Fig. 1).
- [ ] Contour plot of `T(x,y)`.
- [ ] Temperature profile along the diagonal crossing sector A.
- [ ] Solve `T(x,y)` for Fig. 2; contour plot.
- [ ] Energy balance: does the system hold the same total energy in each configuration? Justify with calculations from the results.

---

## 3. Method

---

## 4. Repository layout



---

## 5. Parameters

| Symbol | Value | Unit | Notes |
|--------|-------|------|-------|
| Lx, Ly | 0.10  | m | plate side |
| t      | 0.01  | m | thickness |
| k_A    | 100   | W/m/K | |
| k_B    | 550   | W/m/K | |
| q̇_A   | 5.0e5 | W/m³ | sector A only |
| T_south| 30    | °C | Dirichlet |
| T∞     | 15    | °C | E/W fluid |
| h      | 200   | W/m²/K | E/W convection |
| Nx, Ny | <!-- TODO --> | — | mesh resolution |

---

## 6. How to run

```bash
python TO DO
```

---

## 7. Results

---

## 8. Notes / pitfalls


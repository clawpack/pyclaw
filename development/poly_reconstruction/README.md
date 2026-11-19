# Polynomial reconstruction: generator and checks

Scripts behind the SharpClaw polynomial reconstructions (`lim_type` 4, 5
and 6, selected with `solver.poly_order`). The write-up of the method,
with figures, belongs in the Clawpack documentation.

| File | Purpose |
|---|---|
| `poly_weights.py` | Exact (rational) weights and the generator of `src/pyclaw/sharpclaw/poly.f90`. `src/pyclaw/tests/test_sharpclaw_poly.py` uses it to check that `poly.f90` is up to date. |
| `sine_reconstruction.py` | Reconstructs sin(2πx) from N cell averages with every type and order and prints the error and observed order; `--fortran` compares with the compiled code. |
| `test_poly_reconstruction.py` | Toy tests (accuracy, polynomial exactness, dissipation). Run with `python -m pytest` in this directory; CI does not collect `development/`. |

```bash
python poly_weights.py --check   # exit 1 if poly.f90 is out of date
python poly_weights.py --write   # regenerate poly.f90 after changing the weights
python poly_weights.py --table   # print the weights as Markdown tables
python sine_reconstruction.py --cells 16
```

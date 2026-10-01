#!/usr/bin/env python3
"""Compare MATLAB/Octave ``laplacian_nd.m`` with SciPy ``LaplacianNd``.

The dynamic tests invoke ``laplacian_nd`` through GNU Octave or MATLAB and
compare it with ``scipy.sparse.linalg.LaplacianNd`` for uniform Dirichlet,
Neumann, and periodic boundaries. They cover:

* 1D grids of several sizes for each supported boundary condition.
* 2D grids with distinct axis sizes for each supported boundary condition.
* Two 3D grids and 4D/5D tensor-product grids.
* Dense and sparse matrix equality, plus matrix-free ``matvec`` and ``matmat``
  actions. The SciPy grid shape is reversed to align its C-order indexing with
  MATLAB/Octave's column-major indexing, and its matrix is negated because the
  two implementations use opposite Laplacian signs.
* Smallest eigenvalue equality, eigenvector orthonormality, and eigenpair
  residuals for both implementations.
* Permutation similarity when SciPy receives the unreversed grid shape.
* Pure-boundary entries in ``laplacian_reference.json`` and two SciPy edge
  cases: one-point periodic grids and anisotropic partial-spectrum requests.

Set ``OCTAVE_EXECUTABLE`` or ``MATLAB_EXECUTABLE`` to choose a runner; otherwise
the script searches the PATH and common Windows install locations. Run with
``python test_laplacian_nd_vs_scipy.py`` or
``pytest test_laplacian_nd_vs_scipy.py``.
"""

import json
import os
import shutil
import subprocess
import sys
import unittest
import numpy as np
from scipy.sparse.linalg import LaplacianNd

# Mapping between SciPy boundary conditions and laplacian_nd boundary codes
BC_SCIPY_TO_MATLAB = {
    "dirichlet": "DD",
    "neumann": "NN",
    "periodic": "P",
}

BC_MATLAB_TO_SCIPY = {
    "DD": "dirichlet",
    "NN": "neumann",
    "P": "periodic",
}


def find_octave_or_matlab():
    """Locate GNU Octave or MATLAB executable."""
    env_octave = os.environ.get("OCTAVE_EXECUTABLE")
    if env_octave and os.path.exists(env_octave):
        return env_octave, "octave"

    env_matlab = os.environ.get("MATLAB_EXECUTABLE")
    if env_matlab and os.path.exists(env_matlab):
        return env_matlab, "matlab"

    for candidate in ["octave-cli", "octave"]:
        path = shutil.which(candidate)
        if path:
            return path, "octave"

    path = shutil.which("matlab")
    if path:
        return path, "matlab"

    # Common Windows installation locations
    candidate_patterns = [
        r"J:\Program Files\GNU Octave\Octave-11.3.0\mingw64\bin\octave-cli.exe",
        r"C:\Program Files\GNU Octave\Octave-11.3.0\mingw64\bin\octave-cli.exe",
        r"C:\Program Files\GNU Octave\*\mingw64\bin\octave-cli.exe",
        r"C:\Program Files\MATLAB\*\bin\matlab.exe",
    ]
    import glob
    for pattern in candidate_patterns:
        matches = glob.glob(pattern)
        if matches:
            kind = "octave" if "octave" in matches[0].lower() else "matlab"
            return matches[0], kind

    return None, None


OCTAVE_EXE, RUNNER_KIND = find_octave_or_matlab()
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
LAPLACIAN_DIR = SCRIPT_DIR
FIXTURE_PATH = os.path.join(LAPLACIAN_DIR, "laplacian_reference.json")


def run_laplacian_nd(N, B, M=0):
    """
    Invoke laplacian_nd(N, B, M) via Octave / MATLAB.
    Returns:
        A: (prod(N), prod(N)) np.ndarray
        lam: (M,) np.ndarray
        V: (prod(N), M) np.ndarray
    """
    if not OCTAVE_EXE:
        raise RuntimeError("Neither Octave nor MATLAB executable could be found.")

    N_str = "[" + " ".join(str(n) for n in N) + "]"
    B_str = "{" + ", ".join(f"'{b}'" for b in B) + "}"

    # Build Octave / MATLAB command that outputs JSON
    octave_code = f"""
    addpath('{LAPLACIAN_DIR.replace(os.sep, "/")}');
    [A, lam, V] = laplacian_nd({N_str}, {B_str}, {M});
    out = struct();
    out.A = full(A);
    out.lambda = lam;
    out.V = V;
    printf('__JSON_START__%s__JSON_END__', jsonencode(out));
    """

    if RUNNER_KIND == "octave":
        cmd = [OCTAVE_EXE, "--no-gui", "--no-history", "--no-init-file", "--eval", octave_code]
    else:
        cmd = [OCTAVE_EXE, "-nosplash", "-nodesktop", "-batch", octave_code]

    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    stdout = result.stdout
    start_tag = "__JSON_START__"
    end_tag = "__JSON_END__"
    start_idx = stdout.find(start_tag)
    end_idx = stdout.find(end_tag)
    if start_idx == -1 or end_idx == -1:
        raise RuntimeError(f"Failed to parse JSON output from runner:\nSTDOUT:\n{stdout}\nSTDERR:\n{result.stderr}")

    json_str = stdout[start_idx + len(start_tag):end_idx]
    data = json.loads(json_str)

    order = int(np.prod(N))
    A = np.array(data["A"], dtype=np.float64).reshape((order, order))
    lam = np.array(data["lambda"], dtype=np.float64).reshape(-1) if M > 0 else np.empty(0)
    V = np.array(data["V"], dtype=np.float64).reshape((order, M)) if M > 0 else np.empty((order, 0))
    return A, lam, V


class TestLaplacianNdVsSciPy(unittest.TestCase):
    """Test laplacian_nd.m against scipy.sparse.linalg.LaplacianNd."""

    @classmethod
    def setUpClass(cls):
        if not OCTAVE_EXE:
            print("WARNING: Octave/MATLAB not found. Dynamic runner tests will be skipped.")

    def _check_case(self, N, bc_scipy, M=None, atol=1e-12):
        """Core comparison between laplacian_nd and LaplacianNd."""
        dimension = len(N)
        bc_matlab = [BC_SCIPY_TO_MATLAB[bc_scipy]] * dimension
        order = int(np.prod(N))
        if M is None:
            M = min(order, 4)

        # 1. SciPy LaplacianNd
        # Grid shape in SciPy is reversed because MATLAB uses Fortran (column-major)
        # order where N[0] varies fastest, whereas SciPy uses C (row-major) order.
        grid_shape = tuple(reversed(N))
        lap = LaplacianNd(grid_shape, boundary_conditions=bc_scipy)

        # 2. MATLAB / Octave laplacian_nd
        A_matlab, lam_matlab, V_matlab = run_laplacian_nd(N, bc_matlab, M)

        # 3. Check matrix equality: A_matlab == -lap.toarray()
        A_scipy = -lap.toarray().astype(np.float64)
        np.testing.assert_allclose(
            A_matlab,
            A_scipy,
            atol=atol,
            err_msg=f"Matrix mismatch for N={N}, bc={bc_scipy}",
        )

        # Check sparse format matches
        A_scipy_sparse = -lap.tosparse().toarray().astype(np.float64)
        np.testing.assert_allclose(
            A_matlab,
            A_scipy_sparse,
            atol=atol,
            err_msg=f"Sparse matrix mismatch for N={N}, bc={bc_scipy}",
        )

        # 4. Check operator action: -lap.matvec(x) == A @ x, -lap.matmat(X) == A @ X
        rng = np.random.default_rng(42)
        x = rng.standard_normal(order)
        X = rng.standard_normal((order, 3))
        np.testing.assert_allclose(
            -lap.matvec(x),
            A_matlab @ x,
            atol=atol,
            err_msg=f"matvec mismatch for N={N}, bc={bc_scipy}",
        )
        np.testing.assert_allclose(
            -lap.matmat(X),
            A_matlab @ X,
            atol=atol,
            err_msg=f"matmat mismatch for N={N}, bc={bc_scipy}",
        )

        # 5. Check eigenvalues if M > 0
        if M > 0:
            # SciPy eigenvalues are m largest of Δ (ascending, <= 0)
            # MATLAB eigenvalues are M smallest of -Δ (ascending, >= 0)
            evals_scipy = -lap.eigenvalues(M)[::-1]
            np.testing.assert_allclose(
                lam_matlab,
                evals_scipy,
                atol=atol,
                err_msg=f"Eigenvalues mismatch for N={N}, bc={bc_scipy}, M={M}",
            )

        # 6. Check eigenvectors if M > 0
        if M > 0:
            # Orthonormality of MATLAB eigenvectors
            np.testing.assert_allclose(
                V_matlab.T @ V_matlab,
                np.eye(M),
                atol=atol,
                err_msg=f"MATLAB eigenvectors not orthonormal for N={N}, bc={bc_scipy}",
            )

            # Eigenpair residual for MATLAB eigenvectors: ||A*V - V*diag(lambda)||
            res_matlab = A_matlab @ V_matlab - V_matlab @ np.diag(lam_matlab)
            res_norm_matlab = np.linalg.norm(res_matlab, "fro")
            self.assertLess(
                res_norm_matlab,
                atol * order,
                f"MATLAB eigenvector residual too large: {res_norm_matlab}",
            )

            # SciPy eigenvectors
            V_scipy = lap.eigenvectors(M)[:, ::-1]
            # Orthonormality of SciPy eigenvectors
            np.testing.assert_allclose(
                V_scipy.T @ V_scipy,
                np.eye(M),
                atol=atol,
                err_msg=f"SciPy eigenvectors not orthonormal for N={N}, bc={bc_scipy}",
            )

            # Eigenpair residual for SciPy eigenvectors against A_matlab
            res_scipy = A_matlab @ V_scipy - V_scipy @ np.diag(lam_matlab)
            res_norm_scipy = np.linalg.norm(res_scipy, "fro")
            self.assertLess(
                res_norm_scipy,
                atol * order,
                f"SciPy eigenvector residual too large: {res_norm_scipy}",
            )

    # ---------------- 1D tests ----------------
    def test_1d_dirichlet(self):
        if not OCTAVE_EXE:
            self.skipTest("Octave/MATLAB not found")
        for n in [3, 4, 7]:
            self._check_case([n], "dirichlet", M=min(n, 3))

    def test_1d_neumann(self):
        if not OCTAVE_EXE:
            self.skipTest("Octave/MATLAB not found")
        for n in [3, 5, 8]:
            self._check_case([n], "neumann", M=min(n, 3))

    def test_1d_periodic(self):
        if not OCTAVE_EXE:
            self.skipTest("Octave/MATLAB not found")
        for n in [3, 6, 9]:
            self._check_case([n], "periodic", M=min(n, 3))

    # ---------------- 2D tests ----------------
    def test_2d_dirichlet(self):
        if not OCTAVE_EXE:
            self.skipTest("Octave/MATLAB not found")
        for N in [[2, 3], [3, 3], [4, 2]]:
            self._check_case(N, "dirichlet", M=min(int(np.prod(N)), 5))

    def test_2d_neumann(self):
        if not OCTAVE_EXE:
            self.skipTest("Octave/MATLAB not found")
        for N in [[2, 3], [3, 4], [2, 5]]:
            self._check_case(N, "neumann", M=min(int(np.prod(N)), 5))

    def test_2d_periodic(self):
        if not OCTAVE_EXE:
            self.skipTest("Octave/MATLAB not found")
        for N in [[2, 3], [3, 3], [4, 3]]:
            self._check_case(N, "periodic", M=min(int(np.prod(N)), 5))

    # ---------------- 3D tests ----------------
    def test_3d_all_boundary_conditions(self):
        if not OCTAVE_EXE:
            self.skipTest("Octave/MATLAB not found")
        for bc in ["dirichlet", "neumann", "periodic"]:
            self._check_case([2, 3, 4], bc, M=4)
            self._check_case([3, 2, 2], bc, M=4)

    # ---------------- 4D & 5D tests ----------------
    def test_higher_dimensions(self):
        if not OCTAVE_EXE:
            self.skipTest("Octave/MATLAB not found")
        # 4D case
        for bc in ["dirichlet", "neumann", "periodic"]:
            self._check_case([2, 2, 2, 2], bc, M=4)
        # 5D case
        self._check_case([2, 2, 2, 2, 2], "dirichlet", M=2)
        self._check_case([2, 2, 2, 2, 2], "neumann", M=2)

    # ---------------- Coordinate Permutation Test ----------------
    def test_unreversed_shape_permutation_similarity(self):
        """
        Verify that passing unreversed grid_shape=tuple(N) in SciPy yields a matrix
        that is permutation-similar to laplacian_nd(N):
            P @ A_matlab @ P.T == -lap_unreversed.toarray()
        and that eigenvalues are identical.
        """
        if not OCTAVE_EXE:
            self.skipTest("Octave/MATLAB not found")
        N = [2, 3, 4]
        A_matlab, lam_matlab, _ = run_laplacian_nd(N, ["DD", "DD", "DD"], M=4)

        # Unreversed grid_shape in SciPy
        lap_unreversed = LaplacianNd(tuple(N), boundary_conditions="dirichlet")
        A_scipy_unreversed = -lap_unreversed.toarray().astype(np.float64)

        # Eigenvalues must be identical regardless of coordinate ordering
        evals_scipy = -lap_unreversed.eigenvalues(4)[::-1]
        np.testing.assert_allclose(
            lam_matlab,
            evals_scipy,
            atol=1e-12,
            err_msg="Eigenvalues must be invariant to dimension ordering",
        )

        # Compute the permutation between Fortran and C flattening
        order = int(np.prod(N))
        indices_3d = np.arange(order).reshape(N, order="F")
        perm = indices_3d.ravel(order="C")
        P = np.eye(order)[perm]

        np.testing.assert_allclose(
            P @ A_matlab @ P.T,
            A_scipy_unreversed,
            atol=1e-12,
            err_msg="Matrix must be permutation-similar under reshape permutation",
        )

    # ---------------- Precomputed Fixture validation ----------------
    def test_against_reference_fixture(self):
        """
        Validate all pure BC cases from laplacian_reference.json against SciPy.
        Tests 178 cases across 1D, 2D, and 3D with Dirichlet, Neumann, and Periodic BCs.
        """
        if not os.path.exists(FIXTURE_PATH):
            self.skipTest(f"Reference fixture {FIXTURE_PATH} not found.")

        with open(FIXTURE_PATH, "r", encoding="utf-8") as f:
            fixture = json.load(f)

        tested_count = 0
        periodic_1d_edge_count = 0

        for case in fixture["cases"]:
            B = case["B"]
            # Only test cases where all dimensions have the same pure BC: 'DD', 'NN', or 'P'
            if not all(b == B[0] for b in B) or B[0] not in BC_MATLAB_TO_SCIPY:
                continue

            bc_scipy = BC_MATLAB_TO_SCIPY[B[0]]
            N = [case["N"]] if isinstance(case["N"], int) else list(case["N"])
            M = int(case["M"])
            order = int(np.prod(N))

            expected_A = np.array(case["A"], dtype=np.float64).reshape((order, order))
            expected_lambda = np.array(case["lambda"], dtype=np.float64).reshape(-1)

            # SciPy
            grid_shape = tuple(reversed(N))
            lap = LaplacianNd(grid_shape, boundary_conditions=bc_scipy)

            # Handle 1-point periodic edge case where SciPy toarray() sets -1 instead of 0:
            if bc_scipy == "periodic" and any(n == 1 for n in N):
                periodic_1d_edge_count += 1
                # SciPy's full spectrum eigenvalues match laplacian_nd
                evals_scipy = -lap.eigenvalues()[-M:][::-1]
                np.testing.assert_allclose(
                    expected_lambda,
                    evals_scipy,
                    atol=1e-12,
                    err_msg=f"Fixture eigenvalues mismatch for N={N}, B={B}, M={M}",
                )
                continue

            A_scipy = -lap.toarray().astype(np.float64)

            # Compare matrix
            np.testing.assert_allclose(
                expected_A,
                A_scipy,
                atol=1e-12,
                err_msg=f"Fixture matrix mismatch for N={N}, B={B}",
            )

            # Compare eigenvalues using full spectrum to avoid SciPy's
            # min(tuple) lexicographic bug when M > min(grid_shape)
            evals_scipy = -lap.eigenvalues()[-M:][::-1]
            np.testing.assert_allclose(
                expected_lambda,
                evals_scipy,
                atol=1e-12,
                err_msg=f"Fixture eigenvalues mismatch for N={N}, B={B}, M={M}",
            )
            tested_count += 1

        self.assertGreater(tested_count, 0)
        sys.stderr.write(
            f"\n[Fixture validation] Successfully verified {tested_count} pure-BC cases from laplacian_reference.json against SciPy LaplacianNd "
            f"({periodic_1d_edge_count} 1-point periodic edge cases handled).\n"
        )

    # ---------------- SciPy Edge Case / Quirk Tests ----------------
    def test_scipy_known_edge_cases(self):
        """
        Document and test known SciPy edge cases and quirks:
        1. 1-point periodic grid:
           laplacian_nd(1, {'P'}) produces A=[0] and lambda=[0].
           SciPy's LaplacianNd((1,), boundary_conditions='periodic') has toarray()=[[-1]]
           due to adding 1 once instead of twice for periodic boundaries, while its own
           eigenvalues() returns [0.].
        2. SciPy eigenvalues(m) lexicographical min bug when m > min(grid_shape):
           For grid_shape=(3, 1) and m=2, SciPy does min((3, 1), (2, 2)) -> (2, 2)
           which evaluates out-of-bounds index for dimension 1, whereas lap.eigenvalues()[-m:]
           correctly returns the true eigenvalues [-3., 0.].
        """
        # Edge case 1: 1-point periodic
        if OCTAVE_EXE:
            A_matlab, lam_matlab, _ = run_laplacian_nd([1], ["P"], M=1)
            self.assertEqual(A_matlab[0, 0], 0.0)
            self.assertEqual(lam_matlab[0], 0.0)

        lap_p1 = LaplacianNd((1,), boundary_conditions="periodic")
        # SciPy's own eigenvalues() returns 0.0 (matching laplacian_nd)
        self.assertAlmostEqual(lap_p1.eigenvalues()[0], 0.0, places=12)
        # But SciPy's toarray() returns -1
        self.assertEqual(lap_p1.toarray()[0, 0], -1)

        # Edge case 2: m > min(grid_shape)
        lap_aniso = LaplacianNd((3, 1), boundary_conditions="periodic")
        # Full spectrum correctly yields eigenvalues
        full_evals = lap_aniso.eigenvalues()
        np.testing.assert_allclose(full_evals, [-3.0, -3.0, 0.0], atol=1e-12)


def main():
    print("=" * 70)
    print("Testing BLOPEX laplacian_nd.m vs scipy.sparse.linalg.LaplacianNd")
    print("=" * 70)
    print(f"Runner: {RUNNER_KIND or 'None'} -> {OCTAVE_EXE or 'Not found'}")
    print(f"SciPy:  {__import__('scipy').__version__}")
    print(f"NumPy:  {np.__version__}")
    print("=" * 70)
    unittest.main(verbosity=2)


if __name__ == "__main__":
    main()

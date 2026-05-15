import subprocess
import sys
import textwrap
from pathlib import Path

import pytest
from conftest import TEST_DATA_DIR, original_xponge_repo


XPONGE_REPO = original_xponge_repo()
DATA_POLY_DIR = TEST_DATA_DIR / "poly"
EXPECTED_OUTPUT_NAMES = {
    "poly_LJ.txt",
    "poly_angle.txt",
    "poly_atom_name.txt",
    "poly_atom_type_name.txt",
    "poly_bond.txt",
    "poly_charge.txt",
    "poly_coordinate.txt",
    "poly_dihedral.txt",
    "poly_exclude.txt",
    "poly_mass.txt",
    "poly_nb14.txt",
    "poly_residue.txt",
    "poly_resname.txt",
}


def _snapshot_poly_outputs(path: Path) -> dict[str, str]:
    return {
        item.name: item.read_text()
        for item in sorted(path.iterdir())
        if item.is_file() and item.name in EXPECTED_OUTPUT_NAMES
    }


def _run_original_poly_workflow(output_dir: Path) -> dict[str, str]:
    if not XPONGE_REPO.exists():
        pytest.skip("local Xponge reference repository is not available")
    child_code = textwrap.dedent(
        f"""
        from pathlib import Path
        import sys

        sys.path.insert(0, {str(XPONGE_REPO)!r})

        import Xponge
        from Xponge.forcefield.opls import load_parameter_from_ffitp

        data_dir = Path({str(DATA_POLY_DIR)!r})
        output_dir = Path({str(output_dir)!r})

        def clean_gro_in_place(path: Path) -> None:
            lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
            if not lines:
                return
            last = lines[-1].strip().split()
            if len(last) == 9:
                lines[-1] = "   ".join(last[:3]) + "\\n"
                path.write_text("".join(lines), encoding="utf-8")

        gro_path = output_dir / "ps.gro"
        gro_path.write_text((data_dir / "ps.gro").read_text(encoding="utf-8"), encoding="utf-8")
        clean_gro_in_place(gro_path)
        itp_path = data_dir / "ps.itp"
        psf_path = data_dir / "ps.psf"

        load_parameter_from_ffitp(itp_path.name, str(itp_path.parent))
        system, _ = Xponge.load_molpsf(str(psf_path))
        Xponge.load_gro(str(gro_path), system)
        Xponge.save_sponge_input(system, str(output_dir / "poly"))
        """
    )
    result = subprocess.run(
        [sys.executable, "-c", child_code],
        cwd=XPONGE_REPO,
        check=False,
        text=True,
        capture_output=True,
    )
    if result.returncode != 0:
        pytest.skip(f"local Xponge reference failed: {result.stderr[-500:]}")
    return _snapshot_poly_outputs(output_dir)


def _run_xpongecpp_poly_workflow(output_dir: Path) -> dict[str, str]:
    child_code = textwrap.dedent(
        f"""
        from pathlib import Path

        import XpongeCPP as Xponge

        data_dir = Path({str(DATA_POLY_DIR)!r})
        output_dir = Path({str(output_dir)!r})

        def clean_gro_in_place(path: Path) -> None:
            lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
            if not lines:
                return
            last = lines[-1].strip().split()
            if len(last) == 9:
                lines[-1] = "   ".join(last[:3]) + "\\n"
                path.write_text("".join(lines), encoding="utf-8")

        gro_path = output_dir / "ps.gro"
        gro_path.write_text((data_dir / "ps.gro").read_text(encoding="utf-8"), encoding="utf-8")
        clean_gro_in_place(gro_path)
        itp_path = data_dir / "ps.itp"
        psf_path = data_dir / "ps.psf"

        Xponge.load_parameter_from_ffitp(itp_path.name, str(itp_path.parent))
        system, _ = Xponge.load_molpsf(str(psf_path))
        Xponge.load_gro(str(gro_path), system)
        Xponge.Save_SPONGE_Input(system, prefix="poly", dirname=str(output_dir))
        """
    )
    result = subprocess.run(
        [sys.executable, "-c", child_code],
        check=False,
        text=True,
        capture_output=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            "XpongeCPP poly workflow failed\nSTDOUT:\n"
            + result.stdout
            + "\nSTDERR:\n"
            + result.stderr
        )
    return _snapshot_poly_outputs(output_dir)


def test_poly_non_amber_reference_workflow_runs_with_original_xponge(tmp_path):
    actual = _run_original_poly_workflow(tmp_path)
    assert set(actual) == EXPECTED_OUTPUT_NAMES


def test_poly_non_amber_workflow_matches_original_xponge_reference(tmp_path):
    reference_dir = tmp_path / "reference"
    current_dir = tmp_path / "current"
    reference_dir.mkdir()
    current_dir.mkdir()

    reference = _run_original_poly_workflow(reference_dir)
    current = _run_xpongecpp_poly_workflow(current_dir)

    assert current == reference

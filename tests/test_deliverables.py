import json
from pathlib import Path
import zipfile


ROOT = Path(__file__).resolve().parents[1]


def test_submission_notebook_is_executed_without_errors():
    path = ROOT / "submission" / "BioRydberg_2_0_Pasqal_QoolQit_Contest.ipynb"
    notebook = json.loads(path.read_text(encoding="utf-8"))
    code_cells = [cell for cell in notebook["cells"] if cell["cell_type"] == "code"]
    assert code_cells
    assert all(cell["execution_count"] is not None for cell in code_cells)
    assert not [
        output
        for cell in code_cells
        for output in cell.get("outputs", [])
        if output.get("output_type") == "error"
    ]


def test_powerpoint_contains_exactly_two_slides():
    path = ROOT / "submission" / "BioRydberg_2_0_Pasqal_Contest.pptx"
    with zipfile.ZipFile(path) as archive:
        slides = [
            name
            for name in archive.namelist()
            if name.startswith("ppt/slides/slide") and name.endswith(".xml")
        ]
    assert len(slides) == 2


def test_required_results_exist():
    required = [
        ROOT / "results" / "study_summary.json",
        ROOT / "results" / "tables" / "study_metrics.csv",
        ROOT / "results" / "figures" / "05_baseline_vs_optimized.png",
        ROOT / "submission" / "BioRydberg_2_0_Pasqal_Contest.pdf",
    ]
    assert all(path.is_file() and path.stat().st_size > 0 for path in required)

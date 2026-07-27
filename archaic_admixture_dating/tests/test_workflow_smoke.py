from __future__ import annotations

import json

from archaic_admixture_dating.cli import main


def test_smoke_workflow_completes_end_to_end(tmp_path):
    output = tmp_path / "smoke"
    result = main(
        [
            "run-all",
            "--profile",
            "smoke",
            "--output",
            str(output),
            "--no-resume",
        ]
    )
    assert result == 0
    assert (output / "report" / "report.html").exists()
    state = json.loads((output / "checkpoints" / "run_all.json").read_text(encoding="utf-8"))
    assert state["state"] == "complete"
    report = (output / "report" / "report.md").read_text(encoding="utf-8")
    assert "Interpretation status: **inconclusive/data-limited**" in report
    assert "does **not** establish that Denisovans survived" in report
    assert "Direct late Denisovan admixture is not demonstrated" in report
